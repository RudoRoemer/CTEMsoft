;
; Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
; All rights reserved.
;
; Redistribution and use in source and binary forms, with or without modification, are 
; permitted provided that the following conditions are met:
;
;     - Redistributions of source code must retain the above copyright notice, this list 
;        of conditions and the following disclaimer.
;     - Redistributions in binary form must reproduce the above copyright notice, this 
;        list of conditions and the following disclaimer in the documentation and/or 
;        other materials provided with the distribution.
;     - Neither the names of Marc De Graef, Carnegie Mellon University nor the names 
;        of its contributors may be used to endorse or promote products derived from 
;        this software without specific prior written permission.
;
; THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
; AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE 
; IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE 
; ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE 
; LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
; DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR 
; SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER 
; CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, 
; OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE 
; USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
; ###################################################################
;--------------------------------------------------------------------------
; CTEMsoft2013:EBSDreadHDFdatafile.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDreadHDFdatafile.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Reads the HDF data files produced by the EMMC.f90 and EMEBSDmaster.f90 programs
;
;> @date 03/19/14 MDG 1.0 first attempt 
;> @date 04/02/15 MDG 2.0 modfied from EBSDreaddatafile to cover HDF formatted files
;> @date 05/07/15 MDG 2.1 modified to accommodate changes in path names for the entire EMsoft package
;> @date 10/29/15 MDG 2.2 added support for ECP master files
;> @date 10/31/15 MDG 2.3 removed widget fields and redirected all output to status widget
;--------------------------------------------------------------------------
pro EBSDreadHDFdatafile,MCFILE=MCFILE,MPFILE=MPFILE
;
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

; the next common block contains all the raw data needed to generate the EBSD patterns
common EBSD_rawdata, accum_e, accum_z, mLPNH, mLPSH

  Core_Print,' ',/blank
  EBSDdata.MCMPboth = 0


EMdatapathname = ''

if keyword_set(MPFILE) then begin
  Core_Print,'Reading data file '+EBSDdata.mpfilename

EMdatapathname = Core_getenv(/data)

; first make sure that this is indeed an HDF file
  res = H5F_IS_HDF5(EBSDdata.pathname+'/'+EBSDdata.mpfilename)
  if (res eq 0) then begin
    Core_Print,'  This is not an HDF file ! ',/blank
    goto,skipall
  endif

; ok, so it is an HDF file; let's open it
  file_id = H5F_OPEN(EBSDdata.pathname+'/'+EBSDdata.mpfilename)
  if (file_id eq -1L) then begin
    Core_Print,'  Error opening file',/blank
    goto, skipall
  endif 

  
; open the EMheader group
  group_id = H5G_open(file_id,'EMheader')

;  open and read the ProgramName dataset
  dset_id = H5D_open(group_id,'ProgramName')
  z = H5D_read(dset_id) 
  progname = strtrim(z[0],2)
  H5D_close,dset_id
    Core_Print,' ->File generated by program '+progname+'<-'

  if (progname eq 'EMEBSDmaster.f90') then EBSDdata.EBSDorECP = 0 else EBSDdata.EBSDorECP=1

; open and read the Version dataset
  dset_id = H5D_open(group_id,'Version')
  z = H5D_read(dset_id) 
  scversion = strtrim(z[0],2)
  H5D_close,dset_id
  EBSDdata.scversion = strtrim(scversion,2)
    Core_Print,'      Version identifier : '+scversion 

; close the EMheader group
  H5G_close,group_id

; display the file size in Mb 
  Core_print,'      MP data file size : '+string(float(EBSDdata.mpfilesize)/1024./1024.,FORMAT="(F8.2)")+' Mb'

; open the NMLparameters/EBSDMasterNameList group
  if (EBSDdata.EBSDorECP eq 0) then begin
    group_id = H5G_open(file_id,'NMLparameters/EBSDMasterNameList')
  end else begin
    group_id = H5G_open(file_id,'NMLparameters/ECPMasterNameList')
  endelse

; open and read the energyfile dataset
  dset_id = H5D_open(group_id,'energyfile')
  z = H5D_read(dset_id) 
  res = strtrim(z[0],2)
  H5D_close,dset_id
  finfo = file_info(res)
  EBSDdata.mcfilesize = finfo.size

  spos = strpos(res,'/',/reverse_search)
  dpos = strpos(res,'.',/reverse_search)
  plen = strlen(res)
  EBSDdata.mcpathname = strmid(res,0,spos)
  EBSDdata.mcfilename = strmid(res,spos+1)
  WIDGET_CONTROL, SET_VALUE=EBSDdata.mcfilename, EBSDwidget_s.mcfilename

; outname (to reset the current filename to the correct EMdatapathname)
  dset_id = H5D_open(group_id,'outname')
  z = H5D_read(dset_id)
  res = strtrim(z[0],2)
  H5D_close,dset_id
  EBSDdata.mpfilename = strmid(res,0,strlen(res))

; npx, npy
  dset_id = H5D_open(group_id,'npx')
  z = H5D_read(dset_id) 
  res = long(z[0])
  H5D_close,dset_id
  EBSDdata.mpimx = res
  EBSDdata.mpimy = res
  EBSDdata.Asymsel = -1

; close the group
  H5g_close,group_id

; open the EMData group
  group_id = H5G_open(file_id,'EMData')
 
; numEbins, numset datasets
  if (EBSDdata.EBSDorECP eq 0) then begin
    dset_id = H5D_open(group_id,'numEbins')
    z = H5D_read(dset_id) 
    res = long(z[0])
    H5D_close,dset_id
    EBSDdata.mcenergynumbin = res
  end else begin
    EBSDdata.mcenergynumbin = 1
  endelse

  dset_id = H5D_open(group_id,'numset')
  z = H5D_read(dset_id) 
  res = long(z[0])
  H5D_close,dset_id
  EBSDdata.numset= res

  Core_print,'      Number of energy bins : '+string(EBSDdata.mcenergynumbin,format="(I5)")
  
; energy levels (currently not used anywhere)
  if (EBSDdata.EBSDorECP eq 0) then  begin
    dset_id = H5D_open(group_id,'EkeVs')
    EkeVs = H5D_read(dset_id) 
    H5D_close,dset_id
  end else begin
    dset_id = H5D_open(group_id,'EkeV')
    EkeVs = H5D_read(dset_id) 
    H5D_close,dset_id
  endelse

; atomic numbers for asymmetric unit
  dset_id = H5D_open(group_id,'cell%ATOM_type')
  atnum = H5D_read(dset_id) 
  H5D_close,dset_id
  EBSDdata.atnum(0:EBSDdata.numset-1) = atnum(0:EBSDdata.numset-1)

; Lambert projection type (to be removed in a later version because it is always square)
  dset_id = H5D_open(group_id,'squhex')
  z = H5D_read(dset_id) 
  res = strtrim(z[0],2)
  H5D_close,dset_id
  if (res eq 'hexago') then EBSDdata.mpgridmode = ' Hexagonal' else EBSDdata.mpgridmode = ' Square'
  Core_print,'      Lambert array grid mode : '+EBSDdata.mpgridmode

; and finally the results arrays 
  dset_id = H5D_open(group_id,'mLPNH')
  mLPNH = H5D_read(dset_id) 
  H5D_close,dset_id
  dset_id = H5D_open(group_id,'mLPSH')
  mLPSH = H5D_read(dset_id) 
  H5D_close,dset_id

; close the group
  H5G_close,group_id

; resize the mLPNH/mLPSH arrays to the correct number of dimensions
  if (EBSDdata.EBSDorECP eq 0) then begin
    sz = size(mLPNH)
    if (sz[0] eq 3) then begin
      mLPNH = reform(mLPNH,sz[1],sz[2],sz[3],1)
      mLPSH = reform(mLPSH,sz[1],sz[2],sz[3],1)
    endif
    sz = size(mLPNH,/dimensions)
    Core_Print,'      Size of mLPNH data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)") +' x'+string(sz[2],format="(I5)") +' x'+string(sz[3],format="(I5)") 
  end else begin
    sz = size(mLPNH)
    if (sz[0] eq 2) then begin
      mLPNH = reform(mLPNH,sz[1],sz[2],1)
      mLPSH = reform(mLPSH,sz[1],sz[2],1)
    endif
    sz = size(mLPNH,/dimensions)
    Core_Print,'      Size of mLPNH data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)") +' x'+string(sz[2],format="(I5)") 
  endelse

; and close the file
  H5F_close,file_id

; and initialize the coordinate arrays for the Lambert transformation
  if (EBSDdata.EBSDorECP eq 0) then begin
    Core_LambertS2C,reform(mLPNH[*,*,0,0]),/mp
    Core_LambertS2SP,reform(mLPNH[*,*,0,0]),/mp
  end else begin
    Core_LambertS2C,reform(mLPNH[*,*,0]),/mp
    Core_LambertS2SP,reform(mLPNH[*,*,0]),/mp
  endelse

  WIDGET_CONTROL, EBSDwidget_s.MPbutton, sensitive=1
  WIDGET_CONTROL, EBSDwidget_s.detector, sensitive=1
  EBSDdata.MCMPboth = 1
  Core_print,' '
endif


; read the Monte Carlo data file
if (keyword_set(MCFILE) or (EBSDdata.MCMPboth eq 1)) then begin
  Core_Print,'Reading data file '+EBSDdata.mcfilename
  EBSDdata.Esel = 0

  if (EBSDdata.MCMPboth eq 1) then begin ; get the file size of the MC file
    fname = EBSDdata.pathname+'/'+EBSDdata.mcfilename
    finfo = file_info(fname)
    EBSDdata.mcfilesize = finfo.size
  endif

; first make sure that this is indeed an HDF file
  res = H5F_IS_HDF5(EMdatapathname+EBSDdata.mcpathname+'/'+EBSDdata.mcfilename)
  if (res eq 0) then begin
    Core_Print,'  This is not an HDF file ! ',/blank
    goto,skipall
  endif

; ok, so it is an HDF file; let's open it
  file_id = H5F_OPEN(EMdatapathname+EBSDdata.mcpathname+'/'+EBSDdata.mcfilename)
  if (file_id eq -1L) then begin
    Core_Print,'  Error opening file',/blank
    goto, skipall
  endif 

; open the EMheader group
  group_id = H5G_open(file_id,'EMheader')

;  open and read the ProgramName dataset
  dset_id = H5D_open(group_id,'ProgramName')
  z = H5D_read(dset_id) 
  progname = strtrim(z[0],2)
  H5D_close,dset_id
    Core_Print,' ->File generated by program '+progname+'<-'

; open and read the Version dataset
  dset_id = H5D_open(group_id,'Version')
  z = H5D_read(dset_id) 
  scversion = strtrim(z[0],2)
  H5D_close,dset_id
  EBSDdata.scversion = strtrim(scversion,2)
    Core_Print,'      Version identifier : '+scversion 

; close the EMheader group
  H5G_close,group_id

; display the file size in Mb 
  Core_print,'      MC Data file size : '+string(float(EBSDdata.mcfilesize)/1024./1024.,FORMAT="(F8.2)")+' Mb'

; structure file name
  group_id = H5G_open(file_id,'NMLparameters/MCCLNameList')

; open and read the xtalname dataset
  dset_id = H5D_open(group_id,'xtalname')
  z = H5D_read(dset_id) 
  EBSDdata.xtalname = strtrim(z[0],2)
  H5D_close,dset_id
  Core_print,'      Structure file name : '+EBSDdata.xtalname

; open and read the mode dataset, to distinguish between EBSD and ECP data
  dset_id = H5D_open(group_id,'mode')
  m = H5D_read(dset_id)
  if (m[0] eq 'full') then EBSDdata.EBSDorECP = 0 else EBSDdata.EBSDorECP = 1
  H5D_close,dset_id

; open and read the numsx dataset 
  dset_id = H5D_open(group_id,'numsx')
  EBSDdata.mcimx = (long(H5D_read(dset_id))-1L)/2L
  EBSDdata.mcimy = EBSDdata.mcimx
  H5D_close,dset_id

; open and read the totnum_el dataset
  dset_id = H5D_open(group_id,'totnum_el')
  EBSDdata.mctotale = long(H5D_read(dset_id))
  H5D_close,dset_id

; open and read the EkeV dataset
  dset_id = H5D_open(group_id,'EkeV')
  EBSDdata.mcenergymax = double(H5D_read(dset_id))
  H5D_close,dset_id
  EBSDdata.voltage = EBSDdata.mcenergymax

  if (EBSDdata.EBSDorECP eq 0) then begin
; open and read the Ehistmin dataset
    dset_id = H5D_open(group_id,'Ehistmin')
    EBSDdata.mcenergymin = double(H5D_read(dset_id))
    H5D_close,dset_id

; open and read the Ebinsize dataset
    dset_id = H5D_open(group_id,'Ebinsize')
    EBSDdata.mcenergybinsize = double(H5D_read(dset_id))
    H5D_close,dset_id
  end else begin
    EBSDdata.mcenergymin = EBSDdata.mcenergymax
    EBSDdata.mcenergybinsize = 0.D0
  endelse

; open and read the depthmax dataset
  dset_id = H5D_open(group_id,'depthmax')
  EBSDdata.mcdepthmax = double(H5D_read(dset_id))
  H5D_close,dset_id

; open and read the depthstep dataset
  dset_id = H5D_open(group_id,'depthstep')
  EBSDdata.mcdepthstep = double(H5D_read(dset_id))
  H5D_close,dset_id

; open and read the sig dataset
  if (EBSDdata.EBSDorECP eq 0) then begin
    dset_id = H5D_open(group_id,'sig')
    EBSDdata.mcvangle = double(H5D_read(dset_id))
    H5D_close,dset_id
  end else begin
    dset_id = H5D_open(group_id,'sigstart')
    EBSDdata.mcsigstart = double(H5D_read(dset_id))
    H5D_close,dset_id

    dset_id = H5D_open(group_id,'sigend')
    EBSDdata.mcsigend = double(H5D_read(dset_id))
    H5D_close,dset_id

    dset_id = H5D_open(group_id,'sigstep')
    EBSDdata.mcsigstep = double(H5D_read(dset_id))
    H5D_close,dset_id
  endelse

; open and read the omega dataset
  dset_id = H5D_open(group_id,'omega')
  EBSDdata.mchangle = double(H5D_read(dset_id))
  H5D_close,dset_id

; open and read the MCmode dataset
  dset_id = H5D_open(group_id,'MCmode')
  res = strtrim(H5D_read(dset_id))
  H5D_close,dset_id
  if (res eq 'CSDA') then EBSDdata.mcmode = 'CSDA' else EBSDdata.mcmode = 'DLOS'

  Core_print,'      Lambert dimensions : '+string(2*EBSDdata.mcimx+1,format="(I5)")+' by '+string(2*EBSDdata.mcimy+1,format="(I5)")
  Core_print,'      Incident beam voltage [kV] '+string(EBSDdata.voltage,format="(F7.2)")
  Core_print,'      Monte Carlo mode : '+EBSDdata.mcmode

; close the group
  H5G_close,group_id

; and open the EMData group
  group_id = H5G_open(file_id,'EMData')

; open and read the numEbins dataset
  if (EBSDdata.EBSDorECP eq 0) then begin
    dset_id = H5D_open(group_id,'numEbins')
    EBSDdata.mcenergynumbin= long(H5D_read(dset_id))
    H5D_close,dset_id
    Core_print,'      Sample tilt angles omega/sigma [degrees] : '+string(EBSDdata.mchangle,format="(F7.2)")+'/'+string(EBSDdata.mcvangle,format="(F7.2)")
    Core_print,'      Min/Max energy [keV] : '+string(EBSDdata.mcenergymin,format="(F7.2)")+'/'+string(EBSDdata.mcenergymin,format="(F7.2)")
    Core_print,'      Energy binsize [keV] : '+string(EBSDdata.mcenergybinsize,format="(F7.2)")
    Core_print,'      Number of energy bins : '+string(EBSDdata.mcenergynumbin,format="(I5)")
  end else begin
    dset_id = H5D_open(group_id,'numangle')
    EBSDdata.mcenergynumbin= long(H5D_read(dset_id))
    H5D_close,dset_id
    Core_print,'      Electron energy [keV] : '+string(EBSDdata.mcenergymax,format="(F7.2)")
    Core_print,'      Number of beam tilt angles : '+string(EBSDdata.mcenergynumbin,format="(I5)")
    Core_print,'      Beam tilt range start/end/step [degrees] : '+string(EBSDdata.mcsigstart,format="(F7.2)")+'/'+string(EBSDdata.mcsigend,format="(F7.2)")+'/'+string(EBSDdata.mcsigstep,format="(F7.2)")
  endelse

; open and read the numzbins dataset
  dset_id = H5D_open(group_id,'numzbins')
  EBSDdata.mcdepthnumbins = long(H5D_read(dset_id))
  H5D_close,dset_id

  Core_print,'      Maximum integration depth [nm] : '+string(EBSDdata.mcdepthmax,format="(F7.2)")
  Core_print,'      Integration depth step size [nm] : '+string(EBSDdata.mcdepthstep,format="(F7.2)")
  Core_print,'      Number of depth bins : '+string(EBSDdata.mcdepthnumbins,format="(I5)")
  

; and finally, we read the actual data arrays accum_e and accum_z
  dset_id = H5D_open(group_id,'accum_e')
  accum_e = long(H5D_read(dset_id))
  H5D_close,dset_id

; we do not need accum_z at all in this program, so we won't read it...
;  dset_id = H5D_open(group_id,'accum_z')
;  accum_z = long(H5D_read(dset_id))
;  H5D_close,dset_id

; close the group
  H5G_close,group_id

; total number of BSE electrons
  EBSDdata.mcbse = total(accum_e)
  Core_print,'      Total number of incident electrons : '+string(EBSDdata.mctotale,format="(I12)")
  Core_print,'      Number of BSEs : '+string(EBSDdata.mcbse,format="(I12)")

  sz = size(accum_e,/dimensions)
    Core_Print,'      Size of accum_e data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)")+' x'+string(sz[2],format="(I5)")
;  sz = size(accum_z,/dimensions)
;    Core_Print,'      Size of accum_z data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)") +' x'+string(sz[2],format="(I5)") +' x'+string(sz[3],format="(I5)")

; and close the file
  H5F_close,file_id

; and initialize the coordinate arrays for the Lambert transformation
  Core_LambertS2C,reform(accum_e[0,*,*]),/mc
  Core_LambertS2SP,reform(accum_e[0,*,*]),/mc

; (de)activate buttons
   WIDGET_CONTROL, EBSDwidget_s.MCbutton, sensitive=1
   if (EBSDdata.MCMPboth eq 0) then begin
     WIDGET_CONTROL, EBSDwidget_s.MPbutton, sensitive=0
     WIDGET_CONTROL, EBSDwidget_s.detector, sensitive=0
     WIDGET_CONTROL, SET_VALUE=' ', EBSDwidget_s.mpfilename
   endif
end


  Core_Print,'Completed reading data file(s)',/blank


skipall:

end
