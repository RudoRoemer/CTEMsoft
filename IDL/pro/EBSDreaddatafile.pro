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
; CTEMsoft2013:EBSDreaddatafile.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDreaddatafile.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Reads the data files produced by the CTEMMC.f90 and CTEMEBSDmaster.f90 programs
;
;> @date 03/19/14 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro EBSDreaddatafile,MCFILE=MCFILE,MPFILE=MPFILE
;
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

; the next common block contains all the raw data needed to generate the EBSD patterns
common EBSD_rawdata, accum_e, accum_z, MParray, MParraysum


  Core_Print,' ',/blank
  EBSDdata.MCMPboth = 0

if keyword_set(MPFILE) then begin
  Core_Print,'Reading data file '+EBSDdata.mpfilename

  openu,1,EBSDdata.pathname+'/'+EBSDdata.mpfilename,/f77
; first a string of 132 characters
  progname = bytarr(132)
  readu,1,progname
  progname = strtrim(string(progname))
    Core_Print,' ->File generated by program '+progname+'<-'

; version string
  scversion = bytarr(8)
  readu,1,scversion
  EBSDdata.scversion = strtrim(string(scversion))
    Core_Print,'Version identifier : '+string(scversion) 

; display the file size in Mb 
  WIDGET_CONTROL, SET_VALUE=string(float(EBSDdata.mpfilesize)/1024./1024.,FORMAT="(F8.2)")+' Mb', EBSDwidget_s.mpfilesize

; structure file name
  xtalname = bytarr(132)
  readu,1,xtalname
  EBSDdata.xtalname = strtrim(string(xtalname))
    Core_Print,'Xtalname = ->'+EBSDdata.xtalname+'<-'
  WIDGET_CONTROL, SET_VALUE=EBSDdata.xtalname, EBSDwidget_s.xtalname

; energy file name
  energyname = bytarr(132)
  readu,1,energyname
  res = strtrim(string(energyname))
  finfo = file_info(res)
  EBSDdata.mcfilesize = finfo.size

  spos = strpos(res,'/',/reverse_search)
  dpos = strpos(res,'.',/reverse_search)
  plen = strlen(res)
  EBSDdata.mcpathname = strmid(res,0,spos)
  EBSDdata.mcfilename = strmid(res,spos+1)
    Core_Print,'MC filename = ->'+EBSDdata.mcfilename+'<-'
  WIDGET_CONTROL, SET_VALUE=EBSDdata.mcfilename, EBSDwidget_s.mcfilename

; npx, npy, numEbins, numset
  dims = lonarr(4)
  readu,1,dims
  EBSDdata.mpimx = dims[0]
  EBSDdata.mpimy = dims[1]
  EBSDdata.mcenergynumbin = dims[2]
  EBSDdata.numset= dims[3]
  EBSDdata.Asymsel = -1

  WIDGET_CONTROL, SET_VALUE=string(2*dims[0]+1,format="(I5)"), EBSDwidget_s.mpimx
  WIDGET_CONTROL, SET_VALUE=string(2*dims[1]+1,format="(I5)"), EBSDwidget_s.mpimy
  WIDGET_CONTROL, SET_VALUE=string(dims[2],format="(I5)"), EBSDwidget_s.mcenergynumbin
  
; energy levels
  EkeVs = fltarr(EBSDdata.mcenergynumbin)
  readu,1,EkeVs

; atomic numbers for asymmetric unit
  atnum = lonarr(EBSDdata.numset)
  readu,1,atnum
  EBSDdata.atnum(0:EBSDdata.numset-1) = atnum(0:EBSDdata.numset-1)

; Lambert projection type
  Ltype = bytarr(6)
  readu,1,Ltype
  Ltype = strtrim(string(Ltype))
  if (Ltype eq 'hexago') then EBSDdata.mpgridmode = ' Hexagonal' else EBSDdata.mpgridmode = ' Square'
  WIDGET_CONTROL, SET_VALUE=EBSDdata.mpgridmode, EBSDwidget_s.mpgridmode

; and finally the results array
  MParray = fltarr(2L*EBSDdata.mpimx+1L,2L*EBSDdata.mpimy+1L,EBSDdata.mcenergynumbin,EBSDdata.numset)
  readu,1,MParray
  if (EBSDdata.numset gt 1) then MParraysum = total(MParray,4) else MParraysum = MParray

  sz = size(MParray,/dimensions)
  if (EBSDdata.numset gt 1) then Core_Print,' Size of MParray data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)") +' x'+string(sz[2],format="(I5)") +' x'+string(sz[3],format="(I5)") else Core_Print,' Size of MParray data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)") +' x'+string(sz[2],format="(I5)")

; and close the file
  close,1

; and initialize the coordinate arrays for the Lambert transformation
  Core_LambertS2C,reform(MParray[*,*,0,0]),/mp
  Core_LambertS2SP,reform(MParray[*,*,0,0]),/mp

   WIDGET_CONTROL, EBSDwidget_s.MPbutton, sensitive=1
   WIDGET_CONTROL, EBSDwidget_s.detector, sensitive=1
  EBSDdata.MCMPboth = 1
endif






; read the Monte Carlo data file
if (keyword_set(MCFILE) or (EBSDdata.MCMPboth eq 1)) then begin
  Core_Print,'Reading data file '+EBSDdata.mcfilename
  EBSDdata.Esel = 0

  openu,1,EBSDdata.mcpathname+'/'+EBSDdata.mcfilename,/f77
; first a string of 132 characters
  progname = bytarr(132)
  readu,1,progname
  progname = strtrim(string(progname))
    Core_Print,' ->File generated by program '+progname+'<-'

; version string
  scversion = bytarr(8)
  readu,1,scversion
  EBSDdata.scversion = strtrim(string(scversion))
    Core_Print,'Version identifier : '+string(scversion) 

; display the file size in Mb 
  WIDGET_CONTROL, SET_VALUE=string(float(EBSDdata.mcfilesize)/1024./1024.,FORMAT="(F8.2)")+' Mb', EBSDwidget_s.mcfilesize

; structure file name
  xtalname = bytarr(132)
  readu,1,xtalname
  EBSDdata.xtalname = strtrim(string(xtalname))
    Core_Print,'Xtalname = ->'+EBSDdata.xtalname+'<-'
  WIDGET_CONTROL, SET_VALUE=EBSDdata.xtalname, EBSDwidget_s.xtalname

; six integers parameters (last one is not needed)
  dims = lonarr(6)
  readu,1,dims
  EBSDdata.mcenergynumbin = dims[0]
  EBSDdata.mcdepthnumbins = dims[1]
  EBSDdata.mcimx = (dims[2]-1L)/2L
  EBSDdata.mcimy = (dims[3]-1L)/2L
  EBSDdata.mctotale = dims[4]


  WIDGET_CONTROL, SET_VALUE=string(dims[0],format="(I5)"), EBSDwidget_s.mcenergynumbin
  WIDGET_CONTROL, SET_VALUE=string(dims[1],format="(I5)"), EBSDwidget_s.mcdepthnumbins
  WIDGET_CONTROL, SET_VALUE=string(dims[2],format="(I5)"), EBSDwidget_s.mcimx
  WIDGET_CONTROL, SET_VALUE=string(dims[3],format="(I5)"), EBSDwidget_s.mcimy
  WIDGET_CONTROL, SET_VALUE=string(dims[4],format="(I12)"), EBSDwidget_s.mctotale

; 5 more parameters, all doubles
  dims = dblarr(5)
  readu,1,dims
  EBSDdata.mcenergymax = dims[0]
  EBSDdata.mcenergymin = dims[1]
  EBSDdata.mcenergybinsize = dims[2]
  EBSDdata.mcdepthmax = dims[3]
  EBSDdata.mcdepthstep = dims[4]

  EBSDdata.voltage = EBSDdata.mcenergymax

  WIDGET_CONTROL, SET_VALUE=string(dims[0],format="(F7.2)"), EBSDwidget_s.mcenergymax
  WIDGET_CONTROL, SET_VALUE=string(dims[1],format="(F7.2)"), EBSDwidget_s.mcenergymin
  WIDGET_CONTROL, SET_VALUE=string(dims[2],format="(F7.2)"), EBSDwidget_s.mcenergybinsize
  WIDGET_CONTROL, SET_VALUE=string(dims[3],format="(F7.2)"), EBSDwidget_s.mcdepthmax
  WIDGET_CONTROL, SET_VALUE=string(dims[4],format="(F7.2)"), EBSDwidget_s.mcdepthstep
  WIDGET_CONTROL, SET_VALUE=string(dims[0],format="(F7.2)"), EBSDwidget_s.voltage

; sample tilt angles
  dims = dblarr(2)
  readu,1,dims
  EBSDdata.mcvangle = dims[0]
  EBSDdata.mchangle = dims[1]

  WIDGET_CONTROL, SET_VALUE=string(dims[0],format="(F7.2)"), EBSDwidget_s.mcvangle
  WIDGET_CONTROL, SET_VALUE=string(dims[1],format="(F7.2)"), EBSDwidget_s.mchangle

; Monte Carlo mode  'CSDA' or 'Discrete losses'
  mcm = bytarr(4)
  readu,1,mcm
  mcm = strtrim(string(mcm))
  if (mcm eq 'CSDA') then EBSDdata.mcmode = 'CSDA' else EBSDdata.mcmode = 'DLOS'
  WIDGET_CONTROL, SET_VALUE=EBSDdata.mcmode, EBSDwidget_s.mcmode


; and finally, we read the actual data arrays accum_e and accum_z
  accum_e = lonarr(EBSDdata.mcenergynumbin, 2*EBSDdata.mcimx+1,2*EBSDdata.mcimy+1)
  accum_z = lonarr(EBSDdata.mcenergynumbin, EBSDdata.mcdepthnumbins, 2*EBSDdata.mcimx/10+1,2*EBSDdata.mcimy/10+1)
  readu,1,accum_e
  readu,1,accum_z

; total number of BSE electrons
  EBSDdata.mcbse = total(accum_e)
  WIDGET_CONTROL, SET_VALUE=string(EBSDdata.mcbse,format="(I12)"), EBSDwidget_s.mcbse


  sz = size(accum_e,/dimensions)
    Core_Print,' Size of accum_e data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)")+' x'+string(sz[2],format="(I5)")
  sz = size(accum_z,/dimensions)
    Core_Print,' Size of accum_z data array : '+string(sz[0],format="(I5)")+' x'+string(sz[1],format="(I5)") +' x'+string(sz[2],format="(I5)") +' x'+string(sz[3],format="(I5)")

; and close the file
  close,1

; and initialize the coordinate arrays for the Lambert transformation
  Core_LambertS2C,reform(accum_e[0,*,*]),/mc
  Core_LambertS2SP,reform(accum_e[0,*,*]),/mc

; (de)activate buttons
   WIDGET_CONTROL, EBSDwidget_s.MCbutton, sensitive=1
   if (EBSDdata.MCMPboth eq 0) then begin
     WIDGET_CONTROL, EBSDwidget_s.MPbutton, sensitive=0
     WIDGET_CONTROL, EBSDwidget_s.detector, sensitive=0
     WIDGET_CONTROL, SET_VALUE=' ', EBSDwidget_s.mpfilename
     WIDGET_CONTROL, SET_VALUE=' ', EBSDwidget_s.mpfilesize
     WIDGET_CONTROL, SET_VALUE=' ', EBSDwidget_s.mpimx
     WIDGET_CONTROL, SET_VALUE=' ', EBSDwidget_s.mpimy
     WIDGET_CONTROL, SET_VALUE=' ', EBSDwidget_s.mpgridmode
   endif
end


  Core_Print,'Completed reading data file(s)',/blank


skipall:

end
