;
; Copyright (c) 2013, Marc De Graef/Carnegie Mellon University
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
; CTEMsoft2013:CBEDreaddatafile.pro
;--------------------------------------------------------------------------
;
; PROGRAM: CBEDreaddatafile.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Reads the data file produced by the CTEMlacbed.f90 or CTEMmbcbed.f90 programs
;
;> @date 10/08/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro CBEDreaddatafile,LACBED=LACBED,MBCBED=MBCBED
;
;------------------------------------------------------------
; common blocks
common CBED_widget_common, widget_s
common CBED_data_common, data
common PointGroups, PGTHD, PGTWD, DG

; the next common block contains all the raw data needed to generate the CBED patterns
common CBED_rawdata, gvecs, gmult, gtt, gxy, disks, numHOLZ, HOLZlist
common CBED_current, BFcurrent, DFcurrent, RGBcurrent, mask


  CBEDprint,'Reading data file '+data.dataname

  openu,1,data.pathname+'/'+data.dataname,/f77
; first a pair of strings of 132 characters each
  progname = bytarr(14)
  readu,1,progname
  progname = strtrim(string(progname))
    CBEDprint,' ->File generated by program '+progname+'<-'

  if ( keyword_set(LACBED) and (progname eq 'CTEMmbcbed.f90') ) then begin
    CBEDprint,' This file was generated by the CTEMmbcbed.f90 program; it is not a LACBED file ',/blank
    goto,skipall
  end

  if ( keyword_set(MBCBED) and (progname eq 'CTEMlacbed.f90') ) then begin
    CBEDprint,' This file was generated by the CTEMlacbed.f90 program; it is not an MBCBED file ',/blank
    goto,skipall
  end

; ok, let's read the file; the first part is common to both file formats

; version string
  scversion = bytarr(8)
  readu,1,scversion
  data.scversion = strtrim(string(scversion))
    CBEDprint,'Version identifier : '+string(scversion) 

; finfo = file_info(data.dataname)
; data.filesize = finfo.size
  WIDGET_CONTROL, SET_VALUE=string(data.filesize,FORMAT="(I14)")+' bytes', widget_s.filesize

  dims = lonarr(4)
  readu,1,dims
  data.datadims = long64(dims)
  data.imx = dims[0]
  data.imy = dims[1]
  data.numt = dims[2]
  data.numfam = dims[3]
    CBEDprint,' data dimensions : '+string(dims[0],"(I5)")+ string(dims[1],"(I5)")+ string(dims[2],"(I5)")+ string(dims[3],"(I5)")
  WIDGET_CONTROL, SET_VALUE=string(data.imx,FORMAT="(I5)"), widget_s.imx
  WIDGET_CONTROL, SET_VALUE=string(data.imy,FORMAT="(I5)"), widget_s.imy
; WIDGET_CONTROL, SET_VALUE=string(data.numt,FORMAT="(I14)"), widget_s.numt
  WIDGET_CONTROL, SET_VALUE=string(data.numfam,FORMAT="(I5)"), widget_s.numfam

; create the mask
  mask = shift(dist(data.datadims[0]),data.datadims[0]/2,data.datadims[0]/2)
  mask[where (mask le data.datadims[0]/2)] = 1.0
  mask[where (mask gt data.datadims[0]/2)] = 0.0
  mask = mask gt 0.5

  xtalname = bytarr(132)
  readu,1,xtalname
  data.xtalname = strtrim(string(xtalname))
    CBEDprint,'Xtalname = ->'+data.xtalname+'<-'
  WIDGET_CONTROL, SET_VALUE=data.xtalname, widget_s.xtalname

; accelerating voltage
  voltage = 0.0
  readu,1,voltage
  data.wavelength= 1226.39/sqrt(voltage + 0.97845E-6 * voltage^2)
    CBEDprint,'Wave length = '+string(data.wavelength,FORMAT="(F7.4)")
  WIDGET_CONTROL, SET_VALUE=string(data.wavelength,FORMAT="(F7.4)"), widget_s.wavelength

; beam convergence
  thetac = 0.0
  readu,1,thetac
  data.thetac = thetac
    CBEDprint,'Beam convergence = '+string(data.thetac,FORMAT="(F6.3)")
  WIDGET_CONTROL, SET_VALUE=string(data.thetac,FORMAT="(F6.3)"), widget_s.thetac

; wave vector indices (3 longints)
  wavek = lonarr(3)
  readu,1,wavek
  data.wavek = wavek
  wv = '['+string(data.wavek[0],format="(I2)")+' '+ string(data.wavek[1],format="(I2)")+' '+ string(data.wavek[2],format="(I2)")+']'
    CBEDprint,'Wave vector = '+wv
  WIDGET_CONTROL, SET_VALUE=wv, widget_s.wavek

; number of wave vectors inside the disk
  numk = 0L
  readu,1,numk
  data.numk = numk
    CBEDprint,'Number of k-vectors in disk = '+string(data.numk,FORMAT="(I)")
  WIDGET_CONTROL, SET_VALUE=string(data.numk,FORMAT="(I8)"), widget_s.numk

; first (ga) reflection
  ga = lonarr(3)
  readu,1,ga
  data.ga = ga
  wv = '('+string(data.ga[0],format="(I2)")+' '+ string(data.ga[1],format="(I2)")+' '+ string(data.ga[2],format="(I2)")+')'
    CBEDprint,'Horizontal g-vector = '+wv
  WIDGET_CONTROL, SET_VALUE=wv, widget_s.ga

; maximum HOLZ layer number 
  maxHOLZ = 0L
  readu,1,maxHOLZ
  data.maxHOLZ = maxHOLZ
    CBEDprint,'Maximum HOLZ layer number = '+string(data.maxHOLZ,FORMAT="(I3)")
  WIDGET_CONTROL, SET_VALUE=string(data.maxHOLZ,FORMAT="(I5)"), widget_s.maxHOLZ

; various symmetry group numbers 
  symgroups = lonarr(8)
  readu,1,symgroups
  data.symgroups = symgroups
    CBEDprint,' Crystallographic point group = '+PGTHD[symgroups[0]]
    CBEDprint,' Laue PG                      = '+PGTHD[symgroups[1]]
    CBEDprint,' Diffraction PG               = '+DG[symgroups[2]]
    CBEDprint,' Projection Diff. PG          = '+DG[symgroups[3]]
    CBEDprint,' Bright Field PG              = '+PGTWD[symgroups[4]]
    CBEDprint,' Whole Pattern PG             = '+PGTWD[symgroups[5]]
    CBEDprint,' Dark Field General PG        = '+PGTWD[symgroups[6]]
    CBEDprint,' Dark Field Special PG        = '+PGTWD[symgroups[7]]
  widget_control, set_value=PGTHD[data.symgroups[0]], widget_s.symCPG
  widget_control, set_value=PGTHD[data.symgroups[1]], widget_s.symLPG
  widget_control, set_value=DG[data.symgroups[2]], widget_s.symDPG
  widget_control, set_value=DG[data.symgroups[3]], widget_s.symPDG
  widget_control, set_value=PGTWD[data.symgroups[4]], widget_s.symBFG
  widget_control, set_value=PGTWD[data.symgroups[5]], widget_s.symWPG
  widget_control, set_value=PGTWD[data.symgroups[6]], widget_s.symDFG
  widget_control, set_value=PGTWD[data.symgroups[7]], widget_s.symDFS

; initialize the Whole Pattern symmetry
  CBEDGenerate2DSymmetry,data.symgroups[5]


; starting thickness and thickness increment
  startthick = 0.0
  thickinc = 0.0
  readu,1,startthick,thickinc
  data.startthick = startthick
  data.thickinc = thickinc
    cbedprint,' Starting thickness [nm]   = '+string(startthick,FORMAT="(F6.3)")
    cbedprint,' Thickness increment [nm]  = '+string(thickinc,FORMAT="(F6.3)")
; these are not shown in the main widget, but they are shown in a droplist widget in other areas


; ok, the above items are common to the two data file formats
; here we distinguish between them ... 
  if ( keyword_set(LACBED) or (progname eq 'CTEMlacbed.f90') ) then begin
    CBEDprogressbar,0.0
      CBEDprint,'Allocating memory for diffraction data array',/blank
    disks = fltarr(dims)
    gvecs = lonarr(3,dims[3])
    gmult = lonarr(dims[3])
    gtt = fltarr(dims[3])
    gxy = fltarr(2,dims[3])
    hkl = lonarr(3)
    m = 0L
    g = 0.0
    x = 0.0
    y = 0.0
    dd = fltarr(dims[0],dims[1],dims[2])
    for i=0,dims[3]-1 do begin
      readu,1,hkl,m,g,x,y
      gvecs[0,i] = hkl
      gmult[i] = m
      gtt[i] = g
      gxy[0:1,i] = [x,y]
      readu,1,dd 
      disks[0,0,0,i] = dd
      CBEDprogressbar,100.0*float(i)/float(data.datadims[3])
    endfor
    CBEDprogressbar,100.0
    close,1
; activate the LACBED and CBED buttons
    WIDGET_CONTROL, widget_s.startLACBED, sensitive=1
    WIDGET_CONTROL, widget_s.startCBED, sensitive=1
    WIDGET_CONTROL, widget_s.startMBCBED, sensitive=0
  end; else begin ; this is the MBCBED data format
;end

  CBEDprint,'Completed reading data file',/blank

skipall:

end
