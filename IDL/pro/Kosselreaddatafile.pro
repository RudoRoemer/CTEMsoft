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
; CTEMsoft2013:Kosselreaddatafile.pro
;--------------------------------------------------------------------------
;
; PROGRAM: Kosselreaddatafile.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Reads the data file produced by the CTEMKossel.f90 program
;
;> @date 11/23/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro Kosselreaddatafile,dummy
;
;------------------------------------------------------------
; common blocks
common Kossel_widget_common, widget_s
common Kossel_data_common, data

; the next common block contains all the raw data needed to dsplay the ECpatterns
common Kossel_rawdata, rawdata
common PointGroups, PGTHD, PGTWD, DG



  dummy = 1

  openu,1,data.pathname+'/'+data.dataname,/f77
  progname = bytarr(14)
  readu,1,progname
  progname = strtrim(string(progname))

; if (progname ne 'CTEMKossel.f90') then begin
;   dummy = 0
;   print,''
;   WIDGET_CONTROL, SET_VALUE='The selected file is not a CTEMKossel output file ',widget_s.filename
;   goto,skipall
; end 

; version string
  scversion = bytarr(8)
  readu,1,scversion
  data.scversion = strtrim(string(scversion))

  WIDGET_CONTROL, SET_VALUE=string(float(data.filesize)/1024./1024.,FORMAT="(F8.2)")+' Mb', widget_s.filesize

  dims = lonarr(3)
  readu,1,dims
  data.datadims = long64(dims)
  data.imx = dims[0]
  data.imy = dims[1]
  WIDGET_CONTROL, SET_VALUE=string(data.imx,FORMAT="(I5)"), widget_s.imx
  WIDGET_CONTROL, SET_VALUE=string(data.imy,FORMAT="(I5)"), widget_s.imy
  WIDGET_CONTROL, SET_VALUE=string(data.datadims[2],FORMAT="(I4)"), widget_s.numthick

  xtalname = bytarr(132)
  readu,1,xtalname
  data.xtalname = strtrim(string(xtalname))
  WIDGET_CONTROL, SET_VALUE=data.xtalname, widget_s.xtalname

; accelerating voltage
  voltage = 0.0
  readu,1,voltage
  data.voltage = voltage
  data.wavelength= 1226.39/sqrt(voltage + 0.97845E-6 * voltage^2)
  WIDGET_CONTROL, SET_VALUE=string(data.voltage,FORMAT="(F9.1)"), widget_s.voltage

; beam convergence
  thetac = 0.0
  readu,1,thetac
  data.thetac = thetac
  WIDGET_CONTROL, SET_VALUE=string(data.thetac,FORMAT="(F7.3)"), widget_s.thetac

; ktmax
  ktmax = 0.0
  readu,1,ktmax
  data.ktmax = ktmax

; wave vector indices (3 longints)
  wavek = lonarr(3)
  readu,1,wavek
  data.wavek = wavek
  wv = '['+string(data.wavek[0],format="(I3)")+' '+ string(data.wavek[1],format="(I3)")+' '+ string(data.wavek[2],format="(I3)")+']'
  WIDGET_CONTROL, SET_VALUE=wv, widget_s.wavek

; foil normal indices (3 longints)
  fn = lonarr(3)
  readu,1,fn
  data.fn = fn

; number of wave vectors inside the disk
  numk = 0L
  readu,1,numk
  data.numk = numk
  WIDGET_CONTROL, SET_VALUE=string(data.numk,FORMAT="(I8)"), widget_s.numk

; dmin value (not editable or viewable)
  dmin = 0.0
  readu,1,dmin
  data.dmin = dmin

; first (ga) reflection (not viewable)
  ga = lonarr(3)
  readu,1,ga
  data.ga = ga
  wv = '('+string(data.ga[0],format="(I2)")+' '+ string(data.ga[1],format="(I2)")+' '+ string(data.ga[2],format="(I2)")+')'

; length of ga, used for proper scaling of the Laue center position
  galen = 0.0
  readu,1,galen
  data.galen = galen

; delta and gperp
  delta = 0.0
  readu,1,delta
  data.delta = delta
  gperp = fltarr(3)
  readu,1,gperp
  data.gperp = gperp

; various symmetry group numbers 
  symgroups = lonarr(8)
  readu,1,symgroups
  data.symgroups = symgroups
  widget_control, set_value=PGTHD[data.symgroups[0]], widget_s.symCPG
  widget_control, set_value=PGTWD[data.symgroups[5]], widget_s.symWPG

; starting thickness and thickness increment
  startthick = 0.0
  thickinc = 0.0
  readu,1,startthick,thickinc
  data.startthick = startthick
  data.thickinc = thickinc
; these are not shown in the main widget, but they are shown in a droplist widget in other areas

; and finally, read the data array
  rawdata = fltarr(data.datadims)
  readu,1,rawdata

  close,1

; this shouldn't happen, but in case there are points in any of the patterns
; that have unreasonable intensities, we can replace them by their local average...
  for i=0,data.datadims[2]-1 do begin
    slice = reform(rawdata[*,*,i])
    slice = median(slice,3)
    rawdata[0,0,i] = slice
  endfor

skipall:

end
