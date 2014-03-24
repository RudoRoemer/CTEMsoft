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
; CTEMsoft2013:EBSDshowMC.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDshowMC.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Display MC or MP using either modified Lambert or regular Lambert on a circle
;
;> @note Note that the circular Lambert projections are all precomputed at the time 
;> the data is read from the file...
;
;> @date 03/20/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDshowMC, dummy

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common Image_common, MCimage, MPimage
common EBSD_rawdata, accum_e, accum_z, MParray
common projections, mcxcircle, mcycircle, mpxcircle, mpycircle

wset,EBSDwidget_s.MCdrawID

energy = EBSDdata.mcenergymin + EBSDdata.Esel * EBSDdata.mcenergybinsize
if (EBSDdata.MCLSum eq 0) then begin
  image = reform(accum_e[EBSDdata.Esel,*,*])
  WIDGET_CONTROL, set_value=string(energy,format="(F5.2)"), EBSDwidget_s.MCenergyval
end 

if (EBSDdata.MCLSum eq 1) then begin
  image = reform(total(accum_e,1))
  WIDGET_CONTROL, set_value='sum', EBSDwidget_s.MCenergyval
end

if (EBSDdata.MCLSum eq 2) then begin   ; RGB display, needs some data preparation first
  numEbins = EBSDdata.mcenergynumbin
  energy = EBSDdata.mcenergymin + findgen(EBSDdata.mcenergynumbin) * EBSDdata.mcenergybinsize
  sz = size(accum_e,/dimensions)
  dims2 = sz[1]
  dims3 = sz[2]
  ener = fltarr(dims2,dims3)
  accum_t = float(total(accum_e,1))
  accum_t[where(accum_t eq 0.0)] = 1.0
  for i=0,dims2-1 do for j=0,dims3-1 do ener[i,j] = total(accum_e[0:*,i,j]*energy[0:*])/accum_t[i,j]

  bx = fltarr(dims2,dims3)
  by = fltarr(dims2,dims3)

; scale the energy as an angle between 0 and 270 counterclockwise
  minE = EBSDdata.mcenergymin
  maxE = EBSDdata.mcenergymax
  e = ener
  ener = (ener-minE)/(maxE-minE)
  ener = ener*180. - 180.
  ener *= !dtor
  c = cos(ener)
  s = -sin(ener)

  accum_t /= max(accum_t)

  bx = accum_t*c
  by = accum_t*s

  Core_colorwheel,bx,by,cimage,clegend

  image = cimage
  WIDGET_CONTROL, set_value='sumRGB', EBSDwidget_s.MCenergyval
end

if ((EBSDdata.MCLSmode eq 1) and (EBSDdata.MCLSum ne 2)) then begin
  image = bilinear(image,mcxcircle,mcycircle,missing = 0)
end 

if ((EBSDdata.MCLSmode eq 1) and (EBSDdata.MCLSum eq 2)) then begin
  for i=0,2 do begin
    image = reform(cimage[i,*,*])
    image = bilinear(image,mcxcircle,mcycircle,missing = 0)
    cimage[i,0:*,0:*] = image[0:*,0:*]
  endfor
end 

if (EBSDdata.MCLSum ne 2) then begin
  tvscl, image 
  WIDGET_CONTROL, set_value=string(min(image),format="(F9.1)"), EBSDwidget_s.MCmin
  WIDGET_CONTROL, set_value=string(max(image),format="(F9.1)"), EBSDwidget_s.MCmax
  MCimage = image
end else begin
  tvscl, cimage, true=1
; and add the color legend
  tvscl, clegend, dims2-25, dims3-50, true=1
  WIDGET_CONTROL, set_value='', EBSDwidget_s.MCmin
  WIDGET_CONTROL, set_value='', EBSDwidget_s.MCmax
  MCimage = cimage
  MCimage[0:2,dims2-25:*,dims3-50:*] = clegend[0:*,0:*,0:*]
end

; should we also display the Master Pattern ?
if (EBSDdata.MCMPboth eq 1) then begin
  wset,EBSDwidget_s.MPdrawID
  if (EBSDdata.MCLSum eq 0) then begin
    image = reform(MParray[*,*,EBSDdata.Esel])
  end 

  if (EBSDdata.MCLSum eq 1) then begin
; here, we want to display the energy-weighted master pattern as an example
    image = MParray[0:*,0:*,0] * congrid(reform(float(accum_e[0,0:*,0:*])),2*EBSDdata.MPimx+1,2*EBSDdata.MPimy+1)
    for i=1,EBSDdata.MCenergynumbin-1 do image += MParray[0:*,0:*,i] * congrid(reform(float(accum_e[i,0:*,0:*])),2*EBSDdata.MPimx+1,2*EBSDdata.MPimy+1)
  end

  if (EBSDdata.MCLSum eq 2) then begin   ; RGB display, for now left blank
    MPimage = 0
    erase
    empty
  end

  if ((EBSDdata.MCLSmode eq 1) and (EBSDdata.MCLSum ne 2)) then begin
    image = bilinear(image,mpxcircle,mpycircle,missing = 0)
  end 

; if ((EBSDdata.MCLSmode eq 1) and (EBSDdata.MCLSum eq 2)) then begin
; end 

  if (EBSDdata.MCLSum ne 2) then begin
    tvscl, image 
    WIDGET_CONTROL, set_value=string(min(image),format="(F9.1)"), EBSDwidget_s.MPmin
    WIDGET_CONTROL, set_value=string(max(image),format="(F9.1)"), EBSDwidget_s.MPmax
    MPimage = image
  end 
endif

end

