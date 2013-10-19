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
; CTEMsoft2013:CBEDcircles.pro
;--------------------------------------------------------------------------
;
; PROGRAM: CBEDcircles.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Draw diffraction and Laue limiting circles
;
;> @date 10/15/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro CBEDcircles,dummy

;------------------------------------------------------------
; common blocks
common CBED_widget_common, widget_s
common CBED_data_common, data
common CBED_rawdata, gvecs, gmult, gtt, gxy, disks, numHOLZ, HOLZlist
common CBED_HOLZlists, HOLZvals
common fontstrings, fontstr, fontstrlarge, fontstrsmall
common SYM2D, SYM_MATnum, SYM_direc
common CBEDcirclestuff, CBEDschematic, midx, midy, drad, cang, scl

; the LACBED convergence angle is shown as a red circle of radius 180 pixels
; inside the fixed size display window (401x401)

; the output will be a schematic RGB image
CBEDschematic = bytarr(3,data.detwinx,data.detwiny)

midx = (data.detwinx-1)/2
midy = (data.detwiny-1)/2

; use the Z-buffer to draw the individual components
set_plot,'Z'
device,set_resolution=[data.detwinx,data.detwiny]
erase

; first get the main LACBED convergence circle
drad = 180
numth = 361
th = findgen(numth) * !dtor
ct = cos(th)
st = sin(th)

plots,midx+drad*ct,midy+drad*st,/dev,thick=2

green = tvrd()
erase

; then draw the limiting circle for the Laue center positions
; this depends on the user selected convergence angle
cang =  drad*(data.thetau/data.thetac)
frac = drad - cang
plots,midx+frac*ct,midy+frac*st,/dev,thick=2
red = tvrd()
erase

; finally, draw all the ZOLZ reflections to scale in blue
; (at least, those that fall inside the window)
glen = sqrt(gxy[0,1]^2+gxy[1,1]^2)
scl = drad * (gtt[1] / data.thetac) /glen

for i=0,numHOLZ[0]-1 do begin 
  posxy = scl * CBEDApply2DSymmetryPoint(reform(gxy[0:1,i]))
  for j=0,SYM_MATnum-1 do begin
    px = midx + posxy[0,j]
    py = midy + posxy[1,j]
    plots,px + cang*ct, py + cang*st,/dev,thick=2
  endfor
end

blue = tvrd()


CBEDschematic[0,0:*,0:*] = red
CBEDschematic[1,0:*,0:*] = green
CBEDschematic[2,0:*,0:*] = blue

set_plot,'X'

wset,data.diskdrawID
erase
tvscl,CBEDschematic,true=1




end
