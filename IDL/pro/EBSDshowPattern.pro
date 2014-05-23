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
; CTEMsoft2013:EBSDshowPattern.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDshowPattern.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main routine for display of EBSD patterns
;
;> @date 05/22/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDshowPattern, single=single

; the keyword /single indicates that only one pattern is available 

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common EBSDpatterns, pattern, image, finalpattern

wset,EBSDwidget_s.PatternDrawID
erase
empty

; set the min and max fields
EBSDdata.Patternmin = min(pattern)
EBSDdata.Patternmax = max(pattern)
print,EBSDdata.Patternmin,EBSDdata.Patternmax

WIDGET_CONTROL, set_value=string(EBSDdata.Patternmin,format="(F7.2)"), EBSDwidget_s.Patternmin
WIDGET_CONTROL, set_value=string(EBSDdata.Patternmax,format="(F7.2)"), EBSDwidget_s.Patternmax

; display the pattern
if keyword_set(single) then begin
; first apply the necessary intensity scaling to the current pattern

; what kind of intensity scaling do we need?
  if (EBSDdata.PatternScaling eq 0) then begin  ; this is regular linear scaling
    finalpattern = bytscl(pattern,min=EBSDdata.Patternmin,max=EBSDdata.Patternmax)
  end else begin
    image = pattern^EBSDdata.gammavalue
    finalpattern = bytscl(image)
  end

; then we apply the pattern origin (not sure if this is a useful operation or not)



; finally, we use the binning factor
  if (EBSDdata.detbinning ne 0) then finalpattern = congrid(finalpattern,EBSDdata.detnumsx/2^EBSDdata.detbinning,EBSDdata.detnumsy/2^EBSDdata.detbinning)

; and we display the result
  tv,finalpattern

end else begin
; this is a pre-computed series of patterns, so we can do direct display

end


end
