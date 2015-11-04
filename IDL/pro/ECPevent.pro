;
; Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
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
; EMsoft:ECPevent.pro
;--------------------------------------------------------------------------
;
; PROGRAM: ECPevent.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 10/30/15 MDG 1.0 first version
;--------------------------------------------------------------------------
function ECPevent, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common CommonCore, status, logmode, logunit

common EBSD_rawdata, accum_e, accum_z, mLPNH, mLPSH


if (EBSDdata.eventverbose eq 1) then help,event,/structure


WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

CASE eventval OF

  'ECPFORMAT': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.EBSDformatbgroup
                EBSDdata.imageformat = fix(val[0])
          endcase

  'ECPATTERNORIGIN': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.PatternOrigin
                EBSDdata.PatternOrigin = fix(val[0])
	  	ECPshowPattern,/single
		vals = ['Upper Left','Lower Left','Upper Right','Lower Right']
		  Core_Print, 'Pattern origin set to '+vals[EBSDdata.PatternOrigin]
	endcase
 'ECPATTERNSCALING': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.PatternScaling
                EBSDdata.PatternScaling = fix(val[0])
	  	ECPshowPattern,/single
		vals = ['linear', 'gamma']
		  Core_Print, 'Pattern scaling set to '+vals[EBSDdata.PatternScaling]
	endcase

 'CIRCULARMASK': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.circularmask
                EBSDdata.showcircularmask= fix(val[0])
	  	if (EBSDdata.Pmode eq 0) then EBSDshowPattern,/single else EBSDshowPattern
		vals = ['Off','On']
		  Core_Print, 'Circular mask set to '+vals[EBSDdata.showcircularmask]
	endcase

 'PMODE': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.Pmode
                EBSDdata.Pmode = fix(val[0])
		vals = ['Single Pattern','Angle File','Dictionary']
; next we need to turn on those widgets that belong to the selected mode (sensitivity=1)
		if (EBSDdata.Pmode eq 0) then begin
		  WIDGET_CONTROL, EBSDwidget_s.DisplayECP, sensitive=1
		  WIDGET_CONTROL, EBSDwidget_s.ECPgetanglefilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.PGdroplist, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.EBSDgetdictfilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
		    Core_Print, 'Pattern mode set to '+vals[EBSDdata.Pmode]
		end

		if (EBSDdata.Pmode eq 1) then begin
		  WIDGET_CONTROL, EBSDwidget_s.DisplayECP, sensitive=0
		  WIDGET_CONTROL, EBSDwidget_s.ECPgetanglefilename, sensitive=1
;	  WIDGET_CONTROL, EBSDwidget_s.PGdroplist, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.EBSDgetdictfilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
		    Core_Print, 'Pattern mode set to '+vals[EBSDdata.Pmode]
		end

		if (EBSDdata.Pmode eq 2) then begin
		  WIDGET_CONTROL, EBSDwidget_s.DisplayECP, sensitive=0
		  WIDGET_CONTROL, EBSDwidget_s.ECPgetanglefilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.PGdroplist, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.EBSDgetdictfilename, sensitive=0
;	  if ( (EBSDdata.Ncubochoric ne 0) and (EBSDdata.Dictpointgroup ne 0) and (EBSDdata.EBSDdictfilename ne '') ) then begin
;	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
;	  end else begin
;	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
;	  end
;	    Core_Print, 'Pattern mode set to '+vals[EBSDdata.Pmode]
		    Core_Print, 'Not implemented in this program Release',/blank
		end
	endcase

else: MESSAGE, "Event User Value Not Found"

endcase

return,eventval
end 
