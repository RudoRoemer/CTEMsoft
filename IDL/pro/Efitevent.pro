;
; Copyright (c) 2015, Marc De Graef/Carnegie Mellon University
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
; EMsoft:Efitevent.pro
;--------------------------------------------------------------------------
;
; PROGRAM: Efitevent.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 03/19/14 MDG 1.0 first version
;--------------------------------------------------------------------------
function Efitevent, event

;------------------------------------------------------------
; common blocks
common Efit_widget_common, Efitwidget_s
common Efit_data_common, Efitdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common CommonCore, status, logmode, logunit

WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

CASE eventval OF
        'DEToL' : begin
                Efitdata.detoL = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[0],  'Fit scintillator distance? ')
	endcase
        'DEToOMEGA' : begin
                Efitdata.detoomega = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[1],  'Fit sample omega angle? ')
	endcase
        'DEToXPC' : begin
                Efitdata.detoxpc = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[2],  'Fit pattern center x? ')
	endcase
        'DEToYPC' : begin
                Efitdata.detoypc = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit pattern center y? ')
	endcase
        'DEToGAMMA' : begin
                Efitdata.detogamma = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[4],  'Fit intensity gammma? ') 
	endcase
        'DETophi1' : begin
                Efitdata.detophi1 = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit Euler phi1 angle? ')
	endcase
        'DETophi' : begin
                Efitdata.detophi = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit Euler Phi angle? ')
	endcase
        'DETophi2' : begin
                Efitdata.detophi2 = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit Euler phi2 angle? ')
	endcase
        'DISPLAYOPTION' : begin   ; this comes from the Efit_control widget...
                Efitdata.displayoption = Core_WidgetChoiceEvent( Efitwidget_s.displayoption, 'Set option to ',/value)
        endcase
        'PATTERNFORMAT' : begin
                Efitdata.imageformat = Core_WidgetChoiceEvent( Efitwidget_s.imageformat, 'Set option to ',/value)
        endcase

else: MESSAGE, "Efitevent: Event User Value Not Found"

endcase

return,eventval
end 
