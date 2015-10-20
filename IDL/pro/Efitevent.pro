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
common FitParameters, nFit, fitName, defValue, fitValue, fitStep, fitOnOff, fitManualStep, fitManualUpDown, fitUserLabel, fitStepLabel, fitOnOffLabel, fitUpLabel, fitDownLabel, fitManualStepLabel, fitIterations

common EBSD_EMsoft, MCxtalname, MCmode, nsx, nsy, EkeV, Ehistmin, Ebinsize, depthmax, depthstep, MCsig, MComega, $
                    numEbins, numzbins, accum_e, accum_z, Masterenergyfile, npx, npy, nnE, numset, mLPNH, mLPSH, Masterxtalname, expEBSDpattern, EBSDpattern

WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

CASE eventval OF
        'FITMODE' : begin
                oldmode = Efitdata.fitmode
                Efitdata.fitmode = Core_WidgetChoiceEvent( Efitwidget_s.fitmode,  'Fit mode? ',/value)
                if ((oldmode ne 0) and (Efitdata.fitmode eq 0)) then begin
; reset the fitOnOff parameters to off for all of them...
                  for i=0,nFit-1 do WIDGET_CONTROL, set_value=0, Efitwidget_s.fitOnOff[i]
                endif
                if (Efitdata.fitmode eq 1) then begin
                  for i=0,4 do WIDGET_CONTROL, set_value=1, Efitwidget_s.fitOnOff[i]
                  for i=5,nFit-1 do WIDGET_CONTROL, set_value=0, Efitwidget_s.fitOnOff[i]
                  fitOnOff[0:4] = 1L
                  fitOnOff[5:*] = 0L
                endif
                if (Efitdata.fitmode eq 2) then begin
                  for i=0,4 do WIDGET_CONTROL, set_value=0, Efitwidget_s.fitOnOff[i]
                  for i=5,nFit-1 do WIDGET_CONTROL, set_value=1, Efitwidget_s.fitOnOff[i]
                  fitOnOff[0:4] = 0L
                  fitOnOff[5:*] = 1L
                endif
        endcase

        'CONVCRIT' : begin
                Efitdata.convcrit = Core_WidgetChoiceEvent( Efitwidget_s.convcrit,  'Convergence criterion? ')
        endcase

        'SMOOTHVAL' : begin
                Efitdata.smoothval = Core_WidgetChoiceEvent( Efitwidget_s.smoothval,  'Smoothing parameter? ')
                if ((max(EBSDpattern) gt 0) or (max(expEBSDpattern) gt 0)) then Efit_showpattern
        endcase

        'RAMPONOFF' : begin
                Efitdata.ramponoff = Core_WidgetChoiceEvent( Efitwidget_s.ramponoff,  'Ramp filter? ')
                if ((max(EBSDpattern) gt 0) or (max(expEBSDpattern) gt 0)) then Efit_showpattern
        endcase

        'HIPASSONOFF' : begin
                Efitdata.hipassonoff = Core_WidgetChoiceEvent( Efitwidget_s.hipassonoff,  'Hipass filter? ')
                if ((max(EBSDpattern) gt 0) or (max(expEBSDpattern) gt 0)) then Efit_showpattern
        endcase

        'CIRCULARMASK' : begin
                Efitdata.showcircularmask = Core_WidgetChoiceEvent( Efitwidget_s.circularmask,  'Add circular mask? ')
                if ((max(EBSDpattern) gt 0) or (max(expEBSDpattern) gt 0)) then Efit_showpattern
	endcase

        'BINNING' : begin
                Efitdata.detbinning = Core_WidgetChoiceEvent( Efitwidget_s.detbinning,  'Binning set to  ')
                if ((max(EBSDpattern) gt 0) or (max(expEBSDpattern) gt 0)) then Efit_showpattern
	endcase

        'EBSDEULERCONVENTION' : begin
                Efitdata.EulerConvention = Core_WidgetChoiceEvent( Efitwidget_s.EulerConvention,  'Euler Convention set to ')
	endcase

        'EBSDPATTERNORIGIN' : begin
                Efitdata.PatternOrigin = Core_WidgetChoiceEvent( Efitwidget_s.PatternOrigin,  'Pattern origin set to ')
                if ((max(EBSDpattern) gt 0) or (max(expEBSDpattern) gt 0)) then Efit_showpattern
	endcase

        'DEToL' : begin
                Efitdata.detoL = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[0],  'Fit scintillator distance? ')
                if (Efitdata.detoL eq 0) then fitOnOff[0] = 0 else fitOnOff[0] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DEToOMEGA' : begin
                Efitdata.detoomega = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[1],  'Fit sample omega angle? ')
                if (Efitdata.detoL eq 0) then fitOnOff[1] = 0 else fitOnOff[1] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DEToXPC' : begin
                Efitdata.detoxpc = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[2],  'Fit pattern center x? ')
                if (Efitdata.detoL eq 0) then fitOnOff[2] = 0 else fitOnOff[2] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DEToYPC' : begin
                Efitdata.detoypc = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit pattern center y? ')
                if (Efitdata.detoL eq 0) then fitOnOff[3] = 0 else fitOnOff[3] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DEToGAMMA' : begin
                Efitdata.detogamma = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[4],  'Fit intensity gammma? ') 
                if (Efitdata.detoL eq 0) then fitOnOff[4] = 0 else fitOnOff[4] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DETophi1' : begin
                Efitdata.detophi1 = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit Euler phi1 angle? ')
                if (Efitdata.detoL eq 0) then fitOnOff[5] = 0 else fitOnOff[5] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DETophi' : begin
                Efitdata.detophi = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit Euler Phi angle? ')
                if (Efitdata.detoL eq 0) then fitOnOff[6] = 0 else fitOnOff[6] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DETophi2' : begin
                Efitdata.detophi2 = Core_WidgetChoiceEvent( Efitwidget_s.fitOnOff[3],  'Fit Euler phi2 angle? ')
                if (Efitdata.detoL eq 0) then fitOnOff[7] = 0 else fitOnOff[7] = 1
                if (total(fitOnOff) gt 0) then begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =1
                end else begin
                  WIDGET_CONTROL, Efitwidget_s.mkjson, sensitive =0
                endelse
	endcase

        'DISPLAYOPTION' : begin   ; this comes from the Efit_control widget...
                Efitdata.displayoption = Core_WidgetChoiceEvent( Efitwidget_s.displayoption, 'Set option to ',/value)
                Efit_showpattern
        endcase

        'PATTERNFORMAT' : begin
                Efitdata.imageformat = Core_WidgetChoiceEvent( Efitwidget_s.imageformat, 'Set option to ',/value)
        endcase

else: MESSAGE, "Efitevent: Event User Value Not Found"

endcase

return,eventval
end 
