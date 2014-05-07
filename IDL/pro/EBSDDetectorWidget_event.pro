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
; CTEMsoft2013:EBSDDetectorWidget_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDDetectorWidget_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main event handler for detector widget
;
;> @date 05/07/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDDetectorWidget_event, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall


if (EBSDdata.eventverbose eq 1) then help,event,/structure

; intercept the detector widget movement here 
if (event.id eq EBSDwidget_s.detectorbase) then begin
  EBSDdata.Detectorxlocation = event.x
  EBSDdata.Detectorylocation = event.y-25
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

  CASE eventval OF

 'DETL': EBSDdata.detL = Core_WidgetEvent( EBSDwidget_s.detL,  'Detector distance set to [micron] ', '(F9.2)', /flt)

 'DETTHETA': EBSDdata.dettheta = Core_WidgetEvent( EBSDwidget_s.dettheta, 'Detector angle set to [deg] ', '(F6.2)', /flt)

 'DETDELTA': EBSDdata.detdelta = Core_WidgetEvent( EBSDwidget_s.detdelta, 'Scintillator pixel size set to [micron] ', '(F6.2)', /flt)

 'DETNUMSX': EBSDdata.detnumsx = Core_WidgetEvent( EBSDwidget_s.detnumsx, 'Scintillator number of pixels along x set to [micron] ', '(I6)', /lng)
 'DETNUMSY': EBSDdata.detnumsy = Core_WidgetEvent( EBSDwidget_s.detnumsy, 'Scintillator number of pixels along y set to [micron] ', '(I6)', /lng)

 'DETXPC': EBSDdata.detxpc = Core_WidgetEvent( EBSDwidget_s.detxpc, 'Pattern Center x-coordinate set to [pixels] ', '(F7.2)', /flt)
 'DETYPC': EBSDdata.detypc = Core_WidgetEvent( EBSDwidget_s.detypc, 'Pattern Center y-coordinate set to [pixels] ', '(F7.2)', /flt)

 'DETBINNING': EBSDdata.detbinning = Core_WidgetEvent( EBSDwidget_s.detbinning, 'CCD binning factor set to ', '(I3)', /lng)

 'DETBEAMCURRENT': EBSDdata.detbeamcurrent = Core_WidgetEvent( EBSDwidget_s.detbeamcurrent, 'Beam current set to [A] ', '(D9.2)', /flt)
 'DETDWELLTIME': EBSDdata.detdwelltime = Core_WidgetEvent( EBSDwidget_s.detdwelltime, 'Dwell time set to [s] ', '(D9.2)', /flt)

 'DETphi1': EBSDdata.detphi1 = Core_WidgetEvent( EBSDwidget_s.detphi1, 'Euler angle phi1 set to [deg] ', '(F6.2)', /flt)
 'DETPhi': EBSDdata.detphi = Core_WidgetEvent( EBSDwidget_s.detphi, 'Euler angle Phi set to [deg] ', '(F6.2)', /flt)
 'DETphi2': EBSDdata.detphi2 = Core_WidgetEvent( EBSDwidget_s.detphi2, 'Euler angle phi2 set to [deg] ', '(F6.2)', /flt)

 'EBSDMINENERGYLIST': begin
		EBSDdata.Eminsel = fix(event.index)
		if (EBSDdata.Eminsel gt EBSDdata.Emaxsel) then begin
			EBSDdata.Emaxsel = EBSDdata.Eminsel
			WIDGET_CONTROL, set_droplist_select = EBSDdata.Emaxsel, EBSDwidget_s.EBSDmaxenergylist
		end
		WIDGET_CONTROL, set_droplist_select = EBSDdata.Eminsel, EBSDwidget_s.EBSDminenergylist
	 endcase

 'EBSDMAXENERGYLIST': begin
		EBSDdata.Emaxsel = fix(event.index)
		if (EBSDdata.Emaxsel lt EBSDdata.Eminsel) then begin
			EBSDdata.Eminsel = EBSDdata.Emaxsel
			WIDGET_CONTROL, set_droplist_select = EBSDdata.Eminsel, EBSDwidget_s.EBSDminenergylist
		end
		WIDGET_CONTROL, set_droplist_select = EBSDdata.Emaxsel, EBSDwidget_s.EBSDmaxenergylist
	 endcase


 'CLOSEDETECTOR': begin
; kill the base widget
		WIDGET_CONTROL, EBSDwidget_s.detectorbase, /DESTROY
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 
