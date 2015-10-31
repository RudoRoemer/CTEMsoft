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
; EMsoft:ECPDetectorWidget_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: ECPDetectorWidget_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main event handler for detector widget
;
;> @date 10/30/15 MDG 1.0 first version
;--------------------------------------------------------------------------
pro ECPDetectorWidget_event, event

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

 'DETW': EBSDdata.detW = Core_WidgetEvent( EBSDwidget_s.detW,  'Working distance set to [mm] ', '(F9.2)', /flt)
 'DETRI': EBSDdata.detRi = Core_WidgetEvent( EBSDwidget_s.detRi,  'Detector inner radius set to [mm] ', '(F9.2)', /flt)
 'DETRO': EBSDdata.detRo = Core_WidgetEvent( EBSDwidget_s.detRo,  'Detector outer radius set to [mm] ', '(F9.2)', /flt)

 'DETSAMPLEYTILT': EBSDdata.detsampleytilt = Core_WidgetEvent( EBSDwidget_s.detsampleytilt, 'Sample y-tilt angle set to [deg] ', '(F6.2)', /flt)
 'DETNUMSX': begin
                EBSDdata.detnumsx = Core_WidgetEvent( EBSDwidget_s.detnumsx, 'Scintillator number of pixels along x set to ', '(I4)', /lng)
                EBSDdata.detnumsy = EBSDdata.detnumsx
        endcase

 'DETTHETAC': EBSDdata.detthetac = Core_WidgetEvent( EBSDwidget_s.detthetac, 'Incident beam cone semi-angle to [deg] ', '(F6.2)', /flt)

 'DETBEAMCURRENT': EBSDdata.detbeamcurrent = Core_WidgetEvent( EBSDwidget_s.detbeamcurrent, 'Beam current set to [nA] ', '(F9.2)', /flt)
 'DETDWELLTIME': EBSDdata.detdwelltime = Core_WidgetEvent( EBSDwidget_s.detdwelltime, 'Dwell time set to [mu s] ', '(F9.2)', /flt)

 'DETphi1': EBSDdata.detphi1 = Core_WidgetEvent( EBSDwidget_s.detphi1, 'Euler angle phi1 set to [deg] ', '(F6.2)', /flt)
 'DETPhi': EBSDdata.detphi = Core_WidgetEvent( EBSDwidget_s.detphi, 'Euler angle Phi set to [deg] ', '(F6.2)', /flt)
 'DETphi2': EBSDdata.detphi2 = Core_WidgetEvent( EBSDwidget_s.detphi2, 'Euler angle phi2 set to [deg] ', '(F6.2)', /flt)

 'GETANGLEFILENAME': begin
; display a filesaving widget in the data folder with the file extension filled in
	  filename = DIALOG_PICKFILE(/write,path=EBSDdata.pathname,title='Select angle input file')
	  if (filename ne '') then begin
	    EBSDdata.ECPanglefilename = filename
	    WIDGET_CONTROL, set_value=filename, EBSDwidget_s.ECPanglefilename
	    EBSDreadanglefile,filename
	    WIDGET_CONTROL, EBSDwidget_s.GoAngle, sensitive=1
	  end
	endcase

 'DISPLAYECP': begin
; first we need to make sure that the path to the fortran executables is known... this is stored in the 
; preferences file, but is initially set to 'path_unknown'
	  if (EBSDdata.EMsoftpathname eq 'path_unknown') then begin
                EBSDdata.EMsoftpathname = Core_getenv()+'Build/Bin'
                Core_print,'exeutable path set to '+EBSDdata.EMsoftpathname 
	  end

; is the correct widget up on the screen ?
	  if XRegistered("ECPPatternWidget") then begin
	    if (EBSDdata.currentdisplaywidgetmode ne 0) then WIDGET_CONTROL, EBSDwidget_s.patternbase, /DESTROY
	  end

; first we need to set up the array structures to do a call_external of the SingleECPPatternWrapper
; routine; and then we display the pattern in a new widget

; first, set up the variables and do a call_external
	  status = 0
	  ECPExecute,status,ECpattern,/single

; then we create the EBSDpattern widget and let the user adjust the imaging parameters
	  if (status eq 1) then begin
	    if (XRegistered("EBSDPatternWidget") EQ 0) then ECPatternWidget,ECpattern,/single else ECPshowPattern,ECpattern,/single
	  end

	endcase

 'GOANGLE': begin
; first we need to make sure that the path to the fortran executables is known... this is stored in the 
; preferences file, but is initially set to 'path_unknown'
	  if (EBSDdata.f90exepath eq 'path_unknown') then begin
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
		directoryname = DIALOG_PICKFILE(/write,path=EBSDdata.pathname,/directory,title='Select the f90 executable folder')
	        EBSDdata.f90exepath = directoryname[0]
	  end

; is the correct widget up on the screen ?
	  if XRegistered("EBSDPatternWidget") then begin
	    if (EBSDdata.currentdisplaywidgetmode ne 1) then WIDGET_CONTROL, EBSDwidget_s.patternbase, /DESTROY
	  end

; this does two things.  First of all, the CTEMEBSD program is called with the current
; parameters for the detector and microscope geometry, and the angle file
;
; Then, when the CTEMEBSD program has produced its output file, we create a new widget
; that displays these EBSD patterns; the user can then save selected patterns or all patterns.
; At this point, there is no option to change the imaging parameters; all the settings of the 
; other parts of the widget apply to this pattern calculation

; first, create the nml file and execute the CTEMEBSD program
	  status = 0
	  EBSDExecute,status

; then we create the EBSDpattern widget and let the user adjust the imaging parameters
	  if (status eq 1) then begin
	    if (XRegistered("EBSDPatternWidget") EQ 0) then EBSDPatternWidget else EBSDshowPattern
	  end

	endcase

 'GODICTIONARY': begin
; first we need to make sure that the path to the fortran executables is known... this is stored in the 
; preferences file, but is initially set to 'path_unknown'
	  if (EBSDdata.f90exepath eq 'path_unknown') then begin
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
		directoryname = DIALOG_PICKFILE(/write,path=EBSDdata.pathname,/directory,title='Select the f90 executable folder')
	        EBSDdata.f90exepath = directoryname[0]
	  end

; is the correct widget up on the screen ?
	  if XRegistered("EBSDPatternWidget") then begin
	    if (EBSDdata.currentdisplaywidgetmode ne 1) then WIDGET_CONTROL, EBSDwidget_s.patternbase, /DESTROY
	  end

; this does two things.  First of all, the CTEMEBSD program is called with the current
; parameters for the detector and microscope geometry, and the angle file
;
; Then, when the CTEMEBSD program has produced its output file, we create a new widget
; that displays these EBSD patterns; the user can then save selected patterns or all patterns.
; At this point, there is no option to change the imaging parameters; all the settings of the 
; other parts of the widget apply to this pattern calculation

; first, create the nml file and execute the CTEMEBSD program
	  status = 0
	  EBSDExecute,status

; then we create the EBSDpattern widget and let the user adjust the imaging parameters
	  if (status eq 1) then begin
	    if (XRegistered("EBSDPatternWidget") EQ 0) then EBSDPatternWidget else EBSDshowPattern
	  end

	endcase

 'CLOSEDETECTOR': begin
; kill the base widget
		WIDGET_CONTROL, EBSDwidget_s.detectorbase, /DESTROY
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 
