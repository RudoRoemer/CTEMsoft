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

 'DETNUMSX': EBSDdata.detnumsx = Core_WidgetEvent( EBSDwidget_s.detnumsx, 'Scintillator number of pixels along x set to ', '(I4)', /lng)
 'DETNUMSY': EBSDdata.detnumsy = Core_WidgetEvent( EBSDwidget_s.detnumsy, 'Scintillator number of pixels along y set to ', '(I4)', /lng)

 'DETXPC': EBSDdata.detxpc = Core_WidgetEvent( EBSDwidget_s.detxpc, 'Pattern Center x-coordinate set to [pixels] ', '(F7.2)', /flt)
 'DETYPC': EBSDdata.detypc = Core_WidgetEvent( EBSDwidget_s.detypc, 'Pattern Center y-coordinate set to [pixels] ', '(F7.2)', /flt)

 'DETBEAMCURRENT': EBSDdata.detbeamcurrent = Core_WidgetEvent( EBSDwidget_s.detbeamcurrent, 'Beam current set to [nA] ', '(F7.3)', /flt)
 'DETDWELLTIME': EBSDdata.detdwelltime = Core_WidgetEvent( EBSDwidget_s.detdwelltime, 'Dwell time set to [mu s] ', '(F7.3)', /flt)

 'DETax1': EBSDdata.detax1 = Core_WidgetEvent( EBSDwidget_s.detax1, 'Axis-angle entry 1 set to ', '(F6.2)', /flt)
 'DETax2': EBSDdata.detax2 = Core_WidgetEvent( EBSDwidget_s.detax2, 'Axis-angle entry 2 set to ', '(F6.2)', /flt)
 'DETax3': EBSDdata.detax3 = Core_WidgetEvent( EBSDwidget_s.detax3, 'Axis-angle entry 3 set to ', '(F6.2)', /flt)
 'DETax4': EBSDdata.detax4 = Core_WidgetEvent( EBSDwidget_s.detax4, 'Axis-angle entry 4 set to [deg] ', '(F6.2)', /flt)

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

 'GETEBSDFILENAME': begin
; display a filesaving widget in the data folder with the file extension filled in
	  filename = DIALOG_PICKFILE(/write,path=EBSDdata.pathname,title='enter EBSD output file name ')
	  if (filename ne '') then begin
	    EBSDdata.EBSDpatternfilename = filename
	    WIDGET_CONTROL, set_value=filename, EBSDwidget_s.EBSDpatternfilename
	    WIDGET_CONTROL, EBSDwidget_s.DisplayEBSD, sensitive=1
	  end
	endcase


 'DETphi1': EBSDdata.detphi1 = Core_WidgetEvent( EBSDwidget_s.detphi1, 'Euler angle phi1 set to [deg] ', '(F6.2)', /flt)
 'DETPhi': EBSDdata.detphi = Core_WidgetEvent( EBSDwidget_s.detphi, 'Euler angle Phi set to [deg] ', '(F6.2)', /flt)
 'DETphi2': EBSDdata.detphi2 = Core_WidgetEvent( EBSDwidget_s.detphi2, 'Euler angle phi2 set to [deg] ', '(F6.2)', /flt)

 'GETANGLEFILENAME': begin
; display a filesaving widget in the data folder with the file extension filled in
	  filename = DIALOG_PICKFILE(/write,path=data.pathname,title='Select angle input file')
	  if (filename ne '') then begin
	    EBSDdata.EBSDanglefilename = filename
	    WIDGET_CONTROL, set_value=filename, EBSDwidget_s.EBSDanglefilename
	    EBSDreadanglefile,filename
	  end
	endcase

 'DICTIONARYPG': begin
	  EBSDdata.Dictpointgroup = event.index
	  if ( (EBSDdata.Ncubochoric ne 0) and (EBSDdata.EBSDdictfilename ne '') ) then begin
	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=1
	  end
	endcase

 'NCUBOCHORIC': begin 
	  EBSDdata.Ncubochoric = Core_WidgetEvent( EBSDwidget_s.Ncubochoric,  'Number of smapling points along cube semi-edge set to ', '(I4)', /lng)
	  if ( (EBSDdata.Dictpointgroup ne 0) and (EBSDdata.EBSDdictfilename ne '') ) then begin
	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=1
	  end
	endcase

 'GETDICTFILENAME': begin
; display a filesaving widget 
	  filename = DIALOG_PICKFILE(/write,path=EBSDdata.pathname,title='Set dictionary angle file name ')
	  if (filename ne '') then begin
	    EBSDdata.EBSDdictfilename = filename
	    WIDGET_CONTROL, set_value=filename, EBSDwidget_s.EBSDdictfilename
	  end
	  if ( (EBSDdata.Dictpointgroup ne 0) and (EBSDdata.Ncubochoric ne 0) ) then begin
	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=1
	  end
	endcase

 'GODICT': begin
; first we need to make sure that the path to the fortran executables is known... this is stored in the 
; preferences file, but is initially set to 'path_unknown'
	  if (EBSDdata.f90exepath eq 'path_unknown') then begin
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
  		Core_Print, 'PLEASE SELECT THE PATH TO THE f90 EXECUTABLE FOLDER', /blank
		directoryname = DIALOG_PICKFILE(/write,path=EBSDdata.pathname,/directory,title='Select the f90 executable folder')
	        EBSDdata.f90exepath = directoryname[0]
	  end
print,'executable pathname = ',EBSDdata.f90exepath
	endcase

 'DISPLAYEBSD': begin

; this does two things.  First of all, the CTEMEBSD program is called with the current
; parameters for the detector and microscope geometry, and the single set of Euler angles
; entered in the 'Detector and Pattern Mode Widget'.
;
; Then, when the CTEMEBSD program has produced its output file, we create a new widget
; that displays this raw EBSD pattern; the user can then apply a number of intensity
; corrections, as well as binning, to optimize the pattern quality.  Once this is done, 
; those parameters will be used for the Angle File mode and the Dictionary mode.

; first, create the nml file and execute the CTEMEBSD program
	  status = 0
	  EBSDExecute,status,/single

; then we create the EBSDpattern widget and let the user adjust the imaging parameters
	  if (status eq 1) then begin
	    if (XRegistered("EBSDPatternWidget") EQ 0) then EBSDPatternWidget,/single else EBSDshowPattern,/single
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