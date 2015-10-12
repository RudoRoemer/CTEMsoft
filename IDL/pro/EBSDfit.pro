@Efitpreferences
@EBSDinit
@EBSDcalc

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
; EMsoft:EBSDfit.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDfit.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Electron backscatter diffraction detector parameter fit interface
;
;> @details 10/12/15 a new GUI to interactively determine the best fit parameters
;> for an EBSD pattern.
;
;> @date 10/12/15 MDG 1.0 first attempt at a user-friendly interface
;--------------------------------------------------------------------------
pro EBSDfit,dummy


common CommonCore, status, logmode, logunit


Efitwidget_s = {widgetstruct, $
	base:long(0), $                     	; base widget ID

}

Efitdata = {Efitdatastruct, $
}

; a few font strings (this will need to be redone for Windows systems)
fontstr='-adobe-new century schoolbook-bold-r-normal--14-100-100-100-p-87-iso8859-1'
fontstrlarge='-adobe-new century schoolbook-medium-r-normal--20-140-100-100-p-103-iso8859-1'
fontstrsmall='-adobe-new century schoolbook-medium-r-normal--14-100-100-100-p-82-iso8859-1'

;------------------------------------------------------------
; get the display window size to 80% of the current screen size (but be careful with double screens ... )
; We'll need to guess whether or not the user has a double screen: if the aspect ratio is larger than 16/9,
; then there are likely two screens, so we need to limit ourselves to just the first one...
; This should really become a core function that we can call from all programs.
device,decomposed = 0
device, GET_SCREEN_SIZE = scr

sar = float(scr[0])/float(scr[1])
if (sar gt (1.1*16.0/9.0)) then begin
	scr[0] = scr[0]/2
end
Efitdata.scrdimy = scr[1] * 0.8
Efitdata.scrdimx = scr[0]
Efitdata.xlocation = Efitdata.scrdimx / 8.0
Efitdata.ylocation = Efitdata.scrdimx / 8.0

;------------------------------------------------------------
; does the preferences file exist ?  If not, create it, otherwise read it
;Efitgetpreferences,/noprint

;------------------------------------------------------------
; create the top level widget
Efitwidget_s.base = WIDGET_BASE(TITLE='Electron Backscatter Diffraction Pattern Fit Program', $
                        /ROW, $
                        XSIZE=1220, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EfitDisplay_event', $
                        XOFFSET=Efitdata.xlocation, $
                        YOFFSET=Efitdata.ylocation)

;------------------------------------------------------------
; create the two main columns
; block 1 is the left column, with the logo 
block1 = WIDGET_BASE(Efitwidget_s.base, $
			/FRAME, $
			XSIZE=610, $
			/ALIGN_CENTER, $
			/COLUMN)

Efitwidget_s.logodraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			/ALIGN_CENTER, $
			XSIZE=600, $
			YSIZE=200)

; block 2 is the right column, with the input file widgets
; we're asking for the master pattern file name, which will load everything needed.
block2 = WIDGET_BASE(Efitwidget_s.base, $
			XSIZE=610, $
			/FRAME, $
			/COLUMN)


; fixed parameters
;   Detector tile angle 
;   Scintillator pixel size
;   Euler phi1 convention
;   Number of pixels
;   Beam current
;   Dwell time   
;   Circular mask on/off
;   Binning

; we'll also need some options for background subtraction

; and an option for gamma value scaling ?

; refinable parameters
;   Scintillator distance
;   Detector omega angle 
;   Detector pcx
;   Detector pcy
;   Euler phi1
;   Euler Phi
;   Euler phi2


















;------------------------------------------------------------
;------------------------------------------------------------
; then we have the program message window

Efitwidget_s.status= WIDGET_TEXT(block2, $
			XSIZE=90, $
			YSIZE=22, $
			/SCROLL, $
			VALUE=' ',$
			/ALIGN_LEFT)

; the following is needed by the Core_Print routine
status = Efitwidget_s.status 

;------------------------------------------------------------
;------------------------------------------------------------
; finally we need a couple of control buttons: GoFit and Quit

file11 = WIDGET_BASE(block2, XSIZE=590, /FRAME, /ROW)

Efitwidget_s.mainstop = WIDGET_BUTTON(file11, $
                                VALUE='Quit', $
                                UVALUE='Quit', $
                                EVENT_PRO='EfitDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

Efitwidget_s.gofit = WIDGET_BUTTON(file11, $
                                UVALUE='GOFIT', $
                                VALUE='GoFit', $
                                EVENT_PRO='EfitDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,Efitwidget_s.base,/REALIZE

; realize the draw widgets
WIDGET_CONTROL, Efitwidget_s.logodraw, GET_VALUE=drawID
Efitwidget_s.logodrawID = drawID
;
read_jpeg,'Resources/SEMONRlogo.jpg',logo
wset,Efitwidget_s.logodrawID
tvscl,logo,true=1

; and hand over control to the xmanager
XMANAGER,"EfitDisplay",Efitwidget_s.base,/NO_BLOCK








end ; program
