@EBSDDisplay_event    		; EBSD event handler
@EBSDgetfilename		; select a geometry file
@EBSDreaddatafile		; read geometry and data files
@EBSDMCDisplayWidget		; MC display widget
@EBSDMCDisplayWidget_event	; MC display widget event handler
@EBSDDetectorWidget		; detector widget
@EBSDDetectorWidget_event	; detector widget event handler
@EBSDevent			; event handler for button groups
@EBSDshowMC			; display a Lambert projection image
@EBSDgetpreferences		; read preferences file
@EBSDwritepreferences		; write preferences file
@Core_LambertS2C		; modified Lambert to Lambert projection
@Core_LambertS2SP		; modified Lambert to stereographic projection
@Core_colorwheel		; color representation of energy distribution
@Core_WText			; generate a text widget with a label
@Core_Print			; print messages to status window and log file
@Core_WidgetEvent		; general data handler for various widget events
;
; Copyright (c) 2014, Marc De Graef/Carnegie Mellon University
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
; CTEMsoft2013:EBSDDisplay.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDDisplay.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Electron backscatter diffraction pattern display
;
;> @details 03/19/14 This is the very first implementation of the EBSD visualization GUI.
;> At this stage, we allow the user to display various versions of the Monte Carlo
;> results, as well as the master EBSD pattern in a number of projection modes.
;> The user can then define a series of detector parameters and compute EBSD
;> patterns.  In this version, there are no PSF, scintillator energy dependence,
;> or microscope geometric distortions; only the detector geometry, noise (on/off),
;> binning, and brightness/contrast controls are available.  The other options 
;> will be included in the next version (as will be the dictionary generation option).
;
;> @date 01/27/14 MDG 1.0 first attempt at a user-friendly interface
;> @date 03/17/14 MDG 1.1 main widget rewrite; prepended 'EBSD' to widget_s and data structures for future merge with other IDL routines
;> @date 03/19/14 MDG 1.2 implementation of Monte Carlo and master EBSD widgets
;--------------------------------------------------------------------------
pro EBSDDisplay,dummy
;
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall
common PointGroups, PGTHD, PGTWD, DG

common CommonCore, status, logmode, logunit

common projections, mcxcircle, mcycircle, mpxcircle, mpycircle, mcSPxcircle, mcSPycircle, mpSPxcircle, mpSPycircle 

PGTWD = [ 'none',' 1',' 2',' m',' 2mm',' 4',' 4mm',' 3',' 3m1',' 31m',' 6',' 6mm'] 
PGTHD = ['  ' ,' 1',' -1',' 2',' m',' 2/m',' 222', $
         ' mm2',' mmm',' 4',' -4',' 4/m',' 422', $
         ' 4mm',' -42m','4/mmm',' 3',' -3',' 32', $
         ' 3m',' -3m',' 6',' -6',' 6/m',' 622', $
         ' 6mm',' -6m2',' 6/mmm',' 23',' m3',' 432', $
         ' -43m',' m-3m']
DG = ['  ',' 1',' 1R',' 2',' 2R',' 21R','  mR', $
      ' m',' m1R',' 2mRmR',' 2mm',' 2RmmR',' 2mm1R', $
      ' 4',' 4R',' 41R',' 4mRmR',' 4mm',' 4RmmR', $
      ' 4mm1R',' 3',' 6R',' 3mR',' 3m',' 6RmmR', $
      ' 6',' 31R',' 61R',' 6mRmR',' 6mm',' 3m1R', $
      ' 6mm1R']

!EXCEPT=0

;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("EBSDDisplay") NE 0) then begin
  print,'EBSDDisplay is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; define a few structures (one for widgets, and one for data)
EBSDwidget_s = {widgetstruct, $
	base:long(0), $                     	; base widget ID

	; Monte Carlo widget ids
	mcfilename: long(0), $			; Monte Carlo filename
	mcfilesize: long(0), $			; file size in bytes
	mcenergymin: long(0), $			; minimum energy
	mcenergymax: long(0), $			; maximum energy
	mcenergybinsize: long(0), $		; energy bin size
	mcenergynumbin: long(0), $		; number of energy bins
	voltage: long(0), $			; microscope voltage
	mcdepthmax: long(0), $			; maximum penetration depth
	mcdepthstep: long(0), $			; step size for depth
	mcdepthnumbins: long(0), $		; number of depth bins
	mcimx: long(0), $			; N in 2N+1 (x-pixels)
	mcimy: long(0), $			; N in 2N+1 (y-pixels)
	mctotale: long(0), $			; total number of incident electrons
	mcbse: long(0), $			; total number of bacck-scattered electrons
	mcvangle: long(0), $			; vertical tilt angle
	mchangle: long(0), $			; horizontal tilt angle
	mcmode: long(0), $			; Monte Carlo mode
	mcloadfile: long(0), $                  ; load file button
	mcdisplay: long(0), $                   ; MC display button
	   
	; Master Pattern widget ids
	mpfilename: long(0), $			; file name for master pattern
	mpfilesize: long(0), $			; file size in bytes
	mpimx: long(0), $			; x-pixels N as in (2N+1)
	mpimy: long(0), $			; y-pixels (should be equal to mpimx)
	mpgridmode: long(0), $			; Lambert grid mode
	asymunit: long(0), $			; widget for asymmetric unit position selection
	xtalname: long(0), $			; file name for crystal structure data
	mploadfile: long(0), $                  ; load file button

	; detector parameter widget ids
	detL: long(0), $			; sample scintillator distance [microns]
	dettheta: long(0), $			; detector tilt angle [degrees]
	detdelta: long(0), $			; detector pixel size [microns]
	detnumsx: long(0), $			; number of x-pixels
	detnumsy: long(0), $			; number of y-pixels
	detxpc: long(0), $			; x -pattern center [pixels]
	detypc: long(0), $			; y -pattern center [pixels]
	detbinning: long(0), $			; binning
	detbeamcurrent: long(0), $		; beam current [A]
	detdwelltime: long(0), $		; dwell time [s]
	detphi1: long(0), $			; phi1 Euler angle 
	detphi: long(0), $			; phi Euler angle 
	detphi2: long(0), $			; phi2 Euler angle 
	EulerConvention: long(0), $		; Euler angle convention for phi2 (TSL or HKL)
	BGmode: long(0), $			; background/full display mode
	EBSDminenergylist: long(0), $		; min energy widget
	EBSDmaxenergylist: long(0), $		; max energy widget
	PatternOrigin: long(0), $		; pattern origin widget

	; other collected items
	MCdisplaybase: long(0), $		; Monte Carlo display base
	MPdisplaybase: long(0), $		; Master Pattern & Monte Carlo display base
	detectorbase: long(0), $		; detector base widget
	status:long(0), $                       ; status window
	logfile: long(0), $			; logfile toggle widget ID
	detector:long(0), $                     ; detector widget
	MCbutton:long(0), $                     ; MC button ID
	MCslider:long(0), $                     ; MC slider ID
	MCenergyval: long(0), $			; MC energy value
        MCmin: long(0), $                       ; MC display minimum
        MCmax: long(0), $                       ; MC display maximum
        MPbutton:long(0), $                     ; MP button ID
        MPmin: long(0), $                       ; MP display minimum
        MPmax: long(0), $                       ; MP display maximum
	MCDraw:long(0), $                       ; pattern draw widget
	MCDrawID:long(0), $                     ; pattern draw widget
	MPDraw:long(0), $                       ; pattern draw widget
	MPDrawID:long(0), $                     ; pattern draw widget
	mainstop:long(0), $                     ; program quit button
	logodraw:long(0), $                     ; logo widget 
	logodrawID:long(0),$                    ; logo widget ID
	EBSDformatbgroup: long(0), $		; image file format widget
	MCLambertSelector:long(0),$             ; Lambert Selector widget ID
	MPLambertSelector:long(0),$             ; Lambert Selector widget ID
	MCLambertMode: long(0) $		; Lambert sum or individual image mode
           }

EBSDdata = {EBSDdatastruct, $
	; Monte Carlo parameters first 
	mcfilename: '', $			; Monte Carlo result file name
	mcfilesize: long64(0), $		; Monte Carlo file size [bytes]
	mcenergymin: float(0.0), $		; minimum energy
	mcenergymax: float(0.0), $		; maximum energy
	mcenergybinsize: float(0.0), $		; energy bin size
	mcenergynumbin: long(0), $		; number of energy bins
	voltage: float(0.0), $			; microscope voltage
	mcdepthmax: float(0.0), $		; maximum depth in MC file
	mcdepthstep: float(0.0), $		; depth step size
	mcdepthnumbins: long(0), $		; number of depth bins
	mcimx: long(0), $			; number of pixels along x in modified Lambert map
	mcimy: long(0), $			; same along y
	mctotale: long(0), $			; total number of electrons hitting the sample
	mcbse: long(0), $			; total number of BSE electrons
	mcvangle: float(0.0), $			; vertical sample tilt angle (around TD)
	mchangle: float(0.0), $			; horizontal sample tilt angle (around RD)
	mcmode: '', $				; 'CSDA' (continuous slowing down approximation) or 'DLOS' (discrete losses)
	Esel: long(0), $			; energy selection for slider in MC and MP display routines
	mcprogname: '', $ 			; MC program name
	mcscversion: '', $ 			; source code version number
	mcdataedims: lon64arr(3), $		; dimensions of accum_e
	mcdatazdims: lon64arr(4), $		; dimensions of accum_z

	; then Master Pattern parameters
	mpfilename: '', $ 			; master pattern file name
	mpfilesize: long(0), $			; size (in bytes) of master pattern file
	mpimx: long(0), $			; number of x-pixels in master pattern (N in 2N+1)
	mpimy: long(0), $			; same along y
	numset: long(0), $			; number of positions in asymmetric unit
	Asymsel: long(-1), $			; which asymmetric unit position to display?
	atnum: lonarr(250), $			; number of atoms in asymmetric unit 
	mpenergynumbin: long(0), $		; number of energy bins (may be different from MC file)
	mpgridmode: '', $			; 'hex' or 'squ' for Lambert grid type
	xtalname: '', $				; crystal structure filename
	mpprogname: '', $ 			; Master Pattern program name
	mpscversion: '', $ 			; source code version number
	mpdatadims: lon64arr(3), $		; dimensions of raw data array

	; detector parameters
	detL: float(0), $			; scintillator - sample distance [microns]
	dettheta: float(0), $			; detector tilt angle [degrees]
	detdelta: float(0), $			; scintillator pixel size [microns]
	detnumsx: long(0), $			; number of x-pixels
	detnumsy: long(0), $			; number of y-pixels
	detxpc: float(0), $			; x-pattern center [pixels]
	detypc: float(0), $			; y-pattern center [pixels]
	detbinning: long(0), $			; binning
	detbeamcurrent: float(0), $		; beam current [A]
	detdwelltime: float(0), $		; dwell time [s]
	detphi1: float(0), $			; phi1 Euler angle 
	detphi: float(0), $			; phi Euler angle 
	detphi2: float(0), $			; phi2 Euler angle 
	EulerConvention: long(0), $		; Euler angle convention (TSL = 0, HKL = 2)
	BGmode: long(0), $			; background/full pattern display mode
	Eminsel: long(0), $			; min energy selection 
	Emaxsel: long(0), $			; max energy selection 
	PatternOrigin: long(0), $		; pattern origin indicator

	; then general program parameters
	eventverbose: fix(0), $			; used for event debugging (0=off, 1=on)
	scversion: '', $			; source code version number
	pathname: '', $				; pathname (obviously)
	suffix: '', $				; filename suffix 
	homefolder: '', $			; startup folder of the program
	EBSDroot: 'undefined', $		; current pathname (is stored in preferences file)
	prefname: '~/.EBSDgui.prefs', $		; filename of preferences file (including path)
	nprefs: fix(0), $			; number of preferences in file
	MCLSmode: fix(0), $			; Monte Carlo Lambert Selector tag
	MCLSum: fix(0), $			; Monte Carlo Lambert sum or individual tag
	MPLSmode: fix(0), $			; Master Pattern Lambert Selector tag
        MCMPboth: long(0), $                    ; switch for MC or MC/MP display
        MCmin: long(0), $                       ; min value for MC display
        MCmax: long(0), $                       ; max value for MC display
        MPmin: float(0), $                      ; min value for MP display
        MPmax: float(0), $                      ; max value for MP display
        imageformat: long(0), $                 ; image format 'jpeg', 'tiff', 'bmp'
	logmode: fix(0), $			; keep a log file or not
	logunit: fix(13), $			; logical file unit for log file
	logname: '', $				; filename for log output
	logfileopen: fix(0), $			; log file is open when 1, closed when 0


	; widget location parameters
	xlocation: float(0.0), $		; main widget x-location (can be modified and stored in preferences file)
	ylocation: float(0.0), $		; main widget y-location (can be modified and stored in preferences file)
	EBSDxlocation: float(600.0), $		; image widget x-location (can be modified and stored in preferences file)
	EBSDylocation: float(100.0), $		; image widget y-location 
	MCxlocation: float(600.0), $		; Monte Carlo widget x-location (can be modified and stored in preferences file)
	MCylocation: float(100.0), $		; Monte Carlo widget y-location 
	MPxlocation: float(600.0), $		; Master Pattern widget x-location (can be modified and stored in preferences file)
	MPylocation: float(100.0), $		; Master Pattern widget y-location 
	Detectorxlocation: float(600.0), $	; detector widget x-location (can be modified and stored in preferences file)
	Detectorylocation: float(100.0), $	; detector widget y-location 
        scrdimx:0L, $                           ; display area x size in pixels 
        scrdimy:0L $                            ; display area y size in pixels 
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
EBSDdata.scrdimy = scr[1] * 0.8
EBSDdata.scrdimx = scr[0]
EBSDdata.xlocation = EBSDdata.scrdimx / 8.0
EBSDdata.ylocation = EBSDdata.scrdimx / 8.0

;------------------------------------------------------------
; does the preferences file exist ?  If not, create it, otherwise read it
EBSDgetpreferences,/noprint

;------------------------------------------------------------
; create the top level widget
EBSDwidget_s.base = WIDGET_BASE(TITLE='Electron Backscatter Diffraction Display Program', $
                        /ROW, $
                        XSIZE=1220, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDDisplay_event', $
                        XOFFSET=EBSDdata.xlocation, $
                        YOFFSET=EBSDdata.ylocation)

;------------------------------------------------------------
; create the two main columns
; block 1 is the left column, with logo and MC information
block1 = WIDGET_BASE(EBSDwidget_s.base, $
			/FRAME, $
			XSIZE=610, $
			/ALIGN_CENTER, $
			/COLUMN)

EBSDwidget_s.logodraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			/ALIGN_CENTER, $
			XSIZE=600, $
			YSIZE=200)

; this is the rightmost block; will be filled in later
block2 = WIDGET_BASE(EBSDwidget_s.base, $
			XSIZE=610, $
			/FRAME, $
			/COLUMN)



;------------------------------------------------------------
;------------------------------------------------------------
; this is the block that will display Monte Carlo file information; no user interactions here
; except for a single button.
block11 = WIDGET_BASE(block1, /FRAME, /COLUMN, TITLE='Monte Carlo Trajectory Simulation Parameters')

file1 = WIDGET_BASE(block11, /ROW, XSIZE=600, /ALIGN_RIGHT)
file2 = WIDGET_LABEL(file1, VALUE='Monte Carlo Trajectory Simulation Parameters', font=fontstrlarge, /ALIGN_RIGHT)

EBSDwidget_s.MCbutton= WIDGET_BUTTON(file1, $
                      UVALUE='MCDISPLAY', $
                      VALUE='Display', $
                      EVENT_PRO='EBSDDisplay_event', $
                      SENSITIVE=0, $
		      /ALIGN_RIGHT, $
                      /FRAME)

;---------- file name and size
file1 = WIDGET_BASE(block11, /ROW, XSIZE=600, /ALIGN_CENTER)
EBSDwidget_s.mcfilename = Core_WText(file1,'MC Data File Name', fontstrlarge, 200, 25, 60, 1, EBSDdata.mcfilename)
file1 = WIDGET_BASE(block11, /ROW, /BASE_ALIGN_BOTTOM, /ALIGN_LEFT)
EBSDwidget_s.mcfilesize = Core_WText(file1,'MC Data File Size', fontstrlarge, 200, 25, 30, 1, string(EBSDdata.mcfilesize,FORMAT="(I12)")+' bytes') 

;---------- min and max energy
file1 = WIDGET_BASE(block11, /ROW, XSIZE=600, /ALIGN_CENTER)
EBSDwidget_s.mcenergymin = Core_WText(file1,'Min/Max energy [keV]', fontstrlarge, 200, 25, 20, 1, string(EBSDdata.mcenergymin,format="(F6.2)"))
EBSDwidget_s.mcenergymax = Core_WText(file1,' / ', fontstrlarge, 25, 25, 20, 1, string(EBSDdata.mcenergymax,format="(F6.2)"))

;---------- bin size, number of bins
file1 = WIDGET_BASE(block11, /ROW, XSIZE=600, /ALIGN_CENTER)
EBSDwidget_s.mcenergybinsize = Core_WText(file1,'Energy binsize/#', fontstrlarge, 200, 25, 20, 1, string(EBSDdata.mcenergybinsize,format="(F6.2)"))
EBSDwidget_s.mcenergynumbin = Core_WText(file1,' / ', fontstrlarge, 25, 25, 20, 1, string(EBSDdata.mcenergynumbin,format="(I5)"))

;---------- depth max, and depth stepsize
file1 = WIDGET_BASE(block11, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcdepthmax = Core_WText(file1,'Depth max/step/#', fontstrlarge, 200, 25, 15, 1, string(EBSDdata.mcdepthmax,FORMAT="(F7.1)")) 
EBSDwidget_s.mcdepthstep = Core_WText(file1,' / ', fontstrlarge, 25, 25, 15, 1, string(EBSDdata.mcdepthstep,FORMAT="(F7.1)")) 
EBSDwidget_s.mcdepthnumbins = Core_WText(file1,' / ', fontstrlarge, 25, 25, 15, 1, string(EBSDdata.mcdepthnumbins,FORMAT="(I4)")) 

;---------- microscope voltage and Monte Carlo mode
file1 = WIDGET_BASE(block11, /ROW, /ALIGN_LEFT)
EBSDwidget_s.voltage = Core_WText(file1,'Voltage [kV]', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.voltage,FORMAT="(F7.1)")) 
EBSDwidget_s.mcmode = Core_WText(file1,'   MC Mode ', fontstrlarge, 120, 25, 10, 1, EBSDdata.mcmode) 

;---------- dimensions of modified Lambert projection 
file1 = WIDGET_BASE(block11, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcimx = Core_WText(file1,'Lambert Dimensions', fontstrlarge, 200, 25, 20, 1, string(EBSDdata.mcimx,FORMAT="(I5)")) 
EBSDwidget_s.mcimy = Core_WText(file1,'by', fontstr, 25, 25, 20, 1, string(EBSDdata.mcimy,FORMAT="(I5)")) 

;---------- total and backscattered numbers of electrons
file1 = WIDGET_BASE(block11, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mctotale = Core_WText(file1,'Total #/backscattered e', fontstrlarge, 200, 25, 20, 1, string(EBSDdata.mctotale,FORMAT="(I12)")) 
EBSDwidget_s.mcbse = Core_WText(file1,' / ', fontstrlarge, 25, 25, 20, 1, string(EBSDdata.mcbse,FORMAT="(I12)")) 

;---------- sample tilt angles 
file1 = WIDGET_BASE(block11, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcvangle = Core_WText(file1,'Sample tilt angles (v, h)', fontstrlarge, 200, 25, 20, 1, string(EBSDdata.mcvangle,FORMAT="(F7.2)")) 
EBSDwidget_s.mchangle = Core_WText(file1,', ', fontstrlarge, 25, 25, 20, 1, string(EBSDdata.mchangle,FORMAT="(F7.2)")) 

;------------------------------------------------------------
;------------------------------------------------------------




;------------------------------------------------------------
;------------------------------------------------------------
; this is the block that will display Master Pattern file information; no user interactions here
; except for a few buttons.
block21 = WIDGET_BASE(block2, /COLUMN, XSIZE=600, /ALIGN_RIGHT)
file2 = WIDGET_BASE(block21, /ROW, XSIZE=600, /ALIGN_RIGHT)
file3 = WIDGET_LABEL(file2, VALUE='Master Pattern Simulation Parameters', font=fontstrlarge, /ALIGN_CENTER)

EBSDwidget_s.MPbutton = WIDGET_BUTTON(file2, $
                      UVALUE='MPDISPLAY', $
                      VALUE='Display', $
                      EVENT_PRO='EBSDDisplay_event', $
                      SENSITIVE=0, $
		      /ALIGN_RIGHT, $
                      /FRAME)


;---------- file name and size
file1 = WIDGET_BASE(block21, /ROW, XSIZE=600, /ALIGN_CENTER)
EBSDwidget_s.mpfilename = Core_WText(file1,'MP Data File Name', fontstrlarge, 200, 25, 60, 1, EBSDdata.mpfilename)
file1 = WIDGET_BASE(block21, /ROW, /BASE_ALIGN_BOTTOM, /ALIGN_LEFT)
EBSDwidget_s.mpfilesize = Core_WText(file1,'MP Data File Size', fontstrlarge, 200, 25, 30, 1, string(EBSDdata.mpfilesize,FORMAT="(I12)")+' bytes') 

;---------- dimensions of modified Lambert projection 
file1 = WIDGET_BASE(block21, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mpimx = Core_WText(file1,'Lambert Dimensions', fontstrlarge, 200, 25, 20, 1, string(EBSDdata.mpimx,FORMAT="(I5)")) 
EBSDwidget_s.mpimy = Core_WText(file1,'by', fontstr, 25, 25, 20, 1, string(EBSDdata.mpimy,FORMAT="(I5)")) 

;---------- Lambert grid mode (hexagonal or square)
file1 = WIDGET_BASE(block21, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mpgridmode = Core_WText(file1,'Lambert grid', fontstrlarge, 200, 25, 20, 1, EBSDdata.mpgridmode)

;---------- crystal structure file name
file1 = WIDGET_BASE(block21, /ROW, /ALIGN_LEFT)
EBSDwidget_s.xtalname = Core_WText(file1,'Structure file', fontstrlarge, 200, 25, 20, 1, EBSDdata.xtalname)


;------------------------------------------------------------
;------------------------------------------------------------
; then we have the program message window

EBSDwidget_s.status= WIDGET_TEXT(block2, $
			XSIZE=90, $
			YSIZE=22, $
			/SCROLL, $
			VALUE=' ',$
			/ALIGN_LEFT)

; the following is needed by the Core_Print routine
status = EBSDwidget_s.status 

;------------------------------------------------------------
;------------------------------------------------------------
; finally we need a couple of control buttons: Quit, LoadMC, LoadMaster, Detector

file11 = WIDGET_BASE(block2, XSIZE=590, /FRAME, /ROW)

EBSDwidget_s.mainstop = WIDGET_BUTTON(file11, $
                                VALUE='Quit', $
                                UVALUE='QUIT', $
                                EVENT_PRO='EBSDDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

EBSDwidget_s.mcloadfile = WIDGET_BUTTON(file11, $
                                UVALUE='MCFILE', $
                                VALUE='Load MC file', $
                                EVENT_PRO='EBSDDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

EBSDwidget_s.mploadfile = WIDGET_BUTTON(file11, $
                                UVALUE='MPFILE', $
                                VALUE='Load master file', $
                                EVENT_PRO='EBSDDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

EBSDwidget_s.detector = WIDGET_BUTTON(file11, $
                                UVALUE='DETECTOR', $
                                VALUE='Define detector', $
                                EVENT_PRO='EBSDDisplay_event', $
                                SENSITIVE=0, $
                                /FRAME)

values = ['Off','On']
EBSDwidget_s.logfile= CW_BGROUP(file11, $
			values, $
			/FRAME, $
                        LABEL_LEFT='LogFile', $
			/ROW, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=EBSDdata.logmode, $
                        EVENT_FUNC='EBSDevent', $
			UVALUE='LOGFILE')

; the following is needed by the Core_Print routine
logmode = EBSDdata.logmode
logunit = EBSDdata.logunit

;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,EBSDwidget_s.base,/REALIZE

; realize the draw widgets
WIDGET_CONTROL, EBSDwidget_s.logodraw, GET_VALUE=drawID
EBSDwidget_s.logodrawID = drawID
;
read_jpeg,'Resources/SEMlogo.jpg',logo
wset,EBSDwidget_s.logodrawID
tvscl,logo,true=1

; and hand over control to the xmanager
XMANAGER,"EBSDDisplay",EBSDwidget_s.base,/NO_BLOCK

end

