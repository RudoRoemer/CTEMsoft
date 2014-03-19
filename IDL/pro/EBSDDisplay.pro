@EBSDDisplay_event    		; EBSD event handler
@EBSDgetfilename			; select a geometry file
@EBSDreaddatafile		; read geometry and data files
@EBSDatternWidget		; display widget
@EBSDatternWidget_event		; even handler
@EBSDshow			; show an EBSDattern (with grid)
@EBSDevent			; event handler for button groups
@EBSDgetpreferences		; read preferences file
@EBSDwritepreferences		; write preferences file
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
;> @todo 03/19/14 This is the very first implementation of the EBSD visualization GUI.
;> At this stage, we allow the user to display various versions of the Monte Carlo
;> results, as well as the master EBSD pattern in a number of projection modes.
;> The user can then define a series of detector parameters and compute EBSD
;> patterns.  In this version, there are no PSF, scintillator energy dependence,
;> or microscope geometric distortions; only the detector geometry, noise (on/off)
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
common EBSD_data_common, data
common fontstrings, fontstr, fontstrlarge, fontstrsmall
common PointGroups, PGTHD, PGTWD, DG

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
            base:long(0), $                     ; base widget ID

            EBSDatternbase:long(0), $            ; display window widget ID
            EBSDthicklist:long(0), $             ; integration depth widget
            EBSDgridbgroup:long(0), $            ; grid button group
            EBSDformatbgroup:long(0), $          ; file format button group
            EBSDDraw:long(0), $                  ; pattern draw widget
            EBSDDrawID:long(0), $                ; pattern draw widget ID
            filename:long(0), $                 ; file name field
            xtalname:long(0), $                 ; structure file name field
            filesize:long(0), $                 ; file size field
            mainstop:long(0), $                 ; program quit button
            loadfile:long(0), $                 ; load file button
            symCPG:long(0), $                   ; crystallographic point group
            symWPG:long(0), $                   ; whole pattern symmetry group
            imx:long(0), $                      ; pattern x dimension
            imy:long(0), $                      ; pattern y dimension
            numk:long(0), $                     ; number of wave vectors used 
            numthick:long(0), $                 ; number of thickness values
            thetac:long(0), $                   ; pattern convergence angle
            wavek:long(0), $                    ; zone axis indices 
            blur:float(0.0), $                  ; blur factor widget
            cx:long(0), $                       ; x-coordinate field
            cy:long(0), $                       ; y-coordinate field
            voltage:long(0), $                  ; microscope voltage
            abcdist:long(0), $                  ; lattice parameters
            albegadist:long(0), $               ; lattice parameters (angles)
            logodraw:long(0), $                 ; logo widget 
            logodrawID:long(0) $                ; logo widget ID
           }

EBSDdata = {EBSDdatastruct, $
	; Monte Carlo parameters first 
	mcdataname: '', $			; Monte Carlo result file name
	mcfilesize: long64(0), $		; Monte Carlo file size [bytes]
	mcenergymin: float(0.0), $		; minimum energy
	mcenergymax: float(0.0), $		; maximum energy
	mcenergybinsize: float(0.0), $		; energy bin size
	mcenergynumbins: long(0), $		; number of energy bins
	voltage: float(0.0), $			; microscope voltage
	mcdepthmax: float(0.0), $		; maximum depth in MC file
	mcdepthstep: float(0.0), $		; depth step size
	mcimx: long(0), $			; number of pixels along x in modified Lambert map
	mcimy: long(0), $			; same along y
	mctotale: long(0), $			; total number of electrons hitting the sample
	mcbse: long(0), $			; total number of BSE electrons
	mcvangle: float(0.0), $			; vertical sample tilt angle (around TD)
	mchangle: float(0.0), $			; horizontal sample tilt angle (around RD)
	mcmode: '', $				; 'CSDA' (continuous slowing down approximation) or 'DLOS' (discrete losses)
	; then Master Pattern parameters

	; then general program parameters
	eventverbose: fix(0), $			; used for event debugging (0=off, 1=on)
	dataname: '', $				; filename (without pathname)
	pathname: '', $				; pathname (obviously)
	suffix: '', $				; filename suffix 
	filesize: long64(0), $			; input file size in bytes
	homefolder: '', $			; startup folder of the program
	EBSDroot: 'undefined', $		; current pathname (is stored in preferences file)
	prefname: '~/.EBSDgui.prefs', $		; filename of preferences file (including path)
	nprefs: fix(0), $			; number of preferences in file
	progname: '',$				; program name
	scversion: '',$				; program version

	datadims: lon64arr(3), $		; dimensions of raw data array
	xtalname: '', $				; crystal structure filename
	wavelength: float(0.0), $		; electron wavelength
	ecplegend: long(0), $			; display pattern scale bar toggle (0=do not display, 1=display)
	ecpformat: long(0), $			; pattern output format selector (0=jpeg, 1=tiff, 2=bmp)
	ecpgrid: long(0), $			; grid display selector (0 = off, 1 = on)
	nums: long(0), $			; number of pixels along disk radius (diameter = 2*nums+1)
	scale: float(0.0), $			; scale factor for CBED, [number of pixels per reciprocal nanometer]

	; widget location parameters
	xlocation: float(0.0), $		; main widget x-location (can be modified and stored in preferences file)
	ylocation: float(0.0), $		; main widget y-location (can be modified and stored in preferences file)
	EBSDxlocation: float(600.0), $		; image widget x-location (can be modified and stored in preferences file)
	EBSDylocation: float(100.0), $		; image widget y-location 
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
if (sar.gt.(1.1*16.0/9.0) then begin
	scr[0] = scr[0]/2
end
EBSDdata.scrdimy = scr[1] * 0.8
EBSDdata.scrdimx = scr[0]
EBSDdata.xlocation = EBSDdata.scrdimx / 8.0
EBSDdata.ylocation = EBSDdata.scrdimx / 8.0

;------------------------------------------------------------
; does the preferences file exist ?  If not, create it, otherwise read it
EBSDgetpreferences

;------------------------------------------------------------
; create the top level widget
EBSDwidget_s.base = WIDGET_BASE(TITLE='Electron Backscatter Diffraction Pattern Display Program', $
                        /COLUMN, $
                        XSIZE=700, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDDisplay_event', $
                        XOFFSET=EBSDdata.xlocation, $
                        YOFFSET=EBSDdata.ylocation)

;------------------------------------------------------------
; create the various vertical blocks
; block 1 deals with the input file and displays the EBSDdata.dimensions
block1 = WIDGET_BASE(EBSDwidget_s.base, $
			/FRAME, $
			/ALIGN_CENTER, $
			/COLUMN)

EBSDwidget_s.logodraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			/ALIGN_CENTER, $
			XSIZE=600, $
			YSIZE=200)


;------------------------------------------------------------
;------------------------------------------------------------
; then we need a couple of control buttons: Quit, LoadMC, LoadMaster, Detector



;------------------------------------------------------------
;------------------------------------------------------------
; this is the middle block that will display Monte Carlo file information; no user interactions here
; except for a few buttons.
block1 = WIDGET_BASE(EBSDwidget_s.base, /FRAME, /COLUMN, TITLE='Monte Carlo Trajectory Simulation Parameters')

;----------
; things to display from the MC file:  
;   file name and size 
;   min and max energy, bins size, number of bins
;   microscope voltage
;   depth max, and depth step size
;   number of pixels in modified Lambert projection 
;   total number of electrons in MC run; total number of BSE electrons
;   sample tilt angles
;   MC mode (CSDA, other)


;---------- file name and size
file1 = WIDGET_BASE(block1, /ROW, XSIZE=700, /ALIGN_CENTER)
EBSDwidget_s.mcfilename = Core_WText(file1,'MC Data File Name', fontstrlarge, 200, 25, 77, 1, EBSDdata.mcdataname)
file1 = WIDGET_BASE(block1, /ROW, /BASE_ALIGN_BOTTOM, /ALIGN_LEFT)
EBSDwidget_s.mcfilesize = Core_WText(file1,'MC Data File Size', fontstrlarge, 200, 25, 40, 1, string(EBSDdata.mcfilesize,FORMAT="(I)")+' bytes',/aright) 

;---------- min and max energy
file1 = WIDGET_BASE(block1, /ROW, XSIZE=700, /ALIGN_CENTER)
EBSDwidget_s.mcenergymin = Core_WText(file1,'Min/Max energy', fontstrlarge, 200, 25, 77, 1, string(EBSDdata.mcenergymin,format="(F6.2)"))
EBSDwidget_s.mcenergymax = Core_WText(file1,' / ', fontstrlarge, 3, 25, 77, 1, string(EBSDdata.mcenergymax,format="(F6.2)"))

;---------- bin size, number of bins
file1 = WIDGET_BASE(block1, /ROW, XSIZE=700, /ALIGN_CENTER)
EBSDwidget_s.mcenergybinsize = Core_WText(file1,'Energy binsize/#', fontstrlarge, 200, 25, 77, 1, string(EBSDdata.mcenergybinsize,format="(F6.2)"))
EBSDwidget_s.mcenergynumbin = Core_WText(file1,' / ', fontstrlarge, 3, 25, 77, 1, string(EBSDdata.mcenergynumbin,format="(I5)"))

;---------- microscope voltage
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.voltage = Core_WText(file1,'Voltage [V]', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.voltage,FORMAT="(F7.1)")) 

;---------- depth max, and depth stepsize
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcdepthmax = Core_WText(file1,'Depth max/step [nm]', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mcdepthmax,FORMAT="(F7.1)")) 
EBSDwidget_s.mcdepthstep = Core_WText(file1,' / ', fontstrlarge, 3, 25, 10, 1, string(EBSDdata.mcdepthstep,FORMAT="(F7.1)")) 

;---------- dimensions of modified Lambert projection 
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcimx = Core_WText(file1,'Lambert Dimensions', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mcimx,FORMAT="(I5)")) 
EBSDwidget_s.mcimy = Core_WText(file1,'by', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mcimy,FORMAT="(I5)")) 

;---------- total and backscattered numbers of electrons
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mctotale = Core_WText(file1,'Total #/backscattered e', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mctotale,FORMAT="(I12)")) 
EBSDwidget_s.mcbse = Core_WText(file1,' / ', fontstrlarge, 3, 25, 10, 1, string(EBSDdata.mcbse,FORMAT="(I12)")) 

;---------- sample tilt angles 
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcvangle = Core_WText(file1,'Sample tilt angles (v, h)', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mcvangle,FORMAT="(F7.2)")) 
EBSDwidget_s.mchangle = Core_WText(file1,', ', fontstrlarge, 3, 25, 10, 1, string(EBSDdata.mhvangle,FORMAT="(F7.2)")) 

;---------- Monte Carlo mode 
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mcmode = Core_WText(file1,'MC Mode ', fontstrlarge, 200, 25, 10, 1, EBSDdata.mcmode) 

;---------- and finally we need a button on the lower right to call up the MCDisplay widget


;------------------------------------------------------------
;------------------------------------------------------------




;------------------------------------------------------------
;------------------------------------------------------------
; this is the central block that will display Master Pattern file information; no user interactions here
; except for a few buttons.
block1 = WIDGET_BASE(EBSDwidget_s.base, /FRAME, /COLUMN, TITLE='Master Pattern Simulation Parameters')

;----------
; EBSDdata.things to display from the Master Pattern file:  
;   file name and size
;   number of pixels in modified Lambert projection
;   xtal file name
;   space group
;   etc...

;---------- file name and size
file1 = WIDGET_BASE(block1, /ROW, XSIZE=700, /ALIGN_CENTER)
EBSDwidget_s.mpfilename = Core_WText(file1,'MP Data File Name', fontstrlarge, 200, 25, 77, 1, EBSDdata.mpdataname)
file1 = WIDGET_BASE(block1, /ROW, /BASE_ALIGN_BOTTOM, /ALIGN_LEFT)
EBSDwidget_s.mpfilesize = Core_WText(file1,'MP Data File Size', fontstrlarge, 200, 25, 40, 1, string(EBSDdata.mpfilesize,FORMAT="(I)")+' bytes',/aright) 

;---------- dimensions of modified Lambert projection 
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.mpimx = Core_WText(file1,'Lambert Dimensions', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mpimx,FORMAT="(I5)")) 
EBSDwidget_s.mpimy = Core_WText(file1,'by', fontstrlarge, 200, 25, 10, 1, string(EBSDdata.mpimy,FORMAT="(I5)")) 

;---------- crystal structure file name
file1 = WIDGET_BASE(block1, /ROW, /ALIGN_LEFT)
EBSDwidget_s.xtalname = Core_WText(file1,'Structure file', fontstrlarge, 200, 25, 20, 1, EBSDdata.xtalname)







;----------- next we have a series of parameters that are 
; derived from the input file and can not be changed by
; the user...

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='# of thicknesses', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

EBSDwidget_s.numthick= WIDGET_TEXT(file5, $
			VALUE=string(EBSDdata.datadims[2],format="(I5)"),$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Beam Convergence [mrad]', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

EBSDwidget_s.thetac= WIDGET_TEXT(file5, $
			VALUE=string(EBSDdata.thetac,format="(F6.3)"),$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------


;-------------
;-------------
file6 = WIDGET_BASE(block2, $
			/COLUMN, $
			/ALIGN_LEFT)

;-------------
file7 = WIDGET_BASE(file6, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file7, $
			VALUE='# of k-vectors', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

EBSDwidget_s.numk= WIDGET_TEXT(file7, $
			VALUE=string(EBSDdata.numk,format="(I5)"),$
			XSIZE=20, $
			/ALIGN_LEFT)

;-------------
file7 = WIDGET_BASE(file6, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file7, $
			VALUE='Structure File', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

EBSDwidget_s.xtalname= WIDGET_TEXT(file7, $
			VALUE=EBSDdata.xtalname,$
			XSIZE=20, $
			/ALIGN_LEFT)

;-------------
file7 = WIDGET_BASE(file6, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file7, $
			VALUE='Zone Axis [uvw]', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

wv = '['+string(EBSDdata.wavek[0],format="(I2)")+' '+ string(EBSDdata.wavek[1],format="(I2)")+' '+ string(EBSDdata.wavek[2],format="(I2)")+']'
EBSDwidget_s.wavek= WIDGET_TEXT(file7, $
			VALUE=wv,$
			XSIZE=20, $
			/ALIGN_LEFT)


;----------- next we have the lattice parameters

block2 = WIDGET_BASE(EBSDwidget_s.base, $
			/FRAME, $
			/COLUMN)

file4 = WIDGET_BASE(block2, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file4, $
			VALUE='Lattice parameters ', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

pp = string(EBSDdata.abcdist[0],format="(F10.5)")+', '+ string(EBSDdata.abcdist[1],format="(F10.5)")+', '+ string(EBSDdata.abcdist[2],format="(F10.5)")+' [nm]'
EBSDwidget_s.abcdist= WIDGET_TEXT(file4, $
			VALUE=pp,$
			XSIZE=45, $
			/ALIGN_LEFT)

file4 = WIDGET_BASE(block2, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file4, $
			VALUE='                   ', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

pp = string(EBSDdata.albegadist[0],format="(F10.5)")+', '+ string(EBSDdata.albegadist[1],format="(F10.5)")+', '+ string(EBSDdata.albegadist[2],format="(F10.5)")+' [degrees]'
EBSDwidget_s.albegadist= WIDGET_TEXT(file4, $
			VALUE=pp,$
			XSIZE=45, $
			/ALIGN_LEFT)

;------------------------------------------------------------
block3 = WIDGET_BASE(EBSDwidget_s.base, /FRAME, /ROW)

;-------------
file5 = WIDGET_BASE(block3, /ROW, /ALIGN_LEFT)
EBSDwidget_s.symCPG = Core_WText(file5, 'Crystal PG',fontstrlarge, 230, 25, 10, 1, PGTHD[EBSDdata.symgroups[0]] )
EBSDwidget_s.symWPG = Core_WText(file5, 'Whole Pattern PG',fontstrlarge, 230, 25, 10, 1, PGTWD[EBSDdata.symgroups[5]] )


;------------------------------------------------------------
; block 3 QUIT button, LOAD FILE button
block3 = WIDGET_BASE(EBSDwidget_s.base, $
			XSIZE=650, $
			/FRAME, $
			/ROW)

file11 = WIDGET_BASE(block3, $
			/ROW, $
			/ALIGN_LEFT)

EBSDwidget_s.mainstop = WIDGET_BUTTON(file11, $
                                VALUE='Quit', $
                                UVALUE='QUIT', $
                                EVENT_PRO='EBSDDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

EBSDwidget_s.loadfile = WIDGET_BUTTON(file11, $
                                VALUE='EBSD File', $
                                UVALUE='LOADFILE', $
                                EVENT_PRO='EBSDDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)


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

