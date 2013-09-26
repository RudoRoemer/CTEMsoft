@STEMgetfilename		; select a geometry file
@STEMcomputeCBEDpatterns	; convert the rawdata array to actual CBED patterns
@STEMImageWidget		; image display widget
@STEMImageWidget_event		; image display widget event handler
@STEMCBEDWidget			; CBED widget 
@STEMCBEDWidget_event		; CBED widget event handler
@STEMprint			; appends messages to the status text widget
@STEMimagelegend		; add an image micron marker
@STEMshowCBED			; display a CBED pattern

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
; CTEMsoft2013:STEMDisplay.pro
;--------------------------------------------------------------------------
;
; PROGRAM: CBEDDisplay.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Zone axis CBED computation and display
;
;> @date 09/25/13 MDG 1.0 first attempt at a user-friendly interface
;--------------------------------------------------------------------------
pro CBEDDisplay,dummy
;
;------------------------------------------------------------
; common blocks
common CBED_widget_common, widget_s
common CBED_data_common, data
common fontstrings, fontstr, fontstrlarge, fontstrsmall

!EXCEPT=0

;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("CBEDDisplay") NE 0) then begin
  print,'CBEDDisplay is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; define a few structures (one for widgets, and one for data)
widget_s = {widgetstruct, $
            base:long(0), $                     ; base widget ID
            status:long(0), $                   ; status text widget ID
            progress:long(0), $                 ; progress status bar widget ID
            progressdrawID:long(0), $           ; progress status bar widget draw ID
            imagebase:long(0), $                ; base widget ID for BF/HAADF image display widget
            cbedbase:long(0), $                 ; base widget ID for CBED pattern display widget
            filename:long(0), $                 ; filename widget ID
            xtalname:long(0), $                 ; crystal structure filename widget ID
            geometryfilename:long(0), $         ; geometry filename widget ID
            filesize:long(0), $                 ; filesize widget ID
            pathname:long(0), $                 ; pathname widget ID
            loadfile:long(0), $                 ; pathname widget ID
	    logfile: long(0), $			; logfile toggle widget ID
            imx:long(0), $                 	; number of image pixels along x
            imy:long(0), $                 	; number of image pixels along y
            patx:long(0), $                 	; number of CBED pattern pixels along x
            paty:long(0), $                 	; number of CBED pattern pixels along y
            CBEDzoom:long(0), $                 ; CBED pattern zoom factor
	    imagelegendbgroup:long(0), $	; image legend button group
	    imageformatbgroup:long(0), $	; image format button group
	    cbedlegendbgroup:long(0), $		; cbed legend button group
	    cbedformatbgroup:long(0), $		; cbed format button group
	    cbedmodebgroup:long(0), $		; cbed intensity mode button group
	    saveimage:long(0), $		; save image button
	    savecbed:long(0), $			; save cbed pattern button
	    bfrho: long(0), $			; Bright Field radius in mrad
	    haadfrhoin: long(0), $		; Dark Field inner radius in mrad
	    haadfrhoout: long(0), $		; Dark Field outer radius in mrad
	    detsegm: long(0), $			; number of HAADF detector segments
	    angsegm: long(0), $			; off set angle for first HAADF detector segment
	    patang: long(0), $			; maximum CBED pattern angle (from input file)
            detdraw:long(0), $                  ; detector draw widget ID
            detdrawID:long(0), $                ; detector draw window ID
	    BFdraw: long(0), $			; BF widget
	    BFdrawID: long(0), $		; BF widget ID
	    HAADFdraw: long(0), $		; HAADF widget
	    HAADFdrawID: long(0), $		; HAADF widget ID
	    CBEDdraw: long(0), $		; CBED widget
	    CBEDdrawID: long(0), $		; CBED widget ID
	    BFmin: long(0), $			; BF image minimum intensity
	    BFmax: long(0), $			; BF image maximum intensity
	    HAADFmin: long(0), $		; HAADF image minimum intensity
	    HAADFmax: long(0), $		; HAADF image maximum intensity
	    CBEDmin: long(0), $			; CBED pattern minimum intensity
	    CBEDmax: long(0), $			; CBED atptern maximum intensity
            camlen:long(0), $                   ; camera length field (mm)
            wavelength:long(0), $               ; wave length field (pm)
            wavek:long(0), $                    ; wave vector indices
            numref:long(0), $                   ; number of g-vectors
            numk:long(0), $                     ; number of k-vectors
            thetac:long(0), $                   ; beam divergence angle (mrad)
            aprad:long(0), $               	; aperture radius widget
            dfmode:long(0), $               	; set k/set g selection mode widget
            sectormode:long(0), $               ; single or multiple sector mode
            diffractionmode:long(0), $          ; HAADF or regular dark field diffraction mode
            gosector:long(0), $                 ; compute image for multiple sector mode
            clearsector:long(0), $              ; clear selected sectors and image
            inputbase:long(0), $                ; input base widget ID
            closecbed:long(0), $                ; close CBED widget button
            closeimage:long(0), $               ; close image widget button
            mainstop:long(0), $                 ; stop button
            topbgroup:long(0) $                 ; top level button group widget
           }

data = {datastruct, $
	eventverbose: fix(0), $			; used for event debugging (0=off, 1=on)
	geomname: '', $				; filename (without pathname) for the geometry file
	dataname: '', $				; filename (without pathname)
	pathname: '', $				; pathname (obviously)
	xtalname: '', $				; crystal structure filename (without pathname)
	suffix: '', $				; filename suffix 
	prefname: '~/.STEMgui.prefs', $		; filename of preferences file (including path)
	filesize: long64(0), $			; input file size in bytes
	homefolder: '', $			; startup folder of the program
	STEMroot: 'undefined', $		; current pathname (is stored in preferences file)
	nprefs: fix(0), $			; number of preferences in file
        status:'waiting for input', $           ; current status line
	logmode: fix(0), $			; keep a log file or not
	logunit: fix(13), $			; logical file unit for log file
	logname: '', $				; filename for log output
	logfileopen: fix(0), $			; log file is open when 1, closed when 0
	detwinx: fix(401), $			; detector display window size x
	detwiny: fix(401), $			; detector display window size y
	imx: long(0), $				; number of image pixels along x
	imy: long(0), $				; number of image pixels along y
	patx: long(900), $			; number of pattern pixels along x
	paty: long(900), $			; number of pattern pixels along y
	CBEDzoom: fix(1), $			; zoom factor for CBED pattern
	addlog: float(0.0001), $		; factor to add for logarithmic CBED display
	imagelegend: long(0), $			; display image scale bar toggle (0=do not display, 1=display)
	imageformat: long(0), $			; image output format selector (0=jpeg, 1=tiff, 2=bmp)
	cbedlegend: long(0), $			; display cbed scale bar toggle (0=do not display, 1=display)
	cbedformat: long(0), $			; cbed output format selector (0=jpeg, 1=tiff, 2=bmp)
	cbedmode: long(0), $			; cbed intensity mode toggle (0=normal, 1=logarithmic)
	diffractionmode: long(0), $		; diffraction mode (HAADF=0, regular dark field=1)
	BFrho: float(3.5), $			; Bright Field radius in mm
	HAADFrhoin: float(3.5), $		; Dark Field inner radius in mm
	HAADFrhoout: float(10.0), $		; Dark Field outer radius in mm
	BFmrad: float(0.0), $			; scaled BF radius in mrad
	HAADFimrad: float(0.0), $		; scaled HAADF inner radius in mrad
	HAADFomrad: float(0.0), $		; scaled HAADF outer radius in mrad
	BFmin: float(0.0), $			; BF image minimum intensity
	BFmax: float(0.0), $			; BF image maximum intensity
	HAADFmin: float(0.0), $			; HAADF image minimum intensity
	HAADFmax: float(0.0), $			; HAADF image maximum intensity
	CBEDmin: float(0.0), $			; CBED pattern minimum intensity
	CBEDmax: float(0.0), $			; CBED pattern maximum intensity
	detsegm: fix(1), $			; number of HAADF detector segments
	angsegm: float(0.0), $			; off set angle for first HAADF detector segment
	patang: float(15.0), $			; this is the horizontal half scale of the detector plot in mm
        camlen: float(500), $                   ; camera length field (mm)
        refcamlen: float(1000), $               ; reference camera length to which the others will be scaled (mm)
        sectormode: fix(0), $                   ; single or multiple sector imaging mode
	dfmode: fix(0), $			; set k (0) or set g (1) selection mode
	aprad: float(0.5), $			; aperture radius for DF imaging [mrad]
	apx: float(0.0), $			; aperture x position
	apy: float(0.0), $			; aperture y position
	apminrad: float(0.1), $			; minimal aperture radius for DF imaging [mrad]
	rdisk: float(0.0), $			; disk radius in units of pixels
	wavek: lonarr(3), $			; wave vector indices
	wavelength: float(0.0), $		; wave length [nm] (will be displayed in [pm])
	dfl: float(1.0), $			; pixel size [nm]
	thetac: float(0.0), $			; beam divergence angle [mrad]
	nums: long(0), $			; number of pixels along disk radius (diameter = 2*nums+1)
	bragg: float(0.0), $			; Bragg angle for ga reflection (rad)
	scale: float(0.0), $			; scale factor for CBED, [number of pixels per reciprocal nanometer]
	numref: long(0), $			; number of reflections in CBED pattern
	numk: long(0), $			; number of wave vectors in CBED pattern
	datadims: lon64arr(4), $		; dimensions of rawdata array
	xlocation: float(0.0), $		; main widget x-location (can be modified and stored in preferences file)
	ylocation: float(0.0), $		; main widget y-location (can be modified and stored in preferences file)
	imagexlocation: float(600.0), $		; image widget x-location (can be modified and stored in preferences file)
	imageylocation: float(100.0), $		; image widget y-location 
	cbedxlocation: float(1200.0), $		; cbed widget x-location (can be modified and stored in preferences file)
	cbedylocation: float(100.0), $		; cbed widget y-location 
        scrdimx:0L, $                           ; display area x size in pixels 
        scrdimy:0L $                            ; display area y size in pixels 
        }


; a few font strings
fontstr='-adobe-new century schoolbook-bold-r-normal--14-100-100-100-p-87-iso8859-1'
fontstrlarge='-adobe-new century schoolbook-medium-r-normal--20-140-100-100-p-103-iso8859-1'
fontstrsmall='-adobe-new century schoolbook-medium-r-normal--14-100-100-100-p-82-iso8859-1'

;------------------------------------------------------------
; get the display window size to 80% of the current screen size (but be careful with double screens ... )
device,decomposed = 0
device, GET_SCREEN_SIZE = scr
data.scrdimy = scr[1] * 0.8
data.scrdimx = 0.75 * data.scrdimy   ; doing it this way avoids problems with multiple screens
data.xlocation = data.scrdimx / 8.0
data.ylocation = data.scrdimx / 8.0

;------------------------------------------------------------
; does the preferences file exist ?  If not, create it, otherwise read it
CBEDgetpreferences

;------------------------------------------------------------
; create the top level widget
widget_s.base = WIDGET_BASE(TITLE='Zone Axis CBED Display Program', $
                        /COLUMN, $
                        XSIZE=700, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='CBEDDisplay_event', $
                        XOFFSET=data.xlocation, $
                        YOFFSET=data.ylocation)

;------------------------------------------------------------
; create the various vertical blocks
; block 1 deals with the input file and displays the data dimensions
block1 = WIDGET_BASE(widget_s.base, $
			/FRAME, $
			/COLUMN)

;----------
file1 = WIDGET_BASE(block1, $
			/ROW, $
                        XSIZE=980, $
			/ALIGN_CENTER)

label1 = WIDGET_LABEL(file1, $
			VALUE='Geometry File Name', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.geometryfilename = WIDGET_TEXT(file1, $
			VALUE=data.geomname,$
			XSIZE=77, $
			/ALIGN_LEFT)

file1 = WIDGET_BASE(block1, $
			/ROW, $
                        XSIZE=980, $
			/ALIGN_CENTER)

label2 = WIDGET_LABEL(file1, $
			VALUE='Data File Name', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.filename = WIDGET_TEXT(file1, $
			VALUE=data.dataname,$
			XSIZE=77, $
			/ALIGN_LEFT)

;----------
file3 = WIDGET_BASE(block1, $
			/ROW, $
			/BASE_ALIGN_BOTTOM, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file3, $
			VALUE='Data File Size', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.filesize = WIDGET_TEXT(file3, $
			VALUE=string(data.filesize,FORMAT="(I)")+' bytes', $
			XSIZE=40, $
			/ALIGN_RIGHT)

widget_s.progress = WIDGET_DRAW(file3, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/ALIGN_RIGHT, $
			XSIZE=200, $
			YSIZE=20)


file3 = WIDGET_BASE(block1, $
			/ROW, $
			/ALIGN_LEFT)

label4 = WIDGET_LABEL(file3, $
			VALUE='Image Dimensions', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.imx= WIDGET_TEXT(file3, $
			VALUE=string(data.imx,format="(I5)"),$
			XSIZE=10, $
			/ALIGN_LEFT)

labela = WIDGET_LABEL(file3, $
			VALUE='by', $
			FONT=fontstrlarge, $
			XSIZE=25, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.imy= WIDGET_TEXT(file3, $
			VALUE=string(data.imy,format="(I5)"),$
			XSIZE=10, $
			/ALIGN_LEFT)

;----------- next we have a series of parameters that are 
; derived from the input file and can not be changed by
; the user...

block2 = WIDGET_BASE(widget_s.base, $
			/FRAME, $
			/ROW)

file4 = WIDGET_BASE(block2, $
			/COLUMN, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='# of g-vectors', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.numref= WIDGET_TEXT(file5, $
			VALUE=string(data.numref,format="(I5)"),$
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

widget_s.thetac= WIDGET_TEXT(file5, $
			VALUE=string(data.thetac,format="(F6.3)"),$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Wave Length [pm]', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.wavelength= WIDGET_TEXT(file5, $
			VALUE=string(data.wavelength,format="(F7.4)"),$
			XSIZE=10, $
			/ALIGN_LEFT)



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

widget_s.numk= WIDGET_TEXT(file7, $
			VALUE=string(data.numk,format="(I5)"),$
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

widget_s.xtalname= WIDGET_TEXT(file7, $
			VALUE=data.xtalname,$
			XSIZE=20, $
			/ALIGN_LEFT)

;-------------
file7 = WIDGET_BASE(file6, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file7, $
			VALUE='Wave Vector [uvw]', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

wv = '['+string(data.wavek[0],format="(I2)")+' '+ string(data.wavek[1],format="(I2)")+' '+ string(data.wavek[2],format="(I2)")+']'
widget_s.wavek= WIDGET_TEXT(file7, $
			VALUE=wv,$
			XSIZE=20, $
			/ALIGN_LEFT)



;------------------------------------------------------------
; block 2 controls the detector dimensions and segment number+orientation
block2 = WIDGET_BASE(widget_s.base, $
			/FRAME, $
			/ROW)

cols = WIDGET_BASE(block2, $
			/COLUMN, $
			/ALIGN_LEFT)

cols1 = WIDGET_BASE(cols, $
			/COLUMN, $
			/FRAME, $
			/ALIGN_LEFT)

file5 = WIDGET_BASE(cols1, $
			/ROW, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(file5, $
			VALUE='BF radius (mm)', $
			FONT=fontstrlarge, $
			XSIZE=185, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.BFrho= WIDGET_TEXT(file5, $
			VALUE=string(data.BFrho,format="(F6.2)"),$
			XSIZE=10, $
			/EDITABLE, $
                        EVENT_PRO='STEMDisplay_event', $
			UVALUE='BFRHO', $
			/ALIGN_LEFT)

file5a = WIDGET_BASE(cols1, $
			/ROW, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(file5a, $
			VALUE='HAADF inner radius', $
			FONT=fontstrlarge, $
			XSIZE=185, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.HAADFrhoin= WIDGET_TEXT(file5a, $
			VALUE=string(data.HAADFrhoin,format="(F6.2)"),$
			XSIZE=10, $
			/EDITABLE, $
                        EVENT_PRO='STEMDisplay_event', $
			UVALUE='HAADFRHOIN', $
			/ALIGN_LEFT)

file5b = WIDGET_BASE(cols1, $
			/ROW, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(file5b, $
			VALUE='HAADF outer radius', $
			FONT=fontstrlarge, $
			XSIZE=185, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.HAADFrhoout= WIDGET_TEXT(file5b, $
			VALUE=string(data.HAADFrhoout,format="(F6.2)"),$
			XSIZE=10, $
			/EDITABLE, $
                        EVENT_PRO='STEMDisplay_event', $
			UVALUE='HAADFRHOOUT', $
			/ALIGN_LEFT)

;----------
file7 = WIDGET_BASE(cols1, $
			/ROW, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(file7, $
			VALUE='# detector segments', $
			FONT=fontstrlarge, $
			XSIZE=185, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.detsegm= WIDGET_TEXT(file7, $
			VALUE=string(data.detsegm,format="(I4)"),$
			XSIZE=10, $
			/EDITABLE, $
                        EVENT_PRO='STEMDisplay_event', $
			UVALUE='DETSEGM', $
			/ALIGN_LEFT)

;----------
file8 = WIDGET_BASE(cols1, $
			/ROW, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(file8, $
			VALUE='offset angle (deg) ', $
			FONT=fontstrlarge, $
			XSIZE=185, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.angsegm= WIDGET_TEXT(file8, $
			VALUE=string(data.angsegm,format="(F6.2)"),$
			XSIZE=10, $
			/EDITABLE, $
                        EVENT_PRO='STEMDisplay_event', $
			UVALUE='ANGSEGM', $
			/ALIGN_LEFT)

;----------
file9 = WIDGET_BASE(cols1, $
			/ROW, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(file9, $
			VALUE='camera length (mm) ', $
			FONT=fontstrlarge, $
			XSIZE=185, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.camlen= WIDGET_TEXT(file9, $
			VALUE=string(data.camlen,format="(F8.2)"),$
			XSIZE=10, $
			/EDITABLE, $
                        EVENT_PRO='STEMDisplay_event', $
			UVALUE='CAMLEN', $
			/ALIGN_LEFT)

;----------
cols2 = WIDGET_BASE(cols, $
			/COLUMN, $
			/FRAME, $
			XSIZE=270, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(cols2, $
			VALUE='Sector Selection Mode', $
			FONT=fontstrlarge, $
			XSIZE=250, $
			YSIZE=25, $
			/ALIGN_LEFT)

file10 = WIDGET_BASE(cols2, $
			/ROW, $
			/ALIGN_LEFT)

values = ['Single','Multiple']
widget_s.sectormode= CW_BGROUP(file10, $
			values, $
			/ROW, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=data.sectormode, $
                        EVENT_FUNC='STEMevent', $
			UVALUE='SECTORMODE')

widget_s.gosector = WIDGET_BUTTON(file10, $
			VALUE='Go', $
                        EVENT_PRO='STEMDisplay_event', $
			SENSITIVE = 1, $
			UVALUE='GOSECTOR', $
			/ALIGN_CENTER)

widget_s.clearsector = WIDGET_BUTTON(file10, $
			VALUE='Clear', $
                        EVENT_PRO='STEMDisplay_event', $
			SENSITIVE = 1, $
			UVALUE='CLEARSECTOR', $
			/ALIGN_CENTER)

;----------
cols2 = WIDGET_BASE(cols, $
			/COLUMN, $
			XSIZE=270, $
			/FRAME, $
			/ALIGN_LEFT)

label1 = WIDGET_LABEL(cols2, $
			VALUE='Diffraction Display Mode', $
			FONT=fontstrlarge, $
			XSIZE=250, $
			YSIZE=25, $
			/ALIGN_LEFT)

file10 = WIDGET_BASE(cols2, $
			/ROW, $
			/ALIGN_LEFT)

values = ['HAADF STEM','Regular Dark Field']
widget_s.diffractionmode= CW_BGROUP(file10, $
			values, $
			/ROW, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=data.diffractionmode, $
                        EVENT_FUNC='STEMevent', $
			UVALUE='DIFFRACTIONMODE')

;----------
; next, add the display window that will show the detector geometry 
; in units of the CBED pattern size

widget_s.detdraw = WIDGET_DRAW(block2, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			/BUTTON_EVENTS, $
			EVENT_PRO='STEMDisplay_event', $
			UVALUE = 'SELECTSECTOR', $
			TOOLTIP='Select one or more sectors, then hit GO to display BF/HAADF image', $
;		TOOLTIP='Click on location of aperture center in BF disk', $
;		TOOLTIP='Select a diffracted disk to display dark field image', $
			XSIZE=data.detwinx, $
			YSIZE=data.detwiny)


;----------
; for diffraction mode equal to regular dark field, these are the controls

file1 = WIDGET_BASE(widget_s.base, $
			/ROW, $
			/FRAME, $
			/ALIGN_LEFT)

label5 = WIDGET_LABEL(file1, $
			VALUE='Aperture radius [mrad]', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.aprad= WIDGET_TEXT(file1, $
			VALUE=string(data.aprad,FORMAT="(F6.2)"),$
			XSIZE=20, $
			SENSITIVE=0, $
			/EDITABLE, $
			EVENT_PRO='STEMDisplay_event', $
			UVALUE='APRAD',$
			/ALIGN_LEFT)

values = ['set k vector','select g']
widget_s.dfmode= CW_BGROUP(file1, $
			values, $
			/FRAME, $
                        LABEL_LEFT='mouse click mode', $
			/ROW, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=data.dfmode, $
                        EVENT_FUNC='STEMevent', $
			UVALUE='DFMODE')


;----------
; next, add a text window for program messages

widget_s.status= WIDGET_TEXT(widget_s.base, $
			XSIZE=115, $
			YSIZE=10, $
			/SCROLL, $
			VALUE=' ',$
			/ALIGN_CENTER)

;------------------------------------------------------------
; block 3 QUIT button, LOAD FILE button and progress bar (used for file loading)
block3 = WIDGET_BASE(widget_s.base, $
			XSIZE=650, $
			/FRAME, $
			/ROW)

file11 = WIDGET_BASE(block3, $
			/ROW, $
			/ALIGN_LEFT)

widget_s.mainstop = WIDGET_BUTTON(file11, $
                                VALUE='Quit', $
                                UVALUE='QUIT', $
                                EVENT_PRO='STEMDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

widget_s.loadfile = WIDGET_BUTTON(file11, $
                                VALUE='Load File', $
                                UVALUE='LOADFILE', $
                                EVENT_PRO='STEMDisplay_event', $
                                SENSITIVE=1, $
                                /FRAME)

values = ['Off','On']
widget_s.logfile= CW_BGROUP(file11, $
			values, $
			/FRAME, $
                        LABEL_LEFT='LogFile', $
			/ROW, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=data.logmode, $
                        EVENT_FUNC='STEMevent', $
			UVALUE='LOGFILE')


;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,widget_s.base,/REALIZE

; realize the draw widgets
WIDGET_CONTROL, widget_s.detdraw, GET_VALUE=drawID
widget_s.detdrawID = drawID
WIDGET_CONTROL, widget_s.progress, GET_VALUE=drawID
widget_s.progressdrawID = drawID

; and hand over control to the xmanager
XMANAGER,"CBEDDisplay",widget_s.base,/NO_BLOCK

; init the status text window
CBEDprint,'Zone Axis CBED Display Program [M. De Graef, 2013]',/blank
CBEDprint,'',/blank
CBEDprint,'Please select a display mode and then load a file',/blank


end

