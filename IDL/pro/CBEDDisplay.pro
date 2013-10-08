@CBEDprint			; appends messages to the status text widget
@CBEDgetpreferences		; read the preferences file
@CBEDwritepreferences		; write the preferences file
@CBEDGenerate2DSymmetry		; generate a set of 2D symmetry opertors for a given point group
@CBEDApply2DOperator		; apply a single symmetry operator to a disk image
@CBEDApply2DSymmetry		; apply all the symmetry operators to a disk image and add the results for Eades patterns
@CBEDCompute2DEquivalents 	; apply all the symmetry operators to a disk image and generate all equivalent disks
@CBEDDisplay_event		; main event handler 
@CBEDgetfilename		; select a data file
@CBEDreaddatafile		; read a data file
@CBEDprogressbar		; draws a progress bar during file loading

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
; CTEMsoft2013:CBEDDisplay.pro
;--------------------------------------------------------------------------
;
; PROGRAM: CBEDDisplay.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Zone axis CBED display, used for both MBCBED and LACBED programs
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
common PointGroups, PGTHD, PGTWD, DG

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
            dataname:long(0), $                 ; filename widget ID
            xtalname:long(0), $                 ; crystal structure filename widget ID
            filesize:long(0), $                 ; filesize widget ID
            pathname:long(0), $                 ; pathname widget ID
            loadlacbedfile:long(0), $           ; load LACBED file button ID
            loadmbcbedfile:long(0), $           ; load MBCBED file button ID
	    logfile: long(0), $			; logfile toggle widget ID
            imx:long(0), $                 	; number of image pixels along x
            imy:long(0), $                 	; number of image pixels along y
            patx:long(0), $                 	; number of CBED pattern pixels along x
            paty:long(0), $                 	; number of CBED pattern pixels along y
            CBEDzoom:long(0), $                 ; CBED pattern zoom factor
            maxHOLZ:long(0), $                  ; maximum HOLZ layer number
	    ga:long(0), $			; indices of horizontal g-vector
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
            numfam:long(0), $                   ; number of g-vectors
            numk:long(0), $                     ; number of k-vectors
            numt:long(0), $                     ; number of thickness values
            thetac:long(0), $                   ; beam divergence angle (mrad)
            aprad:long(0), $               	; aperture radius widget
            dfmode:long(0), $               	; set k/set g selection mode widget
	    symCPG: long(0), $			; crystal point group
	    symLPG: long(0), $			; Laue point group
	    symDPG: long(0), $			; diffraction point group
	    symPDG: long(0), $			; projection diffraction point group
	    symBFG: long(0), $			; Bright Field point group
	    symWPG: long(0), $			; Whole Pattern point group
	    symDFG: long(0), $			; General Dark Field point group
	    symDFS: long(0), $			; Special Dark Field point group
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
	scversion: '', $			; version identifier
	eventverbose: fix(0), $			; used for event debugging (0=off, 1=on)
	dataname: '', $				; filename (without pathname)
	pathname: '', $				; pathname (obviously)
	xtalname: '', $				; crystal structure filename (without pathname)
	suffix: '', $				; filename suffix 
	prefname: '~/.CBEDgui.prefs', $		; filename of preferences file (including path)
	filesize: long64(0), $			; input file size in bytes
	homefolder: '', $			; startup folder of the program
	CBEDroot: 'undefined', $		; current pathname (is stored in preferences file)
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
	ga: lonarr(3), $			; indices of horizontal g-vector
	addlog: float(0.0001), $		; factor to add for logarithmic CBED display
	imagelegend: long(0), $			; display image scale bar toggle (0=do not display, 1=display)
	imageformat: long(0), $			; image output format selector (0=jpeg, 1=tiff, 2=bmp)
	cbedlegend: long(0), $			; display cbed scale bar toggle (0=do not display, 1=display)
	cbedformat: long(0), $			; cbed output format selector (0=jpeg, 1=tiff, 2=bmp)
	cbedmode: long(0), $			; cbed intensity mode toggle (0=normal, 1=logarithmic)
	diffractionmode: long(0), $		; diffraction mode (MBCBED=0, LACBED=1)
	maxHOLZ: long(2), $			; maximum HOLZ lyer number in data file
	BFrho: float(3.5), $			; Bright Field radius in mrad
	Eadesrhoin: float(3.5), $		; Eades detector inner radius in mrad
	Eadesrhoout: float(10.0), $		; Eades detector outer radius in mrad
	BFmin: float(0.0), $			; BF image minimum intensity
	BFmax: float(0.0), $			; BF image maximum intensity
	Eadesmin: float(0.0), $			; Eades image minimum intensity
	Eadesmax: float(0.0), $			; Eades image maximum intensity
	CBEDmin: float(0.0), $			; CBED pattern minimum intensity
	CBEDmax: float(0.0), $			; CBED pattern maximum intensity
	patang: float(15.0), $			; this is the horizontal half scale of the detector plot in mm
        camlen: float(500), $                   ; camera length field (mm)
        refcamlen: float(1000), $               ; reference camera length to which the others will be scaled (mm)
	aprad: float(0.5), $			; aperture radius for CBED pattern computation in LACBED mode [mrad]
	apx: float(0.0), $			; aperture x position
	apy: float(0.0), $			; aperture y position
	apminrad: float(0.1), $			; minimal aperture radius for CBED imaging [mrad]
	rdisk: float(0.0), $			; disk radius in units of pixels
	wavek: lonarr(3), $			; wave vector indices
	wavelength: float(0.0), $		; wave length [nm] (will be displayed in [pm])
	dfl: float(1.0), $			; pixel size [nm]
	thetac: float(0.0), $			; beam convergence angle [mrad]
	nums: long(0), $			; number of pixels along disk radius (diameter = 2*nums+1)
	scale: float(0.0), $			; scale factor for CBED, [number of pixels per reciprocal nanometer]
	numfam: long(0), $			; number of reflections in CBED pattern
	numk: long(0), $			; number of wave vectors in CBED pattern
	numt: long(0), $			; number of sample thicknesses
	startthick: float(0), $			; starting thickness
	thickinc: float(0), $			; thickness increment
	symgroups: lonarr(8), $			; symmetry group labels
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
;CBEDgetpreferences

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

label2 = WIDGET_LABEL(file1, $
			VALUE='Data File Name', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.dataname = WIDGET_TEXT(file1, $
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
			VALUE='Disk Dimensions', $
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
			VALUE='# of g-families', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.numfam= WIDGET_TEXT(file5, $
			VALUE=string(data.numfam,format="(I5)"),$
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
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Maximum HOLZ    ', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.maxHOLZ= WIDGET_TEXT(file5, $
			VALUE=string(data.maxHOLZ,format="(I4)"),$
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
			VALUE='Zone axis [uvw]', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

wv = '['+string(data.wavek[0],format="(I2)")+' '+ string(data.wavek[1],format="(I2)")+' '+ string(data.wavek[2],format="(I2)")+']'
widget_s.wavek= WIDGET_TEXT(file7, $
			VALUE=wv,$
			XSIZE=20, $
			/ALIGN_LEFT)

;-------------
file7 = WIDGET_BASE(file6, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file7, $
			VALUE='Horizontal g   ', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

wv = '('+string(data.ga[0],format="(I2)")+' '+ string(data.ga[1],format="(I2)")+' '+ string(data.ga[2],format="(I2)")+')'
widget_s.ga= WIDGET_TEXT(file7, $
			VALUE=wv,$
			XSIZE=20, $
			/ALIGN_LEFT)


;------------------------------------------------------------
; block 3 displays a number of symmetry properties (there are 8 in total)
block3 = WIDGET_BASE(widget_s.base, $
			/FRAME, $
			/ROW)

file4 = WIDGET_BASE(block3, $
			/COLUMN, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Crystal PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

data.symgroups = [3,6,5,4,1,4,3,4]

widget_s.symCPG   = WIDGET_TEXT(file5, $
			VALUE=PGTHD[data.symgroups[0]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Diffraction PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symDPG   = WIDGET_TEXT(file5, $
			VALUE=DG[data.symgroups[2]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Whole Pattern PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symWPG   = WIDGET_TEXT(file5, $
			VALUE=PGTWD[data.symgroups[5]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Dark Field General PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symDFG   = WIDGET_TEXT(file5, $
			VALUE=PGTHD[data.symgroups[6]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
;-------------
file4 = WIDGET_BASE(block3, $
			/COLUMN, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Laue PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symLPG   = WIDGET_TEXT(file5, $
			VALUE=PGTHD[data.symgroups[1]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Projection Diff. PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symPDG   = WIDGET_TEXT(file5, $
			VALUE=DG[data.symgroups[3]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Bright Field PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symBFG   = WIDGET_TEXT(file5, $
			VALUE=PGTWD[data.symgroups[4]],$
			XSIZE=10, $
			/ALIGN_LEFT)

;-------------
file5 = WIDGET_BASE(file4, $
			/ROW, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file5, $
			VALUE='Dark Field Special PG', $
			FONT=fontstrlarge, $
			XSIZE=230, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.symDFS   = WIDGET_TEXT(file5, $
			VALUE=PGTHD[data.symgroups[7]],$
			XSIZE=10, $
			/ALIGN_LEFT)


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
                        EVENT_PRO='CBEDDisplay_event', $
                        SENSITIVE=1, $
                        /FRAME)

widget_s.loadlacbedfile = WIDGET_BUTTON(file11, $
                        VALUE='Load LACBED File', $
                        UVALUE='LOADLACBEDFILE', $
                        EVENT_PRO='CBEDDisplay_event', $
                        SENSITIVE=1, $
                        /FRAME)

widget_s.loadmbcbedfile = WIDGET_BUTTON(file11, $
                        VALUE='Load MBCBED File', $
                        UVALUE='LOADMBCBEDFILE', $
                        EVENT_PRO='CBEDDisplay_event', $
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
;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,widget_s.base,/REALIZE

; realize the draw widgets
;WIDGET_CONTROL, widget_s.detdraw, GET_VALUE=drawID
;widget_s.detdrawID = drawID
WIDGET_CONTROL, widget_s.progress, GET_VALUE=drawID
widget_s.progressdrawID = drawID

; and hand over control to the xmanager
XMANAGER,"CBEDDisplay",widget_s.base,/NO_BLOCK

; init the status text window
CBEDprint,'Zone Axis CBED Display Program [M. De Graef, 2013]',/blank
CBEDprint,'',/blank
CBEDprint,'Please load a file'


end

