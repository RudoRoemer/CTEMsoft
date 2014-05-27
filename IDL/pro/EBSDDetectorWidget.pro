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
; CTEMsoft2013:EBSDDetectorWidget.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDDetectorWidget.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Electron backscatter diffraction detector widget
;
;> @date 03/24/14 MDG 1.0 first attempt at detector widget
;--------------------------------------------------------------------------
pro EBSDDetectorWidget,dummy

;
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall
common PointGroups, PGTHD, PGTWD, DG

common projections, mcxcircle, mcycircle, mpxcircle, mpycircle, mcSPxcircle, mcSPycircle, mpSPxcircle, mpSPycircle 


;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("EBSDDetectorWidget") NE 0) then begin
  print,'EBSDDetectorWidget is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; create the top level widget
EBSDwidget_s.detectorbase = WIDGET_BASE(TITLE='Detector and Pattern Mode Widget', $
                        /ROW, $
                        XSIZE=810, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDDetectorWidget_event', $
                        XOFFSET=EBSDdata.Detectorxlocation, $
                        YOFFSET=EBSDdata.Detectorylocation)

;------------------------------------------------------------
; create the two main columns
block1 = WIDGET_BASE(EBSDwidget_s.detectorbase, $
			XSIZE=350, $
			/ALIGN_TOP, $
			/COLUMN)

block2 = WIDGET_BASE(EBSDwidget_s.detectorbase, $
			XSIZE=440, $
			/ALIGN_TOP, $
			/COLUMN)

;------------------------------------------------------------
;------------------------------------------------------------
file1 = WIDGET_BASE(block1, /COLUMN, /FRAME, YPAD=8, XSIZE=340, /ALIGN_LEFT)
file2 = WIDGET_LABEL(file1, VALUE='Detector & Microscope Geometry', font=fontstrlarge, /ALIGN_LEFT, /FRAME)

; we'll create two columns, one with all the detector geometry parameters,
; the other with pattern variables, including Euler angles, Euler angle
; convention, origin position, display mode (background only or full pattern),
; energy filter window, ...

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detL = Core_WTextE(file2,'Scintillator Distance [micron]', fontstr, 250, 25, 10, 1, string(EBSDdata.detL,format="(F9.2)"),'DETL','EBSDDetectorWidget_event')

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.dettheta = Core_WTextE(file2,'Detector Tilt Angle [deg]', fontstr, 250, 25, 10, 1, string(EBSDdata.dettheta,format="(F6.2)"),'DETTHETA','EBSDDetectorWidget_event')

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detdelta = Core_WTextE(file2,'Scintillator Pixel Size [micron]', fontstr, 250, 25, 10, 1, string(EBSDdata.detdelta,format="(F7.2)"),'DETDELTA','EBSDDetectorWidget_event')

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detnumsx = Core_WTextE(file2,'Number of pixels ', fontstr, 140, 25, 5, 1, string(EBSDdata.detnumsx,format="(I4)"),'DETNUMSX','EBSDDetectorWidget_event')
EBSDwidget_s.detnumsy = Core_WTextE(file2,' by ', fontstr, 30, 25, 5, 1, string(EBSDdata.detnumsy,format="(I4)"),'DETNUMSY','EBSDDetectorWidget_event')

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detxpc = Core_WTextE(file2,'PC [pixels] ', fontstr, 140, 25,  8, 1, string(EBSDdata.detxpc,format="(F7.2)"),'DETXPC','EBSDDetectorWidget_event')
EBSDwidget_s.detypc = Core_WTextE(file2,' , ', fontstr, 20, 25,  8, 1, string(EBSDdata.detypc,format="(F7.2)"),'DETYPC','EBSDDetectorWidget_event')

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detbeamcurrent = Core_WTextE(file2,'Beam current [nA]', fontstr, 140, 25, 10, 1, string(EBSDdata.detbeamcurrent,format="(F7.3)"),'DETBEAMCURRENT','EBSDDetectorWidget_event')

;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detdwelltime = Core_WTextE(file2,'Dwell Time [mu s] ', fontstr, 140, 25, 10, 1, string(EBSDdata.detdwelltime,format="(F7.3)"),'DETDWELLTIME','EBSDDetectorWidget_event')


file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
vals = ['TSL', 'HKL']
EBSDwidget_s.EulerConvention = CW_BGROUP(file2, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'Euler phi2 Convention', $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSDEULERCONVENTION', $
                        SET_VALUE=EBSDdata.EulerConvention)


;------------------------------------------------------------
;file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_CENTER)
;vals = ['Background Only','Full Pattern']
;EBSDwidget_s.BGmode = CW_BGROUP(file2, $
;                        vals, $
;                        /ROW, $
;                        /NO_RELEASE, $
;                        /EXCLUSIVE, $
;                        FONT=fontstr, $
;                        LABEL_LEFT = 'Pattern Mode', $
;                        EVENT_FUNC ='EBSDevent', $
;                        UVALUE='PATTERNMODE', $
;                        SET_VALUE=EBSDdata.BGmode)

;------------------------------------------------------------
; min and max energy settings for filtered imaging
file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_CENTER)
tvals = strarr(EBSDdata.mcenergynumbin)
for i=0,EBSDdata.mcenergynumbin-1 do begin
  th = EBSDdata.mcenergymin + float(i)*EBSDdata.mcenergybinsize
  tvals[i] = string(th,format="(F5.2)")
end

EBSDdata.Eminsel = tvals[0]

EBSDwidget_s.EBSDminenergylist = WIDGET_DROPLIST(file2, $
			EVENT_PRO='EBSDDetectorWidget_event', $
			VALUE=tvals,$
			UVALUE='EBSDMINENERGYLIST', $
			/ALIGN_LEFT)
WIDGET_CONTROL, set_droplist_select=EBSDdata.Eminsel, EBSDwidget_s.EBSDminenergylist

EBSDdata.Emaxsel = tvals[EBSDdata.mcenergynumbin-1]

file3 = WIDGET_LABEL(file2, VALUE='Max ', font=fontstr)
EBSDwidget_s.EBSDmaxenergylist = WIDGET_DROPLIST(file2, $
			EVENT_PRO='EBSDDetectorWidget_event', $
			VALUE=tvals,$
			UVALUE='EBSDMAXENERGYLIST', $
			/ALIGN_LEFT)
WIDGET_CONTROL, set_droplist_select=EBSDdata.Emaxsel, EBSDwidget_s.EBSDmaxenergylist


;------------------------------------------------------------
;------------------------------------------------------------

EBSDwidget_s.EBSDpatternfilename = Core_WText(file1,'EBSD Output File Name', fontstr, 200, 25, 50, 1, EBSDdata.EBSDpatternfilename)

EBSDwidget_s.EBSDgetpatternfilename = WIDGET_BUTTON(file1, $
                      UVALUE='GETEBSDFILENAME', $
                      VALUE='Set Output File Name', $
                      EVENT_PRO='EBSDDetectorWidget_event', $
                      SENSITIVE=1, $
		      /ALIGN_LEFT, $
                      /FRAME)

;------------------------------------------------------------
;------------------------------------------------------------
; and here is the Close button
file1 = WIDGET_BASE(block1, XSIZE=340, /ALIGN_LEFT, /ROW)

EBSDwidget_s.DetectorClose = WIDGET_BUTTON(file1, $
                                UVALUE='CLOSEDETECTOR', $
                                VALUE='Close', $
                                EVENT_PRO='EBSDDetectorWidget_event', $
                                SENSITIVE=1, $
                                /FRAME)

;------------------------------------------------------------
;------------------------------------------------------------
; this box defines the pattern mode and the output file name
file1 = WIDGET_BASE(block2, /COLUMN, /FRAME, YPAD=8, XSIZE=430, /ALIGN_LEFT)
file2 = WIDGET_LABEL(file1, VALUE='Pattern Mode', font=fontstrlarge, /ALIGN_LEFT, /FRAME)

vals = ['Single Pattern','Angle File','Dictionary']
EBSDwidget_s.Pmode = CW_BGROUP(file1, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='PMODE', $
                        SET_VALUE=EBSDdata.Pmode)

;------------------------------------------------------------
;------------------------------------------------------------
file1 = WIDGET_BASE(block2, /COLUMN, /FRAME, YPAD=8, XSIZE=430, /ALIGN_LEFT)
file2 = WIDGET_LABEL(file1, VALUE='Single Pattern Parameters', font=fontstrlarge, /ALIGN_LEFT,/FRAME)


;---------- 
file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_CENTER)
EBSDwidget_s.detphi1 = Core_WTextE(file2,'Euler [deg] phi1', fontstr, 120, 25, 8, 1, string(EBSDdata.detphi1,format="(F6.2)"),'DETphi1','EBSDDetectorWidget_event')
EBSDwidget_s.detphi = Core_WTextE(file2,' Phi', fontstr, 40, 25, 8, 1, string(EBSDdata.detphi,format="(F6.2)"),'DETPhi','EBSDDetectorWidget_event')
EBSDwidget_s.detphi2 = Core_WTextE(file2,' phi2', fontstr, 40, 25, 8, 1, string(EBSDdata.detphi2,format="(F6.2)"),'DETphi2','EBSDDetectorWidget_event')

file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_CENTER)
EBSDwidget_s.detax1 = Core_WTextE(file2,'axis', fontstr, 40, 25, 8, 1, string(EBSDdata.detax1,format="(F6.2)"),'DETax1','EBSDDetectorWidget_event')
EBSDwidget_s.detax2 = Core_WTextE(file2,'', fontstr, 4, 25, 8, 1, string(EBSDdata.detax2,format="(F6.2)"),'DETax2','EBSDDetectorWidget_event')
EBSDwidget_s.detax3 = Core_WTextE(file2,'', fontstr, 4, 25, 8, 1, string(EBSDdata.detax2,format="(F6.2)"),'DETax3','EBSDDetectorWidget_event')
EBSDwidget_s.detax4 = Core_WTextE(file2,'angle [deg]', fontstr, 85, 25, 8, 1, string(EBSDdata.detax4,format="(F6.2)"),'DETax4','EBSDDetectorWidget_event')

file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_CENTER)
EBSDwidget_s.DisplayEBSD = WIDGET_BUTTON(file2, $
                                VALUE='Display Pattern', $
                                UVALUE='DISPLAYEBSD', $
                                EVENT_PRO='EBSDDetectorWidget_event', $
				/ALIGN_LEFT, $
                                SENSITIVE=0, $
                                /FRAME)

vals = ['Off','On']
EBSDwidget_s.circularmask = CW_BGROUP(file2, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
			LABEL_LEFT='Circular Mask', $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='CIRCULARMASK', $
                        SET_VALUE=EBSDdata.showcircularmask)

;------------------------------------------------------------
;------------------------------------------------------------
file1 = WIDGET_BASE(block2, /COLUMN, /FRAME, YPAD=8, XSIZE=430, /ALIGN_LEFT)
file3 = WIDGET_BASE(file1, /ROW, /ALIGN_LEFT)
file2 = WIDGET_LABEL(file3, VALUE='Angle File Parameters', font=fontstrlarge, /ALIGN_LEFT,/FRAME)

EBSDwidget_s.GoAngle = WIDGET_BUTTON(file3, $
                                VALUE='Go', $
                                UVALUE='GOANGLE', $
                                EVENT_PRO='EBSDDetectorWidget_event', $
                                SENSITIVE=0, $
                                /FRAME)


EBSDwidget_s.EBSDanglefilename = Core_WText(file1,'Angle File Name', fontstr, 200, 25, 50, 1, EBSDdata.EBSDanglefilename)

EBSDwidget_s.EBSDgetanglefilename = WIDGET_BUTTON(file1, $
                      UVALUE='GETANGLEFILENAME', $
                      VALUE='Load Angle File', $
                      EVENT_PRO='EBSDDetectorWidget_event', $
                      SENSITIVE=0, $
		      /ALIGN_LEFT, $
                      /FRAME)

file2 = WIDGET_BASE(file1, /ROW, XSIZE=340, /ALIGN_LEFT)
EBSDwidget_s.angletype = Core_WText(file2,'Angle Type', fontstr, 90, 25, 10, 1, EBSDdata.angletype)
EBSDwidget_s.numangles = Core_WText(file2,'# Angles  ', fontstr, 80, 25, 10, 1, string(EBSDdata.numangles,format="(I8)"))

;------------------------------------------------------------
;------------------------------------------------------------
file1 = WIDGET_BASE(block2, /COLUMN, /FRAME, YPAD=8, XSIZE=430, /ALIGN_LEFT)
file3 = WIDGET_BASE(file1, /ROW, /ALIGN_LEFT)
file2 = WIDGET_LABEL(file3, VALUE='Dictionary Parameters', font=fontstrlarge, /ALIGN_LEFT,/FRAME)

EBSDwidget_s.GoDictionary = WIDGET_BUTTON(file3, $
                                VALUE='Go', $
                                UVALUE='GODICTIONARY', $
                                EVENT_PRO='EBSDetectorWidget_event', $
                                SENSITIVE=0, $
                                /FRAME)

PGs = [ 'Select Point Group', '1 (C1) - no symmetry' , '-1 (Ci) - no symmetry ', $ ; [triclinic]
 '2 (C2) - Cyclic', 'm (Cs) - Cyclic', '2/m (C2h) - Cyclic', $ ; [monoclinic]
 '222 (D2) - Dihedral', 'mm2 (C2v) - Dihedral', 'mmm (D2h) - Dihedral', $ ;[orthorhombic]
 '4 (C4) - Cyclic', '-4 (S4) - Cyclic', '4/m (C4h) - Cyclic', '422 (D4) - Dihedral', '4mm (C4v) - Dihedral', '-42m (D2d) - Dihedral', '4/mmm (D4h) - Dihedral', $ ; [tetragonal]
 '3 (C3) - Cyclic', '-3 (C3i) - Cyclic', '32 (D3) - Dihedral', '3m (C3v) - Dihedral', '-3m (D3d) - Dihedral', $ ; [trigonal]
 '6 (C6) - Cyclic', '-6 (C3h) - Cyclic', '6/m (C6h) - Cyclic', '622 (D6) - Dihedral', '6mm (C6v) - Dihedral', '-6m2 (D3h) - Dihedral', '6/mmm (D6h) - Dihedral', $ ; [hexagonal]
 '23 (T) - Tetrahedral', 'm3 (Th) - Tetrahedral', '432 (O) - Octahedral', '-43m (Td) - Tetrahedral', 'm-3m (Oh) - Octahedral'] ; [cubic]

; FZtype
; 0        no symmetry at all
; 1        cyclic symmetry
; 2        dihedral symmetry
; 3        tetrahedral symmetry
; 4        octahedral symmetry
 
; first entry is a space filler
FZtarray = [0, 0,0,1,1,1,2,2,2,1,1,1,2,2,2,2,1,1,2,2,2,1,1,1,2,2,2,2,3,3,4,3,4 ]
FZoarray = [0, 0,0,2,2,2,2,2,2,4,4,4,4,4,4,4,3,3,3,3,3,6,6,6,6,6,6,6,0,0,0,0,0 ]

file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_LEFT)
file3 = WIDGET_LABEL(file2, VALUE='Point Group:', font=fontstr)
EBSDwidget_s.PGdroplist = WIDGET_DROPLIST(file2, $
			EVENT_PRO='EBSDDetectorWidget_event', $
			VALUE=PGs,$
			SENSITIVE = 0, $
			UVALUE='DICTIONARYPG', $
			/ALIGN_LEFT)
WIDGET_CONTROL, set_droplist_select=EBSDdata.Dictpointgroup, EBSDwidget_s.PGdroplist

file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_LEFT)
EBSDwidget_s.Ncubochoric = Core_WTextE(file2,'N (# sampling points)', fontstr, 175, 25, 10, 1, string(EBSDdata.Ncubochoric,format="(I4)"),'NCUBOCHORIC','EBSDDetectorWidget_event')

;file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_LEFT)
EBSDwidget_s.EBSDdictfilename = Core_WText(file1,'Dictionary File Name', fontstr, 160, 25, 50, 1, EBSDdata.EBSDdictfilename)

EBSDwidget_s.EBSDgetdictfilename = WIDGET_BUTTON(file1, $
                      UVALUE='GETDICTFILENAME', $
                      VALUE='Set Dictionary File Name', $
                      EVENT_PRO='EBSDDetectorWidget_event', $
                      SENSITIVE=0, $
		      /ALIGN_LEFT, $
                      /FRAME)

file2 = WIDGET_BASE(file1, /ROW, XSIZE=430, /ALIGN_LEFT)
EBSDwidget_s.GoDict = WIDGET_BUTTON(file2, $
                                VALUE='Create Dictionary', $
                                UVALUE='GODICT', $
                                EVENT_PRO='EBSDDetectorWidget_event', $
                                SENSITIVE=0, $
                                /FRAME)

EBSDwidget_s.NinRFZ = Core_WText(file2,'# points in RFZ', fontstr, 125, 25, 10, 1, string(EBSDdata.NinRFZ,FORMAT="(I10)")) 


;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,EBSDwidget_s.detectorbase,/REALIZE

; and hand over control to the xmanager
XMANAGER,"EBSDDetectorWidget",EBSDwidget_s.detectorbase,/NO_BLOCK


end

