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

common projections, xcircle, ycircle

;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("EBSDDetectorWidget") NE 0) then begin
  print,'EBSDDetectorWidget is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; create the top level widget
EBSDwidget_s.detectorbase = WIDGET_BASE(TITLE='Detector Widget', $
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
			/FRAME, $
			XSIZE=350, $
			/ALIGN_TOP, $
			/COLUMN)

block2 = WIDGET_BASE(EBSDwidget_s.detectorbase, $
			/FRAME, $
			XSIZE=450, $
			/ALIGN_TOP, $
			/COLUMN)

;------------------------------------------------------------
;------------------------------------------------------------
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_RIGHT)
file2 = WIDGET_LABEL(file1, VALUE='Detector & Microscope Geometry', font=fontstrlarge, /ALIGN_RIGHT, /FRAME)

; we'll create two columns, one with all the detector geometry parameters,
; the other with pattern variables, including Euler angles, Euler angle
; convention, origin position, display mode (background only or full pattern),
; energy filter window, ...

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detL = Core_WTextE(file1,'Scintillator Distance [micron]', fontstr, 250, 25, 10, 1, string(EBSDdata.detL,format="(F9.2)"),'DETL','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.dettheta = Core_WTextE(file1,'Detector Tilt Ange [deg]', fontstr, 250, 25, 10, 1, string(EBSDdata.dettheta,format="(F6.2)"),'DETTHETA','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detdelta = Core_WTextE(file1,'Scintillator Pixel Size [micron]', fontstr, 250, 25, 10, 1, string(EBSDdata.detdelta,format="(F7.2)"),'DETDELTA','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detnumsx = Core_WTextE(file1,'Number of pixels ', fontstr, 140, 25, 5, 1, string(EBSDdata.detnumsx,format="(I4)"),'DETNUMSX','EBSDDetectorWidget_event')
EBSDwidget_s.detnumsy = Core_WTextE(file1,' by ', fontstr, 30, 25, 5, 1, string(EBSDdata.detnumsy,format="(I4)"),'DETNUMSY','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detxpc = Core_WTextE(file1,'PC [pixels] ', fontstr, 140, 25,  8, 1, string(EBSDdata.detxpc,format="(F7.2)"),'DETXPC','EBSDDetectorWidget_event')
EBSDwidget_s.detypc = Core_WTextE(file1,' , ', fontstr, 20, 25,  8, 1, string(EBSDdata.detxpc,format="(F7.2)"),'DETXPY','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detbinning = Core_WTextE(file1,'Detector Binning', fontstr, 140, 25, 5, 1, string(EBSDdata.detbinning,format="(I3)"),'DETBINNING','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detbeamcurrent = Core_WTextE(file1,'Beam current [A]', fontstr, 140, 25, 10, 1, string(EBSDdata.detbeamcurrent,format="(D9.2)"),'DETBEAMCURRENT','EBSDDetectorWidget_event')

;---------- 
file1 = WIDGET_BASE(block1, /ROW, XSIZE=340, /ALIGN_CENTER)
EBSDwidget_s.detdwelltime = Core_WTextE(file1,'Dwell Time [s] ', fontstr, 140, 25, 10, 1, string(EBSDdata.detdwelltime,format="(D9.2)"),'DETDWELLTIME','EBSDDetectorWidget_event')


;------------------------------------------------------------
;------------------------------------------------------------


;------------------------------------------------------------
;------------------------------------------------------------
file1 = WIDGET_BASE(block2, /ROW, XSIZE=440, /ALIGN_RIGHT)
file2 = WIDGET_LABEL(file1, VALUE='Pattern Parameters', font=fontstrlarge, /ALIGN_RIGHT,/FRAME)

;---------- 
file1 = WIDGET_BASE(block2, /ROW, XSIZE=440, /ALIGN_CENTER)
EBSDwidget_s.detphi1 = Core_WTextE(file1,'Euler [deg] phi1', fontstr, 120, 25, 8, 1, string(EBSDdata.detphi1,format="(F6.2)"),'DETphi1','EBSDDetectorWidget_event')
EBSDwidget_s.detphi = Core_WTextE(file1,' Phi', fontstr, 40, 25, 8, 1, string(EBSDdata.detphi,format="(F6.2)"),'DETPhi','EBSDDetectorWidget_event')
EBSDwidget_s.detphi2 = Core_WTextE(file1,' phi2', fontstr, 40, 25, 8, 1, string(EBSDdata.detphi2,format="(F6.2)"),'DETphi2','EBSDDetectorWidget_event')

file1 = WIDGET_BASE(block2, /ROW, XSIZE=440, /ALIGN_CENTER)
vals = ['TSL', 'HKL']
EBSDwidget_s.EulerConvention = CW_BGROUP(file1, $
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
file1 = WIDGET_BASE(block2, /ROW, XSIZE=440, /ALIGN_CENTER)
vals = ['Background Only','Full Pattern']
EBSDwidget_s.BGmode = CW_BGROUP(file1, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'Pattern Mode', $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='PATTERNMODE', $
                        SET_VALUE=EBSDdata.BGmode)

;------------------------------------------------------------
; min and max energy settings for filtered imaging
file1 = WIDGET_BASE(block2, /ROW, XSIZE=440, /ALIGN_CENTER)
tvals = strarr(EBSDdata.mcenergynumbin)
for i=0,EBSDdata.mcenergynumbin-1 do begin
  th = EBSDdata.mcenergymin + float(i)*EBSDdata.mcenergybinsize
  tvals[i] = string(th,format="(F5.2)")
end

file2 = WIDGET_LABEL(file1, VALUE='Energy Window:   Min ', font=fontstr)
EBSDwidget_s.EBSDminenergylist = WIDGET_DROPLIST(file1, $
			EVENT_PRO='EBSDDectectorWidget_event', $
			VALUE=tvals,$
			UVALUE='EBSDMINENERGYLIST', $
			/ALIGN_LEFT)
EBSDdata.Eminsel = 0
WIDGET_CONTROL, set_droplist_select=EBSDdata.Eminsel, EBSDwidget_s.EBSDminenergylist

file2 = WIDGET_LABEL(file1, VALUE='Max ', font=fontstr)
EBSDwidget_s.EBSDmaxenergylist = WIDGET_DROPLIST(file1, $
			EVENT_PRO='EBSDDectectorWidget_event', $
			VALUE=tvals,$
			UVALUE='EBSDMAXENERGYLIST', $
			/ALIGN_LEFT)
EBSDdata.Emaxsel = EBSDdata.mcenergynumbin-1
WIDGET_CONTROL, set_droplist_select=EBSDdata.Emaxsel, EBSDwidget_s.EBSDmaxenergylist



;------------------------------------------------------------
file1 = WIDGET_BASE(block2, /ROW, XSIZE=440, /ALIGN_CENTER)
vals = ['UL','LL','UR','LR']
EBSDwidget_s.PatternOrigin = CW_BGROUP(file1, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'Pattern Origin', $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSPATTERNORIGIN', $
                        SET_VALUE=EBSDdata.PatternOrigin)





;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,EBSDwidget_s.detectorbase,/REALIZE

; and hand over control to the xmanager
XMANAGER,"EBSDDetectorWidget",EBSDwidget_s.detectorbase,/NO_BLOCK


end

