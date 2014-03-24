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
; CTEMsoft2013:EBSDMCDisplayWidget.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDMCDisplayWidget.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Display widget for Monte Carlo simulation results
;
;> @date 03/19/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDMCDisplayWidget,dummy
;
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

; the next common block contains all the raw data needed to generate the EBSD patterns
common EBSD_rawdata, accum_e, accum_z, MParray

;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("EBSDMCDisplayWidget") NE 0) then begin
  print,'EBSDMCDisplayWidget is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; create the top level widget
if (EBSDdata.MCMPboth eq 0) then begin
  EBSDwidget_s.MCdisplaybase = WIDGET_BASE(TITLE='Monte Carlo Display Widget', $
                        /COLUMN, $
                        XSIZE=20+max([ 500, 2*EBSDdata.mcimx+1 ]), $
                        /ALIGN_CENTER, $
			MBAR = menubar, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDMCDisplayWidget_event', $
                        XOFFSET=EBSDdata.MCxlocation, $
                        YOFFSET=EBSDdata.MCylocation)
end else begin
  EBSDwidget_s.MCdisplaybase = WIDGET_BASE(TITLE='Master Pattern & Monte Carlo Display Widget', $
                        /COLUMN, $
                        XSIZE=20+max([ 500, 2*EBSDdata.mcimx+1 ]) + max( [ 500, 2*EBSDdata.mpimx+1 ] ), $
                        /ALIGN_CENTER, $
			MBAR = menubar, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDMCDisplayWidget_event', $
                        XOFFSET=EBSDdata.MCxlocation, $
                        YOFFSET=EBSDdata.MCylocation)
end

menu1 = WIDGET_BUTTON(menubar, VALUE='Projection Mode', /MENU)
b1 = WIDGET_BUTTON(menu1, VALUE='Lambert [square]', UVALUE='LAMBERTS')
b2 = WIDGET_BUTTON(menu1, VALUE='Lambert [circle]', UVALUE='LAMBERTC')

menu2 = WIDGET_BUTTON(menubar, VALUE='Display Mode', /MENU)
if (EBSDdata.MCMPboth eq 1) then begin
  b2 = WIDGET_BUTTON(menu2, VALUE='Individual Energy/Master Bin', UVALUE='DISPEBIN')
  b2 = WIDGET_BUTTON(menu2, VALUE='Simple Energy Sum/Weighted Master', UVALUE='DISPESUM')
  b2 = WIDGET_BUTTON(menu2, VALUE='Simple Energy Sum RGB/blank', UVALUE='DISPESUMRGB')
end else begin
  b2 = WIDGET_BUTTON(menu2, VALUE='Individual Energy Bin', UVALUE='DISPEBIN')
  b2 = WIDGET_BUTTON(menu2, VALUE='Simple Energy Sum', UVALUE='DISPESUM')
  b2 = WIDGET_BUTTON(menu2, VALUE='Simple Energy Sum RGB', UVALUE='DISPESUMRGB')
endelse

;------------------------------------------------------------
if (EBSDdata.MCMPboth eq 0) then begin
	block1 = WIDGET_BASE(EBSDwidget_s.MCdisplaybase, $
			/FRAME, $
                        XSIZE=max([ 500, 2*EBSDdata.mcimx+1 ]) , $
			/ALIGN_CENTER, $
			/COLUMN)
end else begin
	block1 = WIDGET_BASE(EBSDwidget_s.MCdisplaybase, $
			/FRAME, $
                        XSIZE=max([ 500, 2*EBSDdata.mcimx+1 ]) + max( [ 500, 2*EBSDdata.mpimx+1 ] ), $
			/ALIGN_CENTER, $
			/COLUMN)
endelse

block2 = WIDGET_BASE(block1, /ROW, /ALIGN_CENTER)

; here's a slider to select the energy window ...
EBSDwidget_s.MCslider = WIDGET_SLIDER(block2, $
			EVENT_PRO='EBSDMCDisplayWidget_event', $
			MINIMUM = 1, $
			MAXIMUM = EBSDdata.mcenergynumbin, $
			SENSITIVE = 1, $
			TITLE = 'Select an energy', $
			XSIZE = 400, $
			VALUE = 1, $
			UVALUE = 'MCSLIDER', $
			/ALIGN_CENTER)

; and right next to it we display the actual energy in a text box
energy = EBSDdata.mcenergymin + EBSDdata.Esel * EBSDdata.mcenergybinsize
EBSDwidget_s.MCenergyval =  WIDGET_TEXT(block2, $
			VALUE=string(energy,format="(F5.2)"), $
			XSIZE=10, $
			YSIZE=1, $
			/ALIGN_RIGHT)
;------------------------------------------------------------
block2 = WIDGET_BASE(block1, /ROW, /ALIGN_CENTER)
block3 = WIDGET_BASE(block2, /COLUMN, /ALIGN_CENTER)

; and here's the MC display window itself
EBSDwidget_s.MCdraw = WIDGET_DRAW(block3, $
                        COLOR_MODEL=2, $
                        RETAIN=2, $
                        /FRAME, $
                        XSIZE=2*EBSDdata.mcimx+1, $
                        YSIZE=2*EBSDdata.mcimy+1)

; and the min-max indicators
block4 = WIDGET_BASE(block3, /ROW, /ALIGN_CENTER)
EBSDwidget_s.MCmin = Core_WText(block4, 'min/max ',fontstr, 75, 25, 15, 1, string(EBSDdata.MCmin,FORMAT="(F9.1)"))
EBSDwidget_s.MCmax = Core_WText(block4, '/',fontstr, 5, 25, 15, 1, string(EBSDdata.MCmax,FORMAT="(F9.1)"))

; and a save button
saveEBSDMC = WIDGET_BUTTON(block4, $
                        VALUE='Save', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDMCDisplayWidget_event', $
                        /FRAME, $
                        UVALUE='SAVEEBSDMC', $
                        /ALIGN_LEFT)



;------------------------------------------------------------
; and the MP window, if needed
if (EBSDdata.MCMPboth eq 1) then begin
        block3 = WIDGET_BASE(block2, /COLUMN, /ALIGN_CENTER)
        EBSDwidget_s.MPdraw = WIDGET_DRAW(block3, $
                        COLOR_MODEL=2, $
                        RETAIN=2, $
                        /FRAME, $
                        XSIZE=2*EBSDdata.mpimx+1, $
                        YSIZE=2*EBSDdata.mpimy+1)

; and the min-max indicators
        block4 = WIDGET_BASE(block3, /ROW, /ALIGN_CENTER)
        EBSDwidget_s.MPmin = Core_WText(block4, 'min/max ',fontstr, 75, 25, 15, 1, string(EBSDdata.MPmin,FORMAT="(F9.1)"))
        EBSDwidget_s.MPmax = Core_WText(block4, '/',fontstr, 5, 25, 15, 1, string(EBSDdata.MPmax,FORMAT="(F9.1)"))

; and a save button
        saveEBSDMP = WIDGET_BUTTON(block4, $
                        VALUE='Save', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDMCDisplayWidget_event', $
                        /FRAME, $
                        UVALUE='SAVEEBSDMP', $
                        /ALIGN_LEFT)
endif


;------------------------------------------------------------
; and finally a close button
file1 = WIDGET_BASE(block1, $
			XSIZE=500, $
			/FRAME, $
			/ROW)

file2 = WIDGET_BUTTON(file1, $
                        VALUE='Close', $
                        UVALUE='CLOSEMC', $
                        EVENT_PRO='EBSDMCDisplayWidget_event', $
                        SENSITIVE=1, $
                        /FRAME)

vals = ['jpeg','tiff','bmp']
EBSDwidget_s.EBSDformatbgroup = CW_BGROUP(file1, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstrlarge, $
                        LABEL_LEFT = 'File Format', $
                        /FRAME, $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSDFORMAT', $
                        SET_VALUE=EBSDdata.imageformat)

;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,EBSDwidget_s.MCdisplaybase,/REALIZE

; realize the draw widget
WIDGET_CONTROL, EBSDwidget_s.MCdraw, GET_VALUE=drawID
EBSDwidget_s.MCdrawID = drawID
if (EBSDdata.MCMPboth eq 1) then begin
	WIDGET_CONTROL, EBSDwidget_s.MPdraw, GET_VALUE=drawID
	EBSDwidget_s.MPdrawID = drawID
endif

; and next we draw the first entry of the accum_e array
wset,EBSDwidget_s.MCdrawID
tvscl,reform(accum_e[0,*,*])

; add a min/max indicator as well as a save button


; and hand over control to the xmanager
XMANAGER,"EBSDMCDisplayWidget",EBSDwidget_s.MCdisplaybase,/NO_BLOCK

end 

