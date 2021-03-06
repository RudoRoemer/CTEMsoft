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
; CTEMsoft2013:KosselPatternWidget.pro
;--------------------------------------------------------------------------
;
; PROGRAM: KosselPatternWidget.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Create the CBED control widget
;
;> @date 10/15/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------

pro KosselPatternWidget,dummy

;------------------------------------------------------------
; common blocks
common Kossel_widget_common, widget_s
common Kossel_data_common, data
common BED_rawdata, rawdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

;------------------------------------------------------------
; create the top level widget
widget_s.KosselPatternbase= WIDGET_BASE(TITLE='Kossel Widget', $
                        /COLUMN, $
                        XSIZE=20 + max([513,data.datadims[0]]), $
                        /ALIGN_CENTER, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='KosselPatternWidget_event', $
                        XOFFSET=data.Kosselxlocation, $
                        YOFFSET=data.Kosselylocation)

;------------------------------------------------------------
; create the various vertical blocks
; block 1 deals with the thickness selection
block1 = WIDGET_BASE(widget_s.KosselPatternbase, /FRAME, /COLUMN)

;----------
file1 = WIDGET_BASE(block1, $
			/ROW, $
                        XSIZE=400, $
			/ALIGN_LEFT)

label2 = WIDGET_LABEL(file1, $
			VALUE='Integration Depth [nm]', $
			FONT=fontstrlarge, $
			XSIZE=200, $
			YSIZE=25, $
			/ALIGN_LEFT)

tvals = strarr(data.datadims[2])
for i=0,data.datadims[2]-1 do begin
  th = data.startthick + float(i)*data.thickinc
  tvals[i] = string(th,format="(F5.1)")
end

widget_s.Kosselthicklist = WIDGET_DROPLIST(file1, $
			EVENT_PRO='KosselPatternWidget_event', $
			VALUE=tvals,$
			UVALUE='KosselTHICKLIST', $
			/ALIGN_LEFT)

data.thicksel = 0
WIDGET_CONTROL, set_droplist_select=data.thicksel, widget_s.Kosselthicklist

;----------
file1 = WIDGET_BASE(block1, $
			/ROW, $
                        XSIZE=500, $
			/ALIGN_LEFT)

label4 = WIDGET_LABEL(file1, $
			VALUE='Blur Radius', $
			FONT=fontstrlarge, $
			XSIZE=110, $
			YSIZE=25, $
			/ALIGN_LEFT)

widget_s.blur=  WIDGET_TEXT(file1, $
		VALUE=string(data.blur,FORMAT="(F6.2)"),$
		XSIZE=10, $
		YSIZE=1, $
		/EDITABLE, $
                EVENT_PRO='KosselPatternWidget_event', $
		UVALUE='BLUR', $
		/ALIGN_LEFT)



;---------- save pattern
file1 = WIDGET_BASE(block1, $
			/ROW, $
			/ALIGN_LEFT)

vals = ['jpeg','tiff','bmp']
widget_s.Kosselformatbgroup = CW_BGROUP(file1, $
			vals, $
			/ROW, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			FONT=fontstrlarge, $
			LABEL_LEFT = 'File Format', $
			/FRAME, $
                        EVENT_FUNC ='Kosselevent', $
			UVALUE='KosselFORMAT', $
			SET_VALUE=data.Kosselformat)


;------------ a few control buttons
block2 = WIDGET_BASE(widget_s.KosselPatternbase, $
			/FRAME, $
			/ROW)

saveMBCBED = WIDGET_BUTTON(block2, $
			VALUE='Save', $
			/NO_RELEASE, $
                        EVENT_PRO='KosselPatternWidget_event', $
			/FRAME, $
			UVALUE='SAVEKossel', $
			/ALIGN_LEFT)


closeKossel = WIDGET_BUTTON(block2, $
			VALUE='Close', $
			/NO_RELEASE, $
                        EVENT_PRO='KosselPatternWidget_event', $
			/FRAME, $
			UVALUE='CLOSEKossel', $
			/ALIGN_RIGHT)

;------------ the actual pattern 
widget_s.Kosseldraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/BUTTON_EVENTS, $
			UVALUE = 'GETCOORDINATES', $
			/ALIGN_CENTER, $
			XSIZE=data.datadims[0], $
			YSIZE=data.datadims[1])


;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,widget_s.KosselPatternbase,/REALIZE

; realize the draw widgets
WIDGET_CONTROL, widget_s.Kosseldraw, GET_VALUE=drawID
widget_s.KosseldrawID = drawID

; and hand over control to the xmanager
XMANAGER,"KosselPatternWidget",widget_s.KosselPatternbase,/NO_BLOCK

Kosselshow

end
