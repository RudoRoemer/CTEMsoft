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
; create the top level widget
EBSDwidget_s.MCdisplaybase = WIDGET_BASE(TITLE='Monte Carlo Display Widget', $
                        /COLUMN, $
                        XSIZE=20+max([ 500, 2*EBSDdata.mcimx+1 ]), $
                        /ALIGN_CENTER, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDMCDisplayWidget_event', $
                        XOFFSET=EBSDdata.MCxlocation, $
                        YOFFSET=EBSDdata.MCylocation)

block1 = WIDGET_BASE(EBSDwidget_s.MCdisplaybase, $
			/FRAME, $
                        XSIZE=max([ 500, 2*EBSDdata.mcimx+1 ]), $
			/ALIGN_CENTER, $
			/COLUMN)

values = ['Lambert', 'Modified Lambert']
EBSDwidget_s.MCLambertSelector = CW_BGROUP(block1, $
			values, $
			/FRAME, $
                        LABEL_LEFT='Projection Mode', $
			/COLUMN, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=EBSDdata.MCLSmode, $
                        EVENT_FUNC='EBSDevent', $
			UVALUE='MCLS')

values = ['Individual', 'Sum']
EBSDwidget_s.MCLambertMode = CW_BGROUP(block1, $
			values, $
			/FRAME, $
                        LABEL_LEFT='Pattern Mode ', $
			/COLUMN, $
			/NO_RELEASE, $
			/EXCLUSIVE, $
			SET_VALUE=EBSDdata.MCLsum, $
                        EVENT_FUNC='EBSDevent', $
			UVALUE='MCLsum')

; here's a slider to select the energy window ...
file1 = WIDGET_SLIDER(block1, $
			EVENT_PRO='EBSDMCDisplayWidget_event', $
			MINIMUM = 1, $
			MAXIMUM = EBSDdata.mcenergynumbin, $
			TITLE = 'Select an energy', $
			XSIZE = 400, $
			VALUE = 1, $
			UVALUE = 'MCSLIDER', $
			/ALIGN_CENTER)

; and here's the display window itself
EBSDwidget_s.MCdraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			XSIZE=2*EBSDdata.mcimx+1, $
			YSIZE=2*EBSDdata.mcimy+1)

;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,EBSDwidget_s.MCdisplaybase,/REALIZE

; realize the draw widget
WIDGET_CONTROL, EBSDwidget_s.MCdraw, GET_VALUE=drawID
EBSDdata.MCdrawID = drawID


; and hand over control to the xmanager
XMANAGER,"EBSDMCDisplayWidget",EBSDwidget_s.MCdisplaybase,/NO_BLOCK

end 

