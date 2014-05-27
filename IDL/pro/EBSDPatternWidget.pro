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
; CTEMsoft2013:EBSDPatternWidget.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDPatternWidget.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main widget for display of the computed EBSD pattern
;
;> @date 05/22/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDPatternWidget, single=single

; the keyword /single indicates that only one pattern is available

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("EBSDPatternWidget") NE 0) then begin
  print,'EBSDPatternWidget is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; create the top level widget
EBSDwidget_s.patternbase = WIDGET_BASE(TITLE='Pattern Display Widget', $
                        /ROW, $
                        XSIZE=max( [600,EBSDdata.detnumsx+100] ), $
                        /ALIGN_CENTER, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='EBSDPatternWidget_event', $
                        XOFFSET=EBSDdata.patternxlocation, $
                        YOFFSET=EBSDdata.patternylocation)


block1 = WIDGET_BASE(EBSDwidget_s.patternbase, $
			XSIZE=600, $
			/ALIGN_TOP, $
			/COLUMN)

if keyword_set(single) then begin
 EBSDdata.currentdisplaywidgetmode = 0
;------------------------------------------------------------
;------------------------------------------------------------
  file1 = WIDGET_BASE(block1, /COLUMN, /FRAME, YPAD=8, XSIZE=600, /ALIGN_LEFT)
  file2 = WIDGET_LABEL(file1, VALUE='Pattern parameters', font=fontstrlarge, /ALIGN_LEFT, /FRAME)


;------------------------------------------------------------
  file2 = WIDGET_BASE(file1, /ROW, XSIZE=600, /ALIGN_CENTER)
  vals = ['1','2','4','8']
  EBSDwidget_s.detbinning = CW_BGROUP(file2, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'Detector Binning', $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSDPATTERNBINNING', $
                        SET_VALUE=EBSDdata.detbinning)

;------------------------------------------------------------
  file2 = WIDGET_BASE(file1, /ROW, XSIZE=600, /ALIGN_CENTER)
  vals = ['UL','LL','UR','LR']
  EBSDwidget_s.PatternOrigin = CW_BGROUP(file2, $
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
  file2 = WIDGET_BASE(file1, /ROW, XSIZE=600, /ALIGN_CENTER)
  vals = ['linear','gamma']
  EBSDwidget_s.PatternScaling = CW_BGROUP(file2, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'Pattern Scaling', $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSPATTERNSCALING', $
                        SET_VALUE=EBSDdata.PatternScaling)


  file2 = WIDGET_BASE(file1, /ROW, XSIZE=600, /ALIGN_CENTER)
; here's a slider to select the gamma setting ...
  EBSDwidget_s.gammaslider = CW_FSLIDER(file2, $
;		EVENT_PRO='EBSDPatternWidget_event', $
			/EDIT, $
			MINIMUM = 0.01, $
			MAXIMUM = 2.00, $
			FORMAT = "(F4.2)", $
			TITLE = 'Gamma Correction Factor', $
			XSIZE = 400, $
			VALUE = EBSDdata.gammavalue, $
			UVALUE = 'GAMMASLIDER')


; Result = CW_FSLIDER( Parent [, /DOUBLE] [, /DRAG] [, /EDIT] [, FORMAT=string] [, /FRAME] [, MAXIMUM=value] [, MINIMUM=value] [, SCROLL=units] [, /SUPPRESS_VALUE] [, TAB_MODE=value] [, TITLE=string] [, UNAME=string] [, UVALUE=value] [, VALUE=initial_value] [, XSIZE=length | {, /VERTICAL [, YSIZE=height]}] )

;------------------------------------------------------------
;------------------------------------------------------------
; and here's the display window itself
EBSDwidget_s.Patterndraw = WIDGET_DRAW(block1, $
                        COLOR_MODEL=2, $
                        RETAIN=2, $
                        /FRAME, $
                        XSIZE=EBSDdata.detnumsx, $
                        YSIZE=EBSDdata.detnumsy)

; and the min-max indicators
block4 = WIDGET_BASE(block1, /ROW, /ALIGN_CENTER)
EBSDwidget_s.Patternmin = Core_WText(block4, 'min/max ',fontstr, 75, 25, 15, 1, string(EBSDdata.Patternmin,FORMAT="(F9.3)"))
EBSDwidget_s.Patternmax = Core_WText(block4, '/',fontstr, 5, 25, 15, 1, string(EBSDdata.Patternmax,FORMAT="(F9.3)"))

; and a save button
block4 = WIDGET_BASE(block1, /ROW, /ALIGN_CENTER)
saveEBSDPattern = WIDGET_BUTTON(block4, $
                        VALUE='Save', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDPatternWidget_event', $
                        /FRAME, $
                        UVALUE='SAVEEBSDPATTERN', $
                        /ALIGN_LEFT)

; and the save format selector
vals = ['jpeg','tiff','bmp']
EBSDwidget_s.EBSDformatbgroup = CW_BGROUP(block4, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'File Format', $
                        /FRAME, $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSDFORMAT', $
                        SET_VALUE=EBSDdata.imageformat)

;------------------------------------------------------------
;------------------------------------------------------------
; and here is the Close button
file1 = WIDGET_BASE(block1, XSIZE=340, /ALIGN_LEFT, /ROW)

EBSDwidget_s.PatternClose = WIDGET_BUTTON(file1, $
                                UVALUE='PATTERNCLOSE', $
                                VALUE='Close', $
                                EVENT_PRO='EBSDPatternWidget_event', $
                                SENSITIVE=1, $
                                /FRAME)


end else begin
; this is not single image mode, so we do not need to display the 
; image scaling parameters and such; we'll have a simpler interface
; that allows the user to browse through the series of images, 
; and save individual ones or all of them (after warning)

 EBSDdata.currentdisplaywidgetmode = 1
;------------------------------------------------------------
;------------------------------------------------------------
  file1 = WIDGET_BASE(block1, /COLUMN, YPAD=8, XSIZE=600, /ALIGN_LEFT)
  file2 = WIDGET_BASE(file1, /ROW, XSIZE=600, /ALIGN_LEFT)

saveEBSDPattern = WIDGET_BUTTON(file2, $
                        VALUE='Previous', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDPatternWidget_event', $
                        /FRAME, $
                        UVALUE='PREVIOUSEBSDPATTERN', $
                        /ALIGN_LEFT)

saveEBSDPattern = WIDGET_BUTTON(file2, $
                        VALUE='Next', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDPatternWidget_event', $
                        /FRAME, $
                        UVALUE='NEXTEBSDPATTERN', $
                        /ALIGN_LEFT)

; display the euler angles/quaternion for the currently displayed pattern
; we'll do this as a simple non-editable text widget

file2 = WIDGET_BASE(block1, /ROW, XSIZE=600, /ALIGN_CENTER)
EBSDwidget_s.angledisplay = Core_WText(file2,'Orientation:', fontstr, 120, 25, 60, 1, ' ')


; then the display window
EBSDwidget_s.Patterndraw = WIDGET_DRAW(block1, $
                        COLOR_MODEL=2, $
                        RETAIN=2, $
                        /FRAME, $
                        XSIZE=EBSDdata.detnumsx, $
                        YSIZE=EBSDdata.detnumsy)

; and the min-max indicators
block4 = WIDGET_BASE(block1, /ROW, /ALIGN_CENTER)
EBSDwidget_s.Patternmin = Core_WText(block4, 'min/max ',fontstr, 75, 25, 15, 1, string(EBSDdata.Patternmin,FORMAT="(F9.3)"))
EBSDwidget_s.Patternmax = Core_WText(block4, '/',fontstr, 5, 25, 15, 1, string(EBSDdata.Patternmax,FORMAT="(F9.3)"))

; a save all button
block4 = WIDGET_BASE(block1, /ROW, /ALIGN_CENTER)
saveEBSDPattern = WIDGET_BUTTON(block4, $
                        VALUE='SaveAll', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDPatternWidget_event', $
                        /FRAME, $
                        UVALUE='SAVEALLEBSDPATTERNS', $
                        /ALIGN_LEFT)

; and a save button
saveEBSDPattern = WIDGET_BUTTON(block4, $
                        VALUE='Save', $
                        /NO_RELEASE, $
                        EVENT_PRO='EBSDPatternWidget_event', $
                        /FRAME, $
                        UVALUE='SAVEEBSDPATTERN', $
                        /ALIGN_LEFT)

; and the save format selector
vals = ['jpeg','tiff','bmp']
EBSDwidget_s.EBSDformatbgroup = CW_BGROUP(block4, $
                        vals, $
                        /ROW, $
                        /NO_RELEASE, $
                        /EXCLUSIVE, $
                        FONT=fontstr, $
                        LABEL_LEFT = 'File Format', $
                        /FRAME, $
                        EVENT_FUNC ='EBSDevent', $
                        UVALUE='EBSDFORMAT', $
                        SET_VALUE=EBSDdata.imageformat)

;------------------------------------------------------------
;------------------------------------------------------------
; and here is the Close button
file1 = WIDGET_BASE(block1, XSIZE=340, /ALIGN_LEFT, /ROW)

EBSDwidget_s.PatternClose = WIDGET_BUTTON(file1, $
                                UVALUE='PATTERNCLOSE', $
                                VALUE='Close', $
                                EVENT_PRO='EBSDPatternWidget_event', $
                                SENSITIVE=1, $
                                /FRAME)
end




;------------------------------------------------------------
; realize the widget structure
WIDGET_CONTROL,EBSDwidget_s.patternbase,/REALIZE

; realize the draw widget
WIDGET_CONTROL, EBSDwidget_s.Patterndraw, GET_VALUE=drawID
EBSDwidget_s.PatternDrawID = drawID

; and display the pattern with the current intensity settings
if (EBSDdata.Pmode eq 0) then EBSDshowPattern,/single else EBSDshowPattern

; and hand over control to the xmanager
XMANAGER,"EBSDPatternWidget",EBSDwidget_s.patternbase,/NO_BLOCK

end

