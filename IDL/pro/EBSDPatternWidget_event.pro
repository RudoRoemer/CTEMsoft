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
; CTEMsoft2013:EBSDPatternWidget_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDPatternWidget_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main widget event handler for display of the computed EBSD pattern
;
;> @date 05/22/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDPatternWidget_event, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall
common EBSDpatterns, pattern, image, finalpattern

if (EBSDdata.eventverbose eq 1) then help,event,/structure

; intercept the detector widget movement here 
if (event.id eq EBSDwidget_s.patternbase) then begin
  EBSDdata.patternxlocation = event.x
  EBSDdata.patternylocation = event.y-25
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

  CASE eventval OF
 'GAMMASLIDER': begin
	  WIDGET_CONTROL, get_value=val, EBSDwidget_s.gammaslider
	  EBSDdata.gammavalue = float(val[0]) 
	  EBSDshowPattern,/single
  	endcase

 'SAVEEBSDPATTERN': begin
; display a filesaving widget in the data folder with the file extension filled in
		delist = ['jpeg','tiff','bmp']
		de = delist[EBSDdata.imageformat]
		filename = DIALOG_PICKFILE(/write,default_extension=de,path=EBSDdata.pathname,title='enter filename without extension')
		im = finalpattern
		case de of
		    'jpeg': write_jpeg,filename,im,quality=100
		    'tiff': write_tiff,filename,reverse(im,2)
		    'bmp': write_bmp,filename,im
		 else: MESSAGE,'unknown file format option'
		endcase
	endcase

 'PATTERNCLOSE': begin
; kill the base widget
	  WIDGET_CONTROL, EBSDwidget_s.patternbase, /DESTROY
	    Core_Print,'EBSD Pattern Widget closed'
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 
