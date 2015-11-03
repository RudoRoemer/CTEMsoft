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
; EMsoft:ECPatternWidget_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: ECPatternWidget_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main widget event handler for display of the computed ECP pattern
;
;> @date 10/30/15 MDG 1.0 first version
;--------------------------------------------------------------------------
pro ECPatternWidget_event, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall
common EBSDpatterns, pattern, image, finalpattern
common EBSD_anglearrays, euler, quaternions
common EBSDmasks, circularmask

common ECPdata, ECPattern

if (EBSDdata.eventverbose eq 1) then help,event,/structure

EMdatapathname = Core_getenv(/data)

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
	  ECPshowPattern,/single
  	endcase

 'SAVEECPATTERN': begin
; display a filesaving widget in the data folder with the file extension filled in
		delist = ['jpeg','tiff','bmp']
		de = delist[EBSDdata.imageformat]
		filename = DIALOG_PICKFILE(/write,default_extension=de,path=EBSDdata.pathname,title='enter filename without extension')
		ECPshowpattern, /nodisplay
		if (EBSDdata.showcircularmask eq 1) then im = finalpattern*byte(circularmask) else im=finalpattern
		case de of
		    'jpeg': write_jpeg,filename,im,quality=100
		    'tiff': write_tiff,filename,reverse(im,2)
		    'bmp': write_bmp,filename,im
		 else: MESSAGE,'unknown file format option'
		endcase
	endcase

 'SAVEALLECPATTERNS': begin
; display a filesaving widget in the data folder with the file extension filled in
		delist = ['jpeg','tiff','bmp']
		de = delist[EBSDdata.imageformat]
		fn = DIALOG_PICKFILE(/write,default_extension=de,path=EBSDdata.pathname,title='enter prefix for image series file name')
		fn = strsplit(fn,'.',/extract)
		for i=0,EBSDdata.numangles-1 do begin
  		  pattern = reform(ECPattern[*,*,i])
		  ECPshowpattern, /single, /nodisplay
		  if (EBSDdata.showcircularmask eq 1) then im = finalpattern*byte(circularmask) else im=finalpattern
		  filename = fn[0]+string(i+1,format="(I5.5)")+'.'+fn[1]
		  case de of
		    'jpeg': write_jpeg,filename,im,quality=75
		    'tiff': write_tiff,filename,reverse(im,2)
		    'bmp': write_bmp,filename,im
		   else: MESSAGE,'unknown file format option'
		  endcase
		end
  		close,1
		  Core_Print,'All image files generated'
	endcase


 'NEXTECPATTERN': begin
	  EBSDdata.currentpatternID += 1
	  if (EBSDdata.currentpatternID ge EBSDdata.numangles) then EBSDdata.currentpatternID = 0
	  i = EBSDdata.currentpatternID
	  if (EBSDdata.angletype eq 'eu') then begin
	    st = string(i+1,format="(I5.5)")+': '+string(euler[0,i],format="(F7.2)")+', '+string(euler[1,i],format="(F7.2)")+', '+string(euler[2,i],format="(F7.2)")
	  end else begin
	    st = string(i+1,format="(I5.5)")+': '+string(quaternions[0,i],format="(F7.2)")+', '+string(quaternionseuler[1,i],format="(F7.2)")+ $
		', '+string(quaternions[2,i],format="(F7.2)")+', '+string(quaternions[3,i],format="(F7.2)")
	  end
	  WIDGET_CONTROL, set_value = st, EBSDwidget_s.angledisplay
	  ECPshowpattern
	endcase

 'PREVIOUSECPATTERN': begin
	  EBSDdata.currentpatternID -= 1
	  if (EBSDdata.currentpatternID lt 0) then EBSDdata.currentpatternID = EBSDdata.numangles-1
	  i = EBSDdata.currentpatternID
	  if (EBSDdata.angletype eq 'eu') then begin
	    st = string(i+1,format="(I5.5)")+': '+string(euler[0,i],format="(F7.2)")+', '+string(euler[1,i],format="(F7.2)")+', '+string(euler[2,i],format="(F7.2)")
	  end else begin
	    st = string(i+1,format="(I5.5)")+': '+string(quaternions[0,i],format="(F7.2)")+', '+string(quaternionseuler[1,i],format="(F7.2)")+ $
		', '+string(quaternions[2,i],format="(F7.2)")+', '+string(quaternions[3,i],format="(F7.2)")
	  end
	  WIDGET_CONTROL, set_value = st, EBSDwidget_s.angledisplay
	  ECPshowpattern
	endcase


 'PATTERNCLOSE': begin
; kill the base widget
	  WIDGET_CONTROL, EBSDwidget_s.patternbase, /DESTROY
	    Core_Print,'ECP Pattern Widget closed'
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 
