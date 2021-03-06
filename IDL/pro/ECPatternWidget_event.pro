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
; CTEMsoft2013:ECPatternWidget_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: ECPatternWidget_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main event handler
;
;> @date 06/13/13 MDG 1.0 first version
;--------------------------------------------------------------------------
pro ECPatternWidget_event, event

;------------------------------------------------------------
; common blocks
common ECP_widget_common, widget_s
common ECP_data_common, data
common ECP_rawdata, rawdata

if (data.eventverbose eq 1) then help,event,/structure

; intercept the image widget movement here 
if (event.id eq widget_s.ECPatternbase) then begin
  data.ECPxlocation = event.x
  data.ECPylocation = event.y-25
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

; IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

  CASE eventval OF
 'GETCOORDINATES': begin
	  if (event.press eq 1B) then begin    ; only act on clicks, not on releases
	    data.cx = (event.x - data.xmid) / data.dgrid
	    data.cy = (event.y - data.xmid) / data.dgrid
	    WIDGET_CONTROL, SET_VALUE=string(data.cx,format="(F6.3)"), widget_s.cx
	    WIDGET_CONTROL, SET_VALUE=string(data.cy,format="(F6.3)"), widget_s.cy
	  end
	endcase

 'ECPTHICKLIST': begin
	  data.thicksel = event.index

; and display the selected ECPattern
          ECPshow
	endcase

 'BLUR':  begin
	    WIDGET_CONTROL, get_value=val,widget_s.blur
	    data.blur= float(val[0])
	    WIDGET_CONTROL, SET_VALUE=string(data.blur,FORMAT="(F6.3)"), widget_s.blur
	    ECPshow
	endcase

 'PATROT':  begin
	    WIDGET_CONTROL, get_value=val,widget_s.patrot
	    data.patrot = float(val[0])
	    WIDGET_CONTROL, SET_VALUE=string(data.patrot,FORMAT="(F6.3)"), widget_s.patrot
	    ECPshow
	endcase
 

 'SAVEECP': begin
; display a filesaving widget in the data folder with the file extension filled in
		delist = ['jpeg','tiff','bmp']
		de = delist[data.ecpformat]
		filename = DIALOG_PICKFILE(/write,default_extension=de,path=data.pathname,title='enter filename without extension')
	        im = tvrd()
		case de of
		  'jpeg': write_jpeg,filename,im,quality=100
		  'tiff': write_tiff,filename,reverse(im,2)
		  'bmp': write_bmp,filename,im
		 else: MESSAGE,'unknown file format option'
		endcase
	  endcase

 'CLOSEECP': begin
; kill the base widget
		WIDGET_CONTROL, widget_s.ECPatternbase, /DESTROY
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 
