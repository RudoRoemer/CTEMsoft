;
; Copyright (c) 2013-2014, Marc De Graef/Carnegie Mellon University
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
; CTEMsoft2013:EBSDMCDisplayWidget_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDMCDisplayWidget_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main event handler for Monte Carlo Display mode
;
;> @date 10/15/13 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDMCDisplayWidget_event, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common EBSD_rawdata, accum_e, accum_z, MParray, MParraysum
common projections, mcxcircle, mcycircle, mpxcircle, mpycircle, mcSPxcircle, mcSPycircle, mpSPxcircle, mpSPycircle 

common Image_common, MCimage, MPimage


if (EBSDdata.eventverbose eq 1) then help,event,/structure

; intercept the image widget movement here 
if (event.id eq EBSDwidget_s.MCdisplaybase) then begin
  EBSDdata.MCxlocation = event.x
  EBSDdata.MCylocation = event.y-25
    Core_Print,' Window moved to location ('+string(fix(EBSDdata.MCxlocation),format="(I4)")+','+string(fix(EBSDdata.MCylocation),format="(I4)")+')'
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

; IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

  CASE eventval OF

  'MCSLIDER':  begin
		WIDGET_CONTROL, get_value=val,EBSDwidget_s.mcslider
		EBSDdata.Esel = fix(val[0]) - 1
		EBSDshowMC
	endcase

  'DISPEBIN': begin
		EBSDdata.MCLSum = 0
		EBSDshowMC
		WIDGET_CONTROL, EBSDwidget_s.MCslider, sensitive=1
	endcase

  'DISPESUM': begin
		EBSDdata.MCLSum = 1
		EBSDshowMC
		WIDGET_CONTROL, EBSDwidget_s.MCslider, sensitive=0
	endcase

  'DISPESUMRGB': begin
		EBSDdata.MCLSum = 2
		EBSDshowMC
		WIDGET_CONTROL, EBSDwidget_s.MCslider, sensitive=0
	endcase

  'DISPWSUM': begin
		EBSDdata.MCLSum = 3
		EBSDshowMC
		WIDGET_CONTROL, EBSDwidget_s.MCslider, sensitive=0
	endcase

  'LAMBERTS': begin
		EBSDdata.MCLSmode = 0
		EBSDshowMC
	endcase

  'LAMBERTC': begin
		EBSDdata.MCLSmode = 1
		EBSDshowMC
	endcase

  'STEREOP': begin
		EBSDdata.MCLSmode = 2
		EBSDshowMC
	endcase

  'ASYMPOS':  begin
		val = event.index
		if (val eq 0) then EBSDdata.Asymsel = -1 else EBSDdata.Asymsel = fix(val)-1
		EBSDshowMC
	endcase

  'SAVEEBSDMC': begin
; display a filesaving widget in the data folder with the file extension filled in
		delist = ['jpeg','tiff','bmp']
		de = delist[EBSDdata.imageformat]
		filename = DIALOG_PICKFILE(/write,default_extension=de,path=EBSDdata.mcpathname,title='enter filename without extension')
		sz = size(MCimage,/dimensions)
		if (sz[0] eq 3) then begin
                 case de of
                    'jpeg': write_jpeg,filename,MCimage,true=1,quality=100
                    'tiff': begin
                        rgb = MCimage 
                        for i=0,2 do begin
                          slice = reform(rgb[i,*,*])
                          rgb[i,0:*,0:*] = reverse(slice,2)
                        end
                        write_tiff,filename,rgb
                      endcase
                    'bmp': write_bmp,filename,MCimage,/RGB
                   else: MESSAGE,'unknown file format option'
                  endcase
		end else begin
		  im = bytscl(MCimage)
		  case de of
		    'jpeg': write_jpeg,filename,im,quality=100
		    'tiff': write_tiff,filename,reverse(im,2)
		    'bmp': write_bmp,filename,im
		   else: MESSAGE,'unknown file format option'
		  endcase
		endelse
	endcase

  'SAVEEBSDMP': begin
; display a filesaving widget in the data folder with the file extension filled in
		delist = ['jpeg','tiff','bmp']
		de = delist[EBSDdata.imageformat]
		filename = DIALOG_PICKFILE(/write,default_extension=de,path=EBSDdata.pathname,title='enter filename without extension')
		im = bytscl(MPimage)
		case de of
		    'jpeg': write_jpeg,filename,im,quality=100
		    'tiff': write_tiff,filename,reverse(im,2)
		    'bmp': write_bmp,filename,im
		 else: MESSAGE,'unknown file format option'
		endcase
	endcase

  'CLOSEMC': begin
; kill the base widget
		WIDGET_CONTROL, EBSDwidget_s.MCdisplaybase, /DESTROY
		  Core_Print,'EBSD Display Widget closed'
	endcase

  else: MESSAGE, "Event User Value Not Found"

  END

endelse


end
