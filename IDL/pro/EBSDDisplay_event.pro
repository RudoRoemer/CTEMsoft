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
; CTEMsoft2013:EBSDDisplay_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDDisplay_event.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Electron backscatter diffraction pattern display event handler
;
;
;> @date 03/19/14 MDG 1.0 initial version
;--------------------------------------------------------------------------
pro EBSDDisplay_event, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata


if (EBSDdata.eventverbose eq 1) then help,event,/structure

if (event.id eq EBSDwidget_s.base) then begin
  EBSDdata.xlocation = event.x
  EBSDdata.ylocation = event.y-25
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value
  
  CASE eventval OF
  	'MCDISPLAY': begin
; create the Monte Carlo display widget
		EBSDMCDisplayWidget
	endcase

  	'MPDISPLAY': begin
; create the Master Pattern display widget
		EBSDMCDisplayWidget
	endcase

  	'MCFILE': begin
; ask the user to select the data file
		EBSDgetfilename,validfile,/MCFILE
; read the data file and populate all the relevant fields
		if (validfile eq 1) then EBSDreaddatafile,/MCFILE
; activate the MC Display button

; and close any other open widgets
	endcase

  	'MPFILE': begin
; ask the user to select the data file
		EBSDgetfilename,validfile,/MPFILE
; read the data file and populate all the relevant fields
		if (validfile eq 1) then EBSDreaddatafile,/MPFILE
; activate both the MC and MP Display buttons

; and close any other open widgets
	endcase

	'DETECTOR': begin
		EBSDDetectorWidget
	endcase

 	'QUIT': begin
		EBSDwritepreferences
; do a general cleanup of potentially open widgets
 		Core_Print,'Quitting program',/blank
		WIDGET_CONTROL, EBSDwidget_s.base, /DESTROY
		!EXCEPT=1
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 



