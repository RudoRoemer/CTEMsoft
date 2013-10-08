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
; CTEMsoft2013:CBEDDisplay_event.pro
;--------------------------------------------------------------------------
;
; PROGRAM: CBEDDisplay_event.pro
;
;> @author Marc De Graef, Carnegie Melon University
;
;> @brief main event handler
;
;> @date 06/13/13 MDG 1.0 first version
;--------------------------------------------------------------------------
pro CBEDDisplay_event, event

;------------------------------------------------------------
; common blocks
common CBED_widget_common, widget_s
common CBED_data_common, data
;common STEM_detectordata, STEMdata, STEMcimage, BFdisk, DFdisk, clickablemap, STEMsectormaps, STEMsectors, STEMsectorranges
;common STEM_CBEDpatterns, CBED, CBEDdisplay
; and this one is used to create the blue channel of the detector plot
;common STEM_circles, th, cth, sth, blue, diskpos
; the next common block contains all the raw data needed to generate the CBED patterns
;common STEM_rawdata, indices, offsets, kperp, rawdata
;common STEM_masks, ktpg, ktpgang, BFmask, HAADFmask, BFindices, HAADFindices, BFcnt, HAADFcnt


if (data.eventverbose eq 1) then help,event,/structure

if (event.id eq widget_s.base) then begin
  data.xlocation = event.x
  data.ylocation = event.y-25
end else begin

  WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value
  
; IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

  CASE eventval OF

  'LOADLACBEDFILE': begin
; loading a new file means that a bunch of variables need to be reset
; ask the user to select an input geometry file
		CBEDgetfilename,/LACBED

; read the geometry file and populate all the relevant fields
		CBEDreaddatafile,/LACBED

; reset the diffraction mode
;	data.diffractionmode = 0
;	WIDGET_CONTROL,SET_VALUE=data.difffractionmode,widget_s.difffractionmode 
;	WIDGET_CONTROL, widget_s.gosector,sensitive=1
;	WIDGET_CONTROL, widget_s.clearsector,sensitive=1
;	WIDGET_CONTROL, widget_s.aprad,sensitive=0
;	WIDGET_CONTROL, get_value=val,widget_s.detsegm
;	data.detsegm = fix(val[0])

; and draw the detector pattern for the current parameters
;	STEMdetectorsetup
	  endcase

  'LOADMBCBEDFILE': begin
; loading a new file means that a bunch of variables need to be reset
; ask the user to select an input geometry file
		CBEDgetfilename,/MBCBED

; read the geometry file and populate all the relevant fields
		CBEDreaddatafile,/MBCBED

; reset the diffraction mode
;	data.diffractionmode = 0
;	WIDGET_CONTROL,SET_VALUE=data.difffractionmode,widget_s.difffractionmode 
;	WIDGET_CONTROL, widget_s.gosector,sensitive=1
;	WIDGET_CONTROL, widget_s.clearsector,sensitive=1
;	WIDGET_CONTROL, widget_s.aprad,sensitive=0
;	WIDGET_CONTROL, get_value=val,widget_s.detsegm
;	data.detsegm = fix(val[0])

; and draw the detector pattern for the current parameters
;	STEMdetectorsetup
		
	  endcase

 'QUIT': begin
; do a general cleanup
		CBEDprint,'Quitting program',/blank
;	if (XRegistered("STEMCBEDWidget") NE 0) then WIDGET_CONTROL, widget_s.cbedbase, /DESTROY
;	if (XRegistered("STEMImageWidget") NE 0) then WIDGET_CONTROL, widget_s.imagebase, /DESTROY
; write the preferences file
;	STEMwritepreferences
; kill a bunch of large variables

; close the log file if it is open
		if (data.logfileopen eq 1) then begin
		  close,data.logunit
		endif
; wait briefly
		wait,2.0
; and finally kill the base widget
		WIDGET_CONTROL, widget_s.base, /DESTROY
		!EXCEPT=1
	endcase

  else: MESSAGE, "Event User Value Not Found"

  endcase

endelse

end 
