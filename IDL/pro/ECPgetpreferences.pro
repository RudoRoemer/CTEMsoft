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
; CTEMsoft2013:ECPgetpreferences.pro
;--------------------------------------------------------------------------
;
; PROGRAM: ECPgetpreferences.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief read the preferences file and initialize all relevant widgets
;
;> @date 11/23/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro ECPgetpreferences,dummy
 
;------------------------------------------------------------
; common blocks
common ECP_widget_common, widget_s
common ECP_data_common, data

; does the preferences file exist ?
rs = file_test(data.prefname)

if (rs eq 1) then begin
  s = ''
  i = 0
  openr,1,data.prefname
  readf,1,i
  data.nprefs = i

; next, do a little loop and read name:value pairs
  for i=0,data.nprefs-1 do begin
    readf,1,s
    spos = strpos(s,'::')
    nm = strmid(s,0,spos)
    val = strmid(s,spos+2)
    case nm of 
; root folder
	'ECProot': data.ECProot=val
; pattern output format
  	'ecpformat': data.ecpformat=fix(val)
; grid on or off ?
  	'ecpgrid': data.ecpgrid=fix(val)

; window locations
  	'xlocation': data.xlocation=float(val)
  	'ylocation': data.ylocation=float(val)
  	'ECPxlocation': data.ECPxlocation=float(val)
  	'ECPylocation': data.ECPylocation=float(val)
;
    else: MESSAGE,'unknown option for preferences file'
    endcase
  endfor

  close,1
end else begin
  s = ''
  cd,current=s
  data.ECProot=s
; prefs file does not exist yet, so let's create it with default values
  ECPwritepreferences
endelse

end

