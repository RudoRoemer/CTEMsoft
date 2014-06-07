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
; CTEMsoft2013:EBSDgetpreferences.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDgetpreferences.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief read the preferences file and initialize all relevant widgets
;
;> @date 06/13/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro EBSDgetpreferences,noprint=noprint
 
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

; does the preferences file exist ?
rs = file_test(EBSDdata.prefname)

if (rs eq 1) then begin
  s = ''
  i = 0
  openr,1,EBSDdata.prefname
  readf,1,i
  EBSDdata.nprefs = i

; next, do a little loop and read name:value pairs
  for i=0,EBSDdata.nprefs-1 do begin
    readf,1,s
    spos = strpos(s,'::')
    nm = strmid(s,0,spos)
    val = strmid(s,spos+2)
    case nm of 
; root folder
	'EBSDroot': EBSDdata.EBSDroot=val
	'f90exepath': EBSDdata.f90exepath=val

; various parameters
  	'detl': EBSDdata.detL = float(val)
  	'dettheta': EBSDdata.dettheta = float(val)
  	'detdelta': EBSDdata.detdelta = float(val)
  	'detnumsx': EBSDdata.detnumsx = long(val)
  	'detnumsy': EBSDdata.detnumsy = long(val)
  	'detxpc': EBSDdata.detxpc = float(val)
  	'detypc': EBSDdata.detypc = float(val)
  	'detbinning': EBSDdata.detbinning = long(val)
  	'detbeamcurrent': EBSDdata.detbeamcurrent = float(val)
  	'detdwelltime': EBSDdata.detdwelltime = float(val)

; window locations
  	'xlocation': EBSDdata.xlocation = float(val)
  	'ylocation': EBSDdata.ylocation = float(val)
  	'EBSDxlocation': EBSDdata.EBSDxlocation = float(val)
  	'EBSDylocation': EBSDdata.EBSDylocation = float(val)
  	'Patternxlocation': EBSDdata.patternxlocation = float(val)
  	'Patternylocation': EBSDdata.patternylocation = float(val)
  	'Detectorxlocation': EBSDdata.Detectorxlocation = float(val)
  	'Detectorylocation': EBSDdata.Detectorylocation = float(val)
  	'MCxlocation': EBSDdata.MCxlocation = float(val)
  	'MCylocation': EBSDdata.MCylocation = float(val)
  	'MPxlocation': EBSDdata.MPxlocation = float(val)
  	'MPylocation': EBSDdata.MPylocation = float(val)
  	'Detectorxlocation': EBSDdata.Detectorxlocation = float(val)
  	'Detectorylocation': EBSDdata.Detectorylocation = float(val)

    else: MESSAGE,'unknown option for preferences file'
    endcase
  endfor

  close,1
end else begin
  s = ''
  cd,current=s
  EBSDdata.EBSDroot=s
; prefs file does not exist yet, so let's create it with default values
  if not keyword_set(noprint) then Core_Print,'Creating preferences file '+EBSDdata.prefname
  if keyword_set(noprint) then EBSDwritepreferences,/noprint else EBSDwritepreferences
endelse

end

