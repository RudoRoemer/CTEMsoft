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
; CTEMsoft2013:EBSDwritepreferences.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDwritepreferences.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief write the preferences file
;
;> @date 06/13/13 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro EBSDwritepreferences,noprint=noprint
 
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

; prefs file
  openw,1,EBSDdata.prefname
  nprefs = 25
  EBSDdata.nprefs = nprefs
  printf,1,nprefs
  printf,1,'EBSDroot::'+EBSDdata.EBSDroot
  printf,1,'EBSDMCroot::'+EBSDdata.EBSDMCroot
  printf,1,'f90exepath::'+EBSDdata.f90exepath

  printf,1,'detl::'+string(EBSDdata.detL,format="(F9.2)")
  printf,1,'dettheta::'+string(EBSDdata.dettheta,format="(F6.2)")
  printf,1,'detdelta::'+string(EBSDdata.detdelta,format="(F6.2)")
  printf,1,'detnumsx::'+string(EBSDdata.detnumsx,format="(I6)")
  printf,1,'detnumsy::'+string(EBSDdata.detnumsy,format="(I6)")
  printf,1,'detxpc::'+string(EBSDdata.detxpc,format="(F7.2)")
  printf,1,'detypc::'+string(EBSDdata.detypc,format="(F7.2)")
  printf,1,'detbinning::'+string(EBSDdata.detbinning,format="(I3)")
  printf,1,'detbeamcurrent::'+string(EBSDdata.detbeamcurrent,format="(D9.2)")
  printf,1,'detdwelltime::'+string(EBSDdata.detdwelltime,format="(D9.2)")

; window locations
  printf,1,'xlocation::'+string(EBSDdata.xlocation,format="(F6.1)")
  printf,1,'ylocation::'+string(EBSDdata.ylocation,format="(F6.1)")
  printf,1,'EBSDxlocation::'+string(EBSDdata.EBSDxlocation,format="(F6.1)")
  printf,1,'EBSDylocation::'+string(EBSDdata.EBSDylocation,format="(F6.1)")
  printf,1,'Detectorxlocation::'+string(EBSDdata.Detectorxlocation,format="(F6.1)")
  printf,1,'Detectorylocation::'+string(EBSDdata.Detectorylocation,format="(F6.1)")
  printf,1,'Patternxlocation::'+string(EBSDdata.patternxlocation,format="(F6.1)")
  printf,1,'Patternylocation::'+string(EBSDdata.patternylocation,format="(F6.1)")
  printf,1,'MCxlocation::'+string(EBSDdata.MCxlocation,format="(F6.1)")
  printf,1,'MCylocation::'+string(EBSDdata.MCylocation,format="(F6.1)")
  printf,1,'MPxlocation::'+string(EBSDdata.MPxlocation,format="(F6.1)")
  printf,1,'MPylocation::'+string(EBSDdata.MPylocation,format="(F6.1)")
; and close the file
  close,1

  if not keyword_set(noprint) then Core_Print,'The preferences file '+EBSDdata.prefname+' was successfully saved '

end

