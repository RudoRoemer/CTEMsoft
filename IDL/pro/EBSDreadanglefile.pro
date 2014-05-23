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
; CTEMsoft2013:EBSDreadanglefile.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDreadanglefile.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief reads an angle input file and determines how many entries there are
;
;> @date 03/19/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDreadanglefile, fname

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

; the file must exist, since it was selected using the file widget
openr,10,fname
angletype = ''
readf,10,angletype
numangles = 0L
readf,10,numangles
; that's really all we need for now ... the f90 EBSD program is the one that 
; will actually read and use all the angles
close,10

EBSDdata.angletype = angletype
if (angletype eq 'eu') then begin
  WIDGET_CONTROL, set_value='Euler', EBSDwidget_s.angletype
end
if (angletype eq 'qu') then begin
  WIDGET_CONTROL, set_value='quaternion', EBSDwidget_s.angletype
end

WIDGET_CONTROL, set_value=string(numangles,FORMAT="(I8)"), EBSDwidget_S.numangles


end
