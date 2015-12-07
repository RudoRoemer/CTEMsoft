;
; Copyright (c) 2015, Marc De Graef/Carnegie Mellon University
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
; EMsoft:Core_getenv.pro
;--------------------------------------------------------------------------
;
; PROGRAM: Core_getenv.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief read the environment variable for the EMsoft path from a json file
;
;> @date 10/08/15 MDG 1.0 initial implementation 
;> @date 11/11/15 MDG 1.1 added functionality for release version library location
;> @date 11/21/15 MDG 1.2 corrected typo in librarylocation variable
;> @date 12/07/15 MDG 1.3 added LibraryLocation variable to EMsoft configuration file to simplify things
;--------------------------------------------------------------------------
function Core_getenv,data=data

common getenv_common, librarylocation

; first parse the json configuration file
result = json_parse('~/.config/EMsoft/EMsoftConfig.json',/toarray)


; then extract either the main path or the data path, depending on the keyword
if keyword_set(data) then begin
  z = result['EMdatapathname']
end else begin
  z = result['EMsoftpathname']
; In the development version of this package, the dynamical library will be located in the 
; the EMsoftpathname/Build/Bin folder.
; In the release version, the user is supposed to place the libEMSoftLib.dylib file
; in the /opt/local/lib/libgcc folder.  The EMsoftLibraryLocation variable in the 
; EMsoftConfig.json file defines where that file is actually located.
  librarylocation = result['EMsoftLibraryLocation']
endelse

return,z
end
