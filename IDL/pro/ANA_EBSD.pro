@ANAEBSD_getpreferences
@ANA_getvendor

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
; EMsoft:ANA_EBSD.pro
;--------------------------------------------------------------------------

; this program is a simple GUI to do a couple of basic analyses
; of EBSD vendor hdf5 files:  
;
; - display of individual patterns with or without background
; - computation of dot product maps
; - our version of PRIAS
; - display of all the datasets from the HDF file
;
;> @date 06/18/15 MDG 1.0 first version
;
;---------------------------------------------------------------------------

pro ANA_EBSD,dummy

; common blocks



!EXCEPT=0

;------------------------------------------------------------
; make sure that this program isn't already running
if (XRegistered("ANA_EBSD") NE 0) then begin
  print,'ANA_EBSD is already running ... (if it is not, please restart your IDL session)'
  return
end

;------------------------------------------------------------
; define a few structures (one for widgets, and one for data)
ANAwidget = {anawidgetstruct, $
             base:long(0), $      
             logodraw:long(0), $      
            }

ANAdata   = {anadatastruct, $
	     scversion: '', $	
	     xlocation: float(0.0), $		; main widget x-location (can be modified and stored in preferences file)
	     ylocation: float(0.0), $		; main widget y-location (can be modified and stored in preferences file)
             scrdimx:0L, $                      ; display area x size in pixels 
             scrdimy:0L $                       ; display area y size in pixels 
            }

; a few font strings
fontstr='-adobe-new century schoolbook-bold-r-normal--14-100-100-100-p-87-iso8859-1'
fontstrlarge='-adobe-new century schoolbook-medium-r-normal--20-140-100-100-p-103-iso8859-1'
fontstrsmall='-adobe-new century schoolbook-medium-r-normal--14-100-100-100-p-82-iso8859-1'

;------------------------------------------------------------
; set the display window size to 80% of the current screen size (but be careful with double screens ... )
device,decomposed = 0
device, GET_SCREEN_SIZE = scr
anadata.scrdimy = scr[1] * 0.8
anadata.scrdimx = 0.75 * anadata.scrdimy   ; doing it this way avoids problems with multiple screens
anadata.xlocation = anadata.scrdimx / 8.0
anadata.ylocation = anadata.scrdimx / 8.0

;------------------------------------------------------------
; does the preferences file exist ?  If not, create it, otherwise read it
;ANA_EBSDgetpreferences

;------------------------------------------------------------
; create the top level widget
anawidget.base = WIDGET_BASE(TITLE='EBSD Vendor HDF5 File Analysis Program', $
                        /ROW, $
                        XSIZE=1300, $
                        /ALIGN_LEFT, $
			/TLB_MOVE_EVENTS, $
			EVENT_PRO='ANAEBSD_event', $
                        XOFFSET=anadata.xlocation, $
                        YOFFSET=anadata.ylocation)

block0 = WIDGET_BASE(anawidget.base, /COLUMN)

block1 = WIDGET_BASE(block0, /FRAME, /COLUMN)

file1 = WIDGET_BASE(block1, /ROW, XSIZE=650, /ALIGN_LEFT)

anawidget.logodraw = WIDGET_DRAW(block1, $
			COLOR_MODEL=2, $
			RETAIN=2, $
			/FRAME, $
			/ALIGN_CENTER, $
			XSIZE=600, $
			YSIZE=200)








; after selecting the input file, the program must analyze the data file to
; figure out what's in it...; this may be different for different vendor files,
; so let's first figure out what the Manufacturer is...





end
