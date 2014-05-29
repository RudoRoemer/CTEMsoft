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
; CTEMsoft2013:EBSDevent.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDevent.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief special event handler for all the CW_BGROUP calls, since CW_BGROUP does not support event_pro
;
;> @date 03/19/14 MDG 1.0 first version
;--------------------------------------------------------------------------
function EBSDevent, event

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common CommonCore, status, logmode, logunit

common EBSD_rawdata, accum_e, accum_z, MParray

if (EBSDdata.eventverbose eq 1) then help,event,/structure


WIDGET_CONTROL, event.id, GET_UVALUE = eventval         ;find the user value

IF N_ELEMENTS(eventval) EQ 0 THEN RETURN,eventval

CASE eventval OF

 'LOGFILE':  begin
; toggle the log mode 
		if (EBSDdata.logmode eq 0) then begin
		   Core_Print,'Turning log mode on',/blank
		 q = systime()
 		 z = strsplit(q,' ',/extract,/regex)
 		 EBSDdata.logname = 'EBSDDisplay'+z[0]+z[1]+z[2]+'_'+z[4]+'_'+z[3]+'.log'
		   Core_Print,'Log file: '+EBSDdata.logname
		 EBSDdata.logmode = 1
		 logmode = 1
		 openw,EBSDdata.logunit,EBSDdata.logname
		 EBSDdata.logfileopen = 1
		end else begin
		   Core_rint,'Turning log mode off',/blank
		 if (EBSDdata.logfileopen eq 1) then begin
		   close,EBSDdata.logunit
		   EBSDdata.logfileopen = 0
		 endif
	    	 EBSDdata.logmode = 0
		 logmode = 0
		endelse
	  endcase

  'EBSDFORMAT': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.EBSDformatbgroup
                EBSDdata.imageformat = fix(val[0])
          endcase


  'MCLS': begin
		WIDGET_CONTROL, get_value=val,EBSDwidget_s.MCLambertSelector
		EBSDdata.MCLSmode= fix(val[0])
	  endcase

  'MPLS': begin
		WIDGET_CONTROL, get_value=val,EBSDwidget_s.MPLambertSelector
		EBSDdata.MPLSmode= fix(val[0])
	  endcase

  'MCLsum': begin
		WIDGET_CONTROL, get_value=val,EBSDwidget_s.MCLambertMode
		EBSDdata.MCLSum =  fix(val[0])
		if (EBSDdata.MCLSum eq 1) then begin
		  EBSDshowMC, reform(total(accum_e,1))
		end else begin
		  EBSDshowMC, reform(accum_e[EBSDdata.Esel,*,*])
		end
	  endcase

  'EBSDEULERCONVENTION': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.EulerConvention
                EBSDdata.EulerConvention = fix(val[0])
	endcase

  'PATTERNMODE': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.BGmode
                EBSDdata.BGmode= fix(val[0])
	endcase

  'EBSPATTERNORIGIN': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.PatternOrigin
                EBSDdata.PatternOrigin = fix(val[0])
	  	EBSDshowPattern,/single
		vals = ['Upper Left','Lower Left','Upper Right','Lower Right']
		  Core_Print, 'Pattern origin set to '+vals[EBSDdata.PatternOrigin]
	endcase
 'EBSPATTERNSCALING': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.PatternScaling
                EBSDdata.PatternScaling = fix(val[0])
	  	EBSDshowPattern,/single
		vals = ['linear', 'gamma']
		  Core_Print, 'Pattern scaling set to '+vals[EBSDdata.PatternScaling]
	endcase

 'EBSDPATTERNBINNING': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.detbinning
                EBSDdata.detbinning= fix(val[0])
	  	EBSDshowPattern,/single
		vals = ['1','2','4','8']
		  Core_Print, 'Pattern binning set to '+vals[EBSDdata.detbinning]
	endcase

 'CIRCULARMASK': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.circularmask
                EBSDdata.showcircularmask= fix(val[0])
	  	if (EBSDdata.Pmode eq 0) then EBSDshowPattern,/single else EBSDshowPattern
		vals = ['Off','On']
		  Core_Print, 'Circular mask set to '+vals[EBSDdata.showcircularmask]
	endcase

 'PMODE': begin
                WIDGET_CONTROL, get_value=val,EBSDwidget_s.Pmode
                EBSDdata.Pmode = fix(val[0])
		vals = ['Single Pattern','Angle File','Dictionary']
; next we need to turn on those widgets that belong to the selected mode (sensitivity=1)
		if (EBSDdata.Pmode eq 0) then begin
		  WIDGET_CONTROL, EBSDwidget_s.DisplayEBSD, sensitive=1
		  WIDGET_CONTROL, EBSDwidget_s.EBSDgetanglefilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.PGdroplist, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.EBSDgetdictfilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
		    Core_Print, 'Pattern mode set to '+vals[EBSDdata.Pmode]
		end

		if (EBSDdata.Pmode eq 1) then begin
		  WIDGET_CONTROL, EBSDwidget_s.DisplayEBSD, sensitive=0
		  WIDGET_CONTROL, EBSDwidget_s.EBSDgetanglefilename, sensitive=1
;	  WIDGET_CONTROL, EBSDwidget_s.PGdroplist, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.EBSDgetdictfilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
		    Core_Print, 'Pattern mode set to '+vals[EBSDdata.Pmode]
		end

		if (EBSDdata.Pmode eq 2) then begin
		  WIDGET_CONTROL, EBSDwidget_s.DisplayEBSD, sensitive=0
		  WIDGET_CONTROL, EBSDwidget_s.EBSDgetanglefilename, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.PGdroplist, sensitive=0
;	  WIDGET_CONTROL, EBSDwidget_s.EBSDgetdictfilename, sensitive=0
;	  if ( (EBSDdata.Ncubochoric ne 0) and (EBSDdata.Dictpointgroup ne 0) and (EBSDdata.EBSDdictfilename ne '') ) then begin
;	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
;	  end else begin
;	    WIDGET_CONTROL, EBSDwidget_s.GoDict, sensitive=0
;	  end
;	    Core_Print, 'Pattern mode set to '+vals[EBSDdata.Pmode]
		    Core_Print, 'Not implemented in this program Release',/blank
		end
	endcase

else: MESSAGE, "Event User Value Not Found"

endcase

return,eventval
end 
