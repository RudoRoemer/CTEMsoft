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
; CTEMsoft2013:EBSDgetfilename.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDgetfilename.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief Display an interface and ask user to select a MC or MP output file
;
;> @date 03/19/14 MDG 1.0 first attempt 
;--------------------------------------------------------------------------
pro EBSDgetfilename,validfile,MCFILE=MCFILE,MPFILE=MPFILE
 
;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata

  validfile = 0

  s = ''
  cd,current = s
  EBSDdata.homefolder = s
  if (EBSDdata.EBSDroot eq 'undefined') then begin
    EBSDdata.EBSDroot = EBSDdata.homefolder
  end 

  rootpath = EBSDdata.EBSDroot

  if keyword_set(MCFILE) then begin
    res=dialog_pickfile(title='Select a valid Monte Carlo data file',path=rootpath,filter='*mc*.*;*MC*.*')
    if (res eq '') then begin
	  Core_Print,'No selection made'
	  goto, skip
    end
; check to make sure that this file is of the correct type
; read the first fnlen=132 characters and compare the non-null characters to
; the program name 'CTEMMC.f90'
    openu,1,res,/f77
    progname = bytarr(132)   
    readu,1,progname
    progname = strtrim(string(progname))
    close,1
    if (progname eq 'CTEMMC.f90') then begin
	validfile = 1
  	finfo = file_info(res)
	EBSDdata.mcfilesize = finfo.size
; find the last folder separator
	spos = strpos(res,'/',/reverse_search)
	dpos = strpos(res,'.',/reverse_search)
	plen = strlen(res)
	EBSDdata.pathname = strmid(res,0,spos)
	EBSDdata.mcfilename = strmid(res,spos+1)
	EBSDdata.suffix = strmid(res,dpos+1)
	EBSDdata.EBSDroot = EBSDdata.pathname

  	WIDGET_CONTROL, SET_VALUE=EBSDdata.mcfilename, EBSDwidget_s.mcfilename

  	Core_Print,' full path '+res
  	Core_Print,' path '+EBSDdata.pathname
  	Core_Print,' data file '+EBSDdata.mcfilename
  	Core_Print,' suffix '+EBSDdata.suffix
    end else begin
  	Core_Print,' This file is not of the correct Monte Carlo data type ',/blank
        goto, skip
    endelse
  end else begin
    res=dialog_pickfile(title='Select a valid Master Pattern data file',path=rootpath,filter='*Master*.*;*master*.*;*MASTER*.*')
    if (res eq '') then begin
	  Core_Print,'No selection made'
	  goto, skip
    end
; check to make sure that this file is of the correct type
; read the first fnlen=132 characters and compare the non-null characters to
; the program name 'CTEMEBSDmaster.f90'
    openu,1,res,/f77
    progname = bytarr(132)   
    readu,1,progname
    progname = strtrim(string(progname))
    close,1
    if (progname eq 'CTEMEBSDmaster.f90') then begin
	validfile = 1
  	finfo = file_info(res)
	EBSDdata.mpfilesize = finfo.size
; find the last folder separator
	spos = strpos(res,'/',/reverse_search)
	dpos = strpos(res,'.',/reverse_search)
	plen = strlen(res)
	EBSDdata.pathname = strmid(res,0,spos)
	EBSDdata.mpfilename = strmid(res,spos+1)
	EBSDdata.suffix = strmid(res,dpos+1)
	EBSDdata.EBSDroot = EBSDdata.pathname

  	WIDGET_CONTROL, SET_VALUE=EBSDdata.mpfilename, EBSDwidget_s.mpfilename

  	Core_Print,' full path '+res
  	Core_Print,' path '+EBSDdata.pathname
  	Core_Print,' data file '+EBSDdata.mpfilename
  	Core_Print,' suffix '+EBSDdata.suffix
    end else begin
  	Core_Print,' This file is not of the correct Master Pattern data type ',/blank
        goto, skip
    endelse
  endelse

skip:
end

