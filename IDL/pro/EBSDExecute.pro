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
; CTEMsoft2013:EBSDExecute.pro
;--------------------------------------------------------------------------
;
; PROGRAM: EBSDExecute.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main routine for creation of nml file and execution of CTEMEBSD code
;
;> @date 05/22/14 MDG 1.0 first version
;--------------------------------------------------------------------------
pro EBSDExecute, status, single=single

; the keyword /single indicates that only one pattern should be computed

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common EBSDpatterns, pattern, image, finalpattern
common EBSD_anglearrays, euler, quaternions
common EBSDmasks, circularmask


; check whether the mask needs to be recomputed or not
s = size(circularmask)
dbin = 2^EBSDdata.detbinning
sm = min( [EBSDdata.detnumsx/dbin, EBSDdata.detnumsy/dbin] )
if (s[0] ne sm) then begin
  d = shift(dist(sm),sm/2,sm/2)
  d[where(d le sm/2)] = 1.0
  d[where(d gt sm/2)] = 0.0
  circularmask = fltarr(EBSDdata.detnumsx/dbin, EBSDdata.detnumsy/dbin)
  if (sm eq EBSDdata.detnumsx/dbin) then begin
    dm = (EBSDdata.detnumsy/dbin - sm)/2
    circularmask[0,dm] = d
  end else begin
    dm = (EBSDdata.detnumsx/dbin - sm)/2
    circularmask[dm,0] = d
  end
endif

v = double([EBSDdata.detax1,EBSDdata.detax2,EBSDdata.detax3])
axang = double([EBSDdata.detax1,EBSDdata.detax2,EBSDdata.detax3,EBSDdata.detax4])
euang = double([EBSDdata.detphi1,EBSDdata.detphi,EBSDdata.detphi2])

if keyword_set(single) then begin
; generate an angle input file
  phi1 = EBSDdata.detphi1
  Phi = EBSDdata.detphi
  phi2 = EBSDdata.detphi2
  openw,10,EBSDdata.pathname+'/'+'tmpangle.txt'
  printf,10,'eu'
  printf,10,'1'
  s = string(phi1,format="(F8.3)")+' '+ string(Phi,format="(F8.3)")+' '+ string(phi2,format="(F8.3)")
  printf,10,s
  close,10
    Core_Print, 'temporary single angle file created'
end

; if the keyword is not set, then we expect the angle file to already exist

; then generate the proper nml file
  openw,10,EBSDdata.pathname+'/'+'CTEMEBSDtmp.nml'
  printf,10,'&EBSDdata'
  printf,10,'L = '+string(EBSDdata.detL,FORMAT="(F10.2)")
  printf,10,'thetac = '+string(EBSDdata.dettheta,FORMAT="(F6.2)")
  printf,10,'delta = '+string(EBSDdata.detdelta,FORMAT="(F7.3)")
  printf,10,'numsx = '+string(EBSDdata.detnumsx,FORMAT="(I4)")
  printf,10,'numsy = '+string(EBSDdata.detnumsy,FORMAT="(I4)")
  printf,10,'xpc = '+string(EBSDdata.detxpc,FORMAT="(F8.3)")
  printf,10,'ypc = '+string(EBSDdata.detypc,FORMAT="(F8.3)")
  th = EBSDdata.mcenergymin + EBSDdata.Eminsel*EBSDdata.mcenergybinsize
  printf,10,'energymin = '+string(th,FORMAT="(F8.3)")
  th = EBSDdata.mcenergymin + EBSDdata.Emaxsel*EBSDdata.mcenergybinsize
  printf,10,'energymax = '+string(th,FORMAT="(F8.3)")
  if (EBSDData.EulerConvention eq 0) then begin
    printf,10,'eulerconvention = ''tsl'''
  end else begin
    printf,10,'eulerconvention = ''hkl'''
  end 
  printf,10,'FZfile = '''+EBSDdata.pathname+'/'+EBSDdata.mpfilename+''''
  printf,10,'energyfile = '''+EBSDdata.pathname+'/'+EBSDdata.mcfilename+''''
  printf,10,'datafile = '''+EBSDdata.EBSDpatternfilename+''''
  printf,10,'beamcurrent = '+string(EBSDdata.detbeamcurrent,FORMAT="(F9.2)")
  printf,10,'dwelltime = '+string(EBSDdata.detdwelltime,FORMAT="(F9.2)")
; and here are the imaging parameters that are only needed if we are NOT in single image mode
  if not keyword_set(single) then begin
; angle file
    if (EBSDdata.Pmode eq 1) then printf,10,'anglefile = '''+EBSDdata.EBSDanglefilename+''''
    if (EBSDdata.Pmode eq 2) then printf,10,'anglefile = '''+EBSDdata.EBSDdictfilename+''''
; binning parameter
    printf,10,'binning = '+string(2^EBSDdata.detbinning,format="(I1)")
; intensity scaling mode
    if (EBSDdata.PatternScaling eq 0) then begin
      printf,10,'scalingmode = ''lin'''
    end else begin
      printf,10,'scalingmode = ''gam'''
      printf,10,'gammavalue = '+string(EBSDdata.gammavalue,FORMAT="(F6.3)")
    end
; place holder for Pattern Origin setting
  end else begin
    printf,10,'anglefile = '''+EBSDdata.pathname+'/'+'tmpangle.txt'''
  end
; do we need to do an additional axis-angle pair rotation to all patterns ?
  if ( (axang[3] ne 0.D0) and (max(abs(v)) ne 0.D0) ) then begin
      s = string(axang[0],format="(F8.3)")+','+ string(axang[1],format="(F8.3)") +','+ string(axang[2],format="(F8.3)") +','+ string(axang[3],format="(F8.3)")
      printf,10,'axisangle = '+s
  end
  printf,10,'/'
  close,10

  if (EBSDdata.numangles gt 1000) then begin
    Core_Print,'',/blank
    Core_Print,'You are computing more than 1000 EBSPs; this will take a while...'
    Core_Print,'The program will not provide any further updates until the run has been completed.'
    Core_Print,'',/blank
  endif

; and finally run the CTEMEBSD program
  cmd = EBSDdata.f90exepath+'CTEMEBSD '+EBSDdata.pathname+'/'+'CTEMEBSDtmp.nml'
  spawn,cmd, cmdoutput, cmderroutput

  Core_Print,'CTEMEBSD program reported the following output '
  len = size(cmdoutput,/dimensions)
  for i=0,len[0]-1 do Core_Print,cmdoutput[i]

  if (cmderroutput[0] ne '') then begin
    Core_Print,'CTEMEBSD program reported an error !!!'
    len = size(cmderroutput,/dimensions)
    for i=0,len[0]-1 do Core_Print,cmderroutput[i]
  end else status = 1 

if keyword_set(single) then begin
  cmd = '/bin/rm '+EBSDdata.pathname+'/'+'tmpangle.txt'
  spawn, cmd
end

  cmd = '/bin/rm '+EBSDdata.pathname+'/'+'CTEMEBSDtmp.nml'
  spawn, cmd

; next, we need to load the pattern if we are in single mode
if keyword_set(single) then begin
  openu,1,EBSDdata.EBSDpatternfilename,/f77
  nsx = 0L
  nsy = 0L
  neu = 0L
  readu,1,nsx, nsy, neu
  Core_Print,' Dimensions read from pattern file : '+string(nsx,FORMAT="(I5)")+string(nsy,FORMAT="(I5)")+string(neu,FORMAT="(I8)")
  pattern = fltarr(nsx,nsy)
  readu,1,pattern
  close,1
end

; if we are not in single mode, then we need to load the angle file
if not keyword_set(single) then begin
  EBSDreadanglefile, EBSDdata.EBSDanglefilename, /list
end


end

