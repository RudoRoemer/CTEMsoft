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

if keyword_set(single) then begin
; first determine whether or not these euler angles need to be further rotated
; this would be the case if detax4 is not zero
  v = double([EBSDdata.detax1,EBSDdata.detax2,EBSDdata.detax3])
  if ( (EBSDdata.detax4 ne 0.D0) and (max(abs(v)) ne 0.D0) ) then begin
; first create the rotation quaternion
    om2 = EBSDdata.detax4*!dtor*0.5D0
    v = v/sqrt(total(v^2))
    qr = [cos(om2),sin(om2)*v] ; unit quaternion
; then convert the Euler angles to a quaternion
    ee = 0.5D0*[EBSDdata.detphi1,EBSDdata.detphi,EBSDdata.detphi2]*!dtor
    cPhi = cos(ee[1])
    sPhi = sin(ee[1])
    cm = cos(ee[0]-ee[2])
    sm = sin(ee[0]-ee[2])
    cp = cos(ee[0]+ee[2])
    sp = sin(ee[0]+ee[2])
    q = !pi*cp
    qe = [-q, sPhi*cm, sPhi*sm, cPhi*sp]
; next, multiply the two quaternions
    q = [0.D0,Core_quatmult(qr,qe)]  ; to account for 0 start of arrays in IDL
; convert this quaternion back to Euler angles
    q03 = q[1]^2+q[4]^2
    q12 = q[2]^2+q[3]^2
    chi = sqrt(q03*q12)

    if (chi eq 0.D0) then begin
      if (q12 eq 0.D0) then  begin
        Phi = 0.D0
        phi2 = 0.D0   		; arbitrarily due to degeneracy
        phi1 = atan(-2.D0*q[1]*q[4],q[1]^2-q[4]^2)/!dtor
      end else begin
        Phi = !dpi
        phi2 = 0.D0   		; arbitrarily due to degeneracy
        phi1 = atan(2.D0*q[2]*q[3],q[2]^2-q[3]^2)/!dtor
      end 
    end else begin 		; this is not a special degenerate case
      Phi = atan( 2.D0*chi, q03-q12 )/!dtor
      chi = 1.D0/chi
      phi1 = atan( (-q[1]*q[3]+q[2]*q[4])*chi, (-q[1]*q[2]-q[3]*q[4])*chi )/!dtor
      phi2 = atan( (q[1]*q[3]+q[2]*q[4])*chi, (-q[1]*q[2]+q[3]*q[4])*chi )/!dtor
    end 
      s = string(phi1,format="(F8.3)")+' '+ string(Phi,format="(F8.3)")+' '+ string(phi2,format="(F8.3)")
      Core_Print,'Euler angles rotated to '+s
  end else begin
    phi1 = EBSDdata.detphi1
    Phi = EBSDdata.detphi
    phi2 = EBSDdata.detphi2
  end

; generate an angle input file
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
  printf,10,'anglefile = '''+EBSDdata.pathname+'/'+'tmpangle.txt'''
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
  printf,10,'/'
  close,10

; and finally run the CTEMEBSD program
  cmd = EBSDdata.f90exepath+'CTEMEBSD '+EBSDdata.pathname+'/'+'CTEMEBSDtmp.nml'
  spawn,cmd, cmdoutput, cmderroutput

  Core_Print,'CTEMEBSD program reported the following output '
  len = size(cmdoutput,/dimensions)
  for i=0,len[0]-1 do Core_Print,cmdoutput[i]

  if (cmderroutput ne '') then begin
    Core_Print,'CTEMEBSD program reported an error !!!'
    len = size(cmderroutput,/dimensions)
    for i=0,len[0]-1 do Core_Print,cmderroutput[i]
  end else status = 1 

if keyword_set(single) then begin
  cmd = '/bin/rm '+EBSDdata.pathname+'/'+'tmpangle.txt'
  spawn, cmd
end


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


end

