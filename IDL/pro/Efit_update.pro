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
; EMsoft:Efit_update.pro
;--------------------------------------------------------------------------
;
; PROGRAM: Efit_update.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief called by Efit_amoeba; returns dot product or mutual information
;
;> @date 10/20/15 MDG 1.0 first version
;--------------------------------------------------------------------------
function Efit_update, param

;------------------------------------------------------------
; common blocks
common Efit_widget_common, Efitwidget_s
common Efit_data_common, Efitdata
common fontstrings, fontstr, fontstrlarge, fontstrsmall

common CommonCore, status, logmode, logunit
common FitParameters, nFit, fitName, defValue, fitValue, fitStep, fitOnOff, fitManualStep, fitManualUpDown, fitUserLabel, fitStepLabel, fitOnOffLabel, fitUpLabel, fitDownLabel, fitManualStepLabel, fitIterations


common EBSD_EMsoft, MCxtalname, MCmode, nsx, nsy, EkeV, Ehistmin, Ebinsize, depthmax, depthstep, MCsig, MComega, $
                    numEbins, numzbins, accum_e, accum_z, Masterenergyfile, npx, npy, nnE, numset, mLPNH, mLPSH, Masterxtalname, expEBSDpattern, EBSDpattern

common cancelcommon, cancel
common Efitdisplaycommon, mask, maskready, expvector

nset = size(param,/dimensions)
nset = nset[0]

; set up the ipar and fpar arrays; all integers must be long64 !!!!
; none of the integer parameters are refinable, so this part is always the same
nipar = long(8)
ipar = lon64arr(nipar)
ipar[0] = long64(1) ; long(recompute) ; 1 if rgx, rgy, rgz detector arrays need to be computed, 0 if not (arrays will have save status)
if ((fitIterations gt 0L) and (total(fitOnOff[0:4]) eq 0)) then begin
  ipar[0] = long64(0) 
endif
ipar[1] = long64(Efitdata.detnumsx)
ipar[2] = long64(Efitdata.detnumsy)
ipar[3] = long64(numEbins)
ipar[4] = long64(nsx)
ipar[5] = long64(nsy)
ipar[6] = long64((nsx-1)/2)
ipar[7] = long64(npx)

nfpar = long(13)
fpar = fltarr(nfpar)

; loop through all potential fitting parameters
; this is a little tricky, since there are many potential combinations of fitting parameters
np = 0
if (fitOnOff[0] eq 1) then begin
  fpar[6] = param[np]
  np += 1
end else fpar[6] = fitValue[0]

if (fitOnOff[1] eq 1) then begin
  fpar[4] = param[np]
  np += 1
end else fpar[4] = fitValue[1]

if (fitOnOff[2] eq 1) then begin
  fpar[0] = param[np]
  np += 1
end else fpar[0] = fitValue[2]

if (fitOnOff[3] eq 1) then begin
  fpar[1] = param[np]
  np += 1
end else fpar[1] = fitValue[3]

if (fitOnOff[4] eq 1) then begin
  Efitdata.detgamma = param[np]
  np += 1
end else Efitdata.detgamma = fitValue[4]

; here are the remaining constants
fpar[2] = Efitdata.detdelta
fpar[3] = Efitdata.detMCsig
fpar[5] = Efitdata.dettheta
fpar[7] = Efitdata.detbeamcurrent
fpar[8] = Efitdata.detdwelltime

; and finally convert the Euler angle triplet (in degrees) to a quaternion
; they are either refinable or not
if (fitOnOff[5] eq 1) then begin
    phi1 = param[np]
    np += 1
end else phi1 = fitValue[5]
if (fitOnOff[6] eq 1) then begin
    phi = param[np]
    np += 1
end else phi = fitValue[6]
if (fitOnOff[7] eq 1) then begin
    phi2 = param[np]
    np += 1
end else phi2 = fitValue[7]
Efitdata.quaternion = Core_eu2qu( [phi1, phi, phi2] )
fpar[9:12] = Efitdata.quaternion[0:3]

; initialize the simulated pattern array
EBSDpattern = fltarr(Efitdata.detnumsx,Efitdata.detnumsy)

callname = 'SingleEBSDPatternWrapper'

res = call_external(Efitdata.EMsoftpathname+'Build/Bin/libEMSoftLib.dylib', callname, $
      nipar, nfpar, long(numEbins), ipar[4], long(npx), ipar[1], ipar[2], ipar, fpar, float(accum_e), mLPNH, mLPSH, EBSDpattern, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)

if (res ne 1.0) then begin
  Core_print,'SingleEBSDPatternWrapper return code = '+string(res,format="(F4.1)")
end 

; apply the correct intensity scaling
Epat = EBSDpattern^Efitdata.detgamma

; try a high pass filter
if (Efitdata.hipassonoff eq 1) then begin
    hipass = DIGITAL_FILTER(Efitdata.hipasscutoff,1.0,50,7)
    Epat = Convol(Epat,hipass)
endif

; remove any ramp
if (Efitdata.ramponoff eq 1) then begin
    Eav = mean(Epat)
    ramp = sfit(Epat,1)
    Epat = Epat/ramp
    Epat = Epat * Eav/mean(Epat)
endif

if (Efitdata.smoothval ne 0) then begin
    Epat = smooth(Epat,2*Efitdata.smoothval+1)
endif

npixels = Efitdata.detnumsx*Efitdata.detnumsy

; any gradient or laplacian operators to be applied to the patterns ?
if (Efitdata.preproc ne 0) then begin
  if (Efitdata.preproc eq 1) then begin
    v = reform(roberts(EBSDpattern)*mask,npixels)
    vl = sqrt(total(v*v))
    simvector = v/vl
; store the laplacian of the experimental pattern in a normalized column vector
    v = reform(roberts(float(expEBSDpattern))*mask,npixels)
    vl = sqrt(total(v*v))
    expvector = v/vl
  endif
end else begin
  v = reform(float(EBSDpattern),npixels)
  vl = sqrt(total(v*v))
  simvector = v/vl
  v = reform(float(expEBSDpattern),npixels)
  vl = sqrt(total(v*v))
  expvector = v/vl
endelse



if (Efitdata.ConvCrit eq 0) then begin   ; dot product approach
; and get the dot product
  cval = total(simvector * expvector)
end else begin    ; mutual information approach
  V = bytarr(2,npixels)
  V[0,0:*] = bytscl(expvector[0:*]) 
  V[1,0:*] = bytscl(simvector[0:*]) 
  bs = 1
  h = Core_histnd(V,bs,/normalize)
  cval = Core_mind(h)
endelse

; then show the pattern ...
wset,Efitdata.drawID
Efit_showpattern

; update the iterations counter
fitIterations += 1L

; and update the parameter value widgets
np = 0
for i=0,nFit-1 do begin
  if (fitOnOff[i] eq 1) then begin
    fitValue[i] = param[np]
    np += 1
  endif
endfor
for i=0,nFit-1 do WIDGET_CONTROL, set_value=string(fitValue[i],format="(F9.2)"), Efitwidget_s.fitValue[i]

; and finally, put the cval value in the appropriate progress box in the main widget
WIDGET_CONTROL, set_value=string(cval,FORMAT="(F12.6)"), Efitwidget_s.progress

; and we return the negative value of the convergence parameters, since Efit_amoeba finds a minimum, not a maximum
return, -cval
end
