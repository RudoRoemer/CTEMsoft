;
; Copyright (c) 2013-2015, Marc De Graef/Carnegie Mellon University
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
; EMsoft:ECPExecute.pro
;--------------------------------------------------------------------------
;
; PROGRAM: ECPExecute.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief main routine for creation of nml file and execution of CTEMEBSD code
;
;> @date 05/22/14 MDG 1.0 first version
;> @date 04/14/15 MDG 1.1 added HDF5 support
;> @date 10/30/15 MDG 2.0 copied from EBSDExecute.pro and adapted for ECP
;--------------------------------------------------------------------------
pro ECPExecute, status, ECpattern, single=single

; the keyword /single indicates that only one pattern should be computed

;------------------------------------------------------------
; common blocks
common EBSD_widget_common, EBSDwidget_s
common EBSD_data_common, EBSDdata
common EBSDpatterns, pattern, image, finalpattern
common EBSD_anglearrays, euler, quaternions
common EBSDmasks, circularmask

common EBSD_rawdata, accum_e, accum_z, mLPNH, mLPSH


; check whether the mask needs to be recomputed or not
s = size(circularmask)
sm = EBSDdata.detnumsx
if (s[0] ne sm) then begin
  d = shift(dist(sm),sm/2,sm/2)
  d[where(d le sm/2)] = 1.0
  d[where(d gt sm/2)] = 0.0
  circularmask = fltarr(EBSDdata.detnumsx, EBSDdata.detnumsy)
  dm = (EBSDdata.detnumsx - sm)/2
  circularmask[dm,0] = d
endif

if (EBSDdata.numangles gt 1000) then begin
    Core_Print,'',/blank
    Core_Print,'You are computing more than 1000 EBSPs; this will take a while...'
    Core_Print,'The program will not provide any further updates until the run has been completed.'
    Core_Print,'',/blank
endif

; next, generate the ipar and fpar parameter arrays used for call_external

; ipar(1) = 1 if GetVectorsCone detector arrays need to be computed, 0 if not (arrays will have save status)
; ipar(2) = enl%npix
; ipar(3) = enl%numsy
; ipar(4) = enl%numEbins
; ipar(5) = enl%nsx
; ipar(6) = enl%nsy
; ipar(7) = nx = (enl%nsx-1)/2
; ipar(8) = enl%npx

nipar = long(8)
ipar = lon64arr(nipar)

ipar[0] = long64(1)    ; will need to be modified 
ipar[1] = long64(EBSDdata.detnumsx)
ipar[2] = long64(EBSDdata.detnumsy)
ipar[3] = long64(EBSDdata.mcenergynumbin)
ipar[4] = long64(EBSDdata.mcimx)
ipar[5] = long64(EBSDdata.mcimy)
ipar[6] = long64((EBSDdata.mcimx-1)/2)
ipar[7] = long64(EBSDdata.mpimx)

; fpar(1) = enl%thetac
; fpar(2) = enl%sampletilt
; fpar(3) = enl%workingdistance
; fpar(4) = enl%Rin
; fpar(5) = enl%Rout
; fpar(6) = enl%sigstart
; fpar(7) = enl%sigend
; fpar(8) = enl%sigstep
; fpar(9-12) =  quaternion for requested Euler angles

nfpar = long(12)
fpar = fltarr(nfpar)

fpar[0] = EBSDdata.detthetac
fpar[1] = EBSDdata.detsampleytilt
fpar[2] = EBSDdata.detW
fpar[3] = EBSDdata.detRi
fpar[4] = EBSDdata.detRo
fpar[5] = EBSDdata.mcsigstart
fpar[6] = EBSDdata.mcsigend
fpar[7] = EBSDdata.mcsigstep

; and here is the quaternion that represents the Euler angle triplet
quaternion = Core_eu2qu( [EBSDdata.detphi1, EBSDdata.detphi, EBSDdata.detphi2] )
fpar[8:11] = quaternion[0:3]

; initialize the simulated pattern array
ECpattern = fltarr(EBSDdata.detnumsx,EBSDdata.detnumsy)

callname = 'SingleECPatternWrapper'

res = call_external(EBSDdata.EMsoftpathname+'/libEMSoftLib.dylib', callname, $
      nipar, nfpar, long(EBSDdata.mcenergynumbin), ipar[4], long(EBSDdata.mpimx), ipar[1], ipar[2], ipar, fpar, float(accum_e), mLPNH, mLPSH, ECpattern, /F_VALUE, /VERBOSE, /SHOW_ALL_OUTPUT)

if (res ne 1.0) then begin
  Core_print,'SingleECPPatternWrapper return code = '+string(res,format="(F4.1)")
end 

end

