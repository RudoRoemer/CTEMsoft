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
; EMsoft:Efitinit.pro
;--------------------------------------------------------------------------
;
; PROGRAM: Efitinit.pro
;
;> @author Marc De Graef, Carnegie Mellon University
;
;> @brief initialize all arrays for EBSD pattern sumulation (mostly from EBSDmod.f90)
;
;> @details 10/12/15 a new GUI to interactively determine the best fit parameters
;> for an EBSD pattern; this routine is an IDL copy of the EMEBSD initialization code
;
;> @date 10/12/15 MDG 1.0 first attempt at a user-friendly interface
;> @date 10/13/15 MDG 1.1 updated with corrected structure names
;> @date 10/17/15 MDG 1.2 removed array computations; moved to callable f90 code in EMdymod.f90
;--------------------------------------------------------------------------
pro Efitinit,dummy

common Efit_widget_common, Efitwidget_s
common Efit_data_common, Efitdata

common EBSD_EMsoft, MCxtalname, MCmode, nsx, nsy, EkeV, Ehistmin, Ebinsize, depthmax, depthstep, MCsig, MComega, $
                    numEbins, numzbins, accum_e, accum_z, Masterenergyfile, npx, npy, nnE, numset, mLPNH, mLPSH, Masterxtalname, expEBSDpattern, EBSDpattern


Core_Print,'Loading data files'
Core_Print,'   Opening Master pattern file'

; first read the master file
file_id = H5F_OPEN(strtrim(Efitdata.pathname+'/'+Efitdata.mpfilename,2))

group1_id = H5G_OPEN(file_id,'NMLparameters')
group2_id = H5G_OPEN(group1_id,'EBSDMasterNameList')

dset_id = H5D_OPEN(group2_id,'energyfile')
energyfile = H5D_READ(dset_id)
Efitdata.energyfilename = energyfile
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'npx')
npx = H5D_READ(dset_id)
npx = npx[0]
npy = npx
H5D_close,dset_id

H5G_close,group2_id
H5G_close,group1_id

group2_id = H5G_OPEN(file_id,'EMData')

dset_id = H5D_OPEN(group2_id,'numEbins')
nnE = H5D_READ(dset_id)
nnE = nnE[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'numset')
numset = H5D_READ(dset_id)
numset = numset[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'mLPNH')
mLPNH = H5D_READ(dset_id)
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'mLPSH')
mLPSH = H5D_READ(dset_id)
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'xtalname')
Masterxtalname = H5D_READ(dset_id)
H5D_close,dset_id

H5G_close,group2_id
H5F_close,file_id

Core_Print,'   -> Master Pattern file read'
Core_Print,'   Opening Monte Carlo file'

; then the Monte Carlo file
file_id = H5F_OPEN(strtrim(Efitdata.EMdatapathname+'/'+Efitdata.energyfilename,2))

group1_id = H5G_OPEN(file_id,'NMLparameters')
group2_id = H5G_OPEN(group1_id,'MCCLNameList')

dset_id = H5D_OPEN(group2_id,'xtalname')
MCxtalname = H5D_READ(dset_id)
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'MCmode')
MCmode = H5D_READ(dset_id)
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'numsx')
nsx = H5D_READ(dset_id)
nsx = nsx[0]
nsx = (nsx-1)/2
nsy = nsx
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'EkeV')
EkeV = H5D_READ(dset_id)
EkeV = EkeV[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'Ehistmin')
Ehistmin = H5D_READ(dset_id)
Ehistmin = Ehistmin[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'Ebinsize')
Ebinsize = H5D_READ(dset_id)
Ebinsize = Ebinsize[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'depthmax')
depthmax = H5D_READ(dset_id)
depthmax = depthmax[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'depthstep')
depthstep = H5D_READ(dset_id)
depthstep = depthstep[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'sig')
MCsig = H5D_READ(dset_id)
MCsig = MCsig[0]
Efitdata.detMCsig = MCsig
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'omega')
MComega = H5D_READ(dset_id)
MComega = MComega[0]
H5D_close,dset_id

H5G_close,group2_id
H5G_close,group1_id

group2_id = H5G_OPEN(file_id,'EMData')

dset_id = H5D_OPEN(group2_id,'numEbins')
numEbins = H5D_READ(dset_id)
numEbins = numEbins[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'numzbins')
numzbins = H5D_READ(dset_id)
numzbins = numzbins[0]
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'accum_e')
accum_e = H5D_READ(dset_id)
H5D_close,dset_id

dset_id = H5D_OPEN(group2_id,'accum_z')
accum_z = H5D_READ(dset_id)
H5D_close,dset_id

H5G_close,group2_id
H5F_close,file_id

Core_Print,'   -> Monte Carlo file read'



end
