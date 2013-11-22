pro tester,fname
;

openu,1,fname,/f77

; get the file type (first four characters of the file)
ftp = bytarr(4)
readu,1,ftp
ftp = string(ftp)
print,' file type = ',ftp

; first a pair of strings of 132 characters each
dataname = bytarr(132)
readu,1,dataname
print,'dataname = ',string(dataname)

finfo = file_info(fname)
print,'file size = ',finfo.size

xtalname = bytarr(132)
readu,1,xtalname
print,'Xtalname = ->'+string(xtalname)+'<-'

; wave vector indices (3 longints)
wavek = lonarr(3)
readu,1,wavek
wv = '['+string(wavek[0],format="(I2)")+' '+ string(wavek[1],format="(I2)")+' '+ string(wavek[2],format="(I2)")+']'
print,'Wave vector = '+wv


; Bragg angle of first (ga) reflection
bragg = 0.0
readu,1,bragg
print,'Primary Bragg angle = '+string(bragg,FORMAT="(F6.3)")

; number of pixels along disk radius
nums = 0L
readu,1,nums
print,'Number of pixels along disk radius = '+string(nums,FORMAT="(I)")

; wave length
mLambda = 0.0
readu,1,mLambda
print,'Wave length = '+string(mlambda*1000.0,FORMAT="(F7.4)")

; beam convergence
thetac = 0.0
readu,1,thetac
print,'Beam convergence = '+string(thetac,FORMAT="(F6.3)")

; pixel size 
dfl = 1.0
readu,1,dfl
print,'Pixel Size [nm] = '+string(dfl,FORMAT="(F6.3)")

; reflections
numref = 0L
readu,1,numref
print,'Number of reflections = '+string(numref,FORMAT="(I4)")
hkl = lonarr(3)
qx = 0.0
qy = 0.0
indices = lonarr(3,numref)
offsets = fltarr(2,numref)
for i=0,numref-1 do begin
  readu,1,hkl
  readu,1,qx,qy
  indices[0:2,i] = hkl
  offsets[0:1,i] = [qx,qy]
endfor

; wavevectors
numk = 0L
readu,1,numk
print,'Number of wave vectors = '+string(numk,FORMAT="(I6)")
ki = 0L
kj = 0L
kperp=fltarr(2,numk)
for i=0,numk-1 do begin
  readu,1,ki,kj
  kperp[0:1,i] = [ki,kj]
endfor

; read the dimensions of the output array
datadims = lonarr(3)
readu,1,datadims
print,'Data array has dimensions '+string(datadims[0],FORMAT="(I8)")+string(datadims[1],FORMAT="(I8)")+ $
	string(datadims[2],FORMAT="(I8)")

rawdata = fltarr(datadims)
readu,1,rawdata

close,1

stop

end
