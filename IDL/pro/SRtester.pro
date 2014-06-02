pro SRtester,fname
;
; test program to read the output from the CTEMSRdefect program before the 
; GUI is completed...

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
wavek = fltarr(3)
readu,1,wavek
wv = '['+string(wavek[0],format="(F5.2)")+' '+ string(wavek[1],format="(F5.2)")+' '+ string(wavek[2],format="(F5.2)")+']'
print,'Foil normal = '+wv
print,wavek

; length of fundamental reflection
leng = 0.0
readu,1,leng
print,' |g| = '+string(leng,FORMAT="(F6.3)")

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

; # reflections
numref = 0L
readu,1,numref
print,'Number of reflections = '+string(numref,FORMAT="(I4)")

; read the dimensions of the output array
datadims = lonarr(4)
print,'datadims=',datadims
readu,1,datadims
print,'Data array has dimensions '+string(datadims[0],FORMAT="(I8)")+string(datadims[1],FORMAT="(I8)")+ $
	string(datadims[2],FORMAT="(I8)")+string(datadims[3],FORMAT="(I8)")

rawdata = fltarr(datadims)
readu,1,rawdata

close,1

stop

end
