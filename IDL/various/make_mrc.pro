@MRCstuff

pro make_mrc,sinogram,outname
;
; take a sinogram and dump it to a .mrc file.
;
dims = size(sinogram,/dimensions)

nth = dims[0]
th = findgen(nth) - float(nth/2)

mi = min(sinogram,max=ma)
print,'extreme values = ',mi,ma

sinogram = uint(sinogram)
mi = min(sinogram,max=ma)
print,'extreme uint values = ',mi,ma
me = mean(sinogram)

; and then convert these slices into a single .mrc file
MRCinit,M,F,/create

; set the M variable
M.nx = long(dims[0])
M.ny = long(dims[1])
M.mx = M.nx
M.my = M.ny
M.nz = long(nth)
M.mz = M.nz
M.amin = float(mi)
M.amax = float(ma)
M.amean = float(me)
M.xlen = M.nx
M.ylen = M.ny
M.zlen = M.nz


; set the tilt angles (meaningless in this case, but needs to be done)
for i=0,nth-1 do begin
  F[i].b_tilt = th[i]
  F[i].defocus = 0.0
  F[i].pixelsize = 1.0
  F[i].magnification = 10000.0
  F[i].voltage = 200000
  F[i].mean_int = mean(sinogram[*,*,i])
endfor

; and write the file
spawn,'/usr/bin/touch '+outname
openu,1,outname
writeu,1,M,F
for i=0,nth-1 do writeu,1,reform(sinogram[0:*,0:*,i])
close,1

end
