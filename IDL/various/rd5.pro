@make_mrc
@read_ang
@prep_mrc_file

pro rd5,dummy

; version 5 of a program to create an mrc file of a 3D quatenrion stereographic projection 


; Titanium data set UCSB

fname = 'Ti64_2D_EBSD_Set2_Strain1/strain1_1.ang'

; this is the number of lines in this file
; and the first 150 are not part of the Euler angle list
lc = 11829198L
skip = 150
lc = lc - skip

; when you first run this program, comment this line
goto,skipit

; first we determine how many alpha and beta phase entries there are 
nalpha = 0L
nbeta = 0L

openr,1,fname
line=''
for i=1,skip do begin
  readf,1,line
endfor

for i=1L,lc do begin
  readf,1,line
  z = strsplit(line,' ',/extract)
  if (fix(z[7]) eq 1) then nalpha += 1 else nbeta += 1
endfor
close,1

;skipit:
; once we've doen that we can simply hardcode those values
nalpha = 11482427L
nbeta = 346621L

; then we create the two euler angle files
eu_alpha = fltarr(3,nalpha)
eu_beta = fltarr(3,nbeta)

openr,1,fname
line=''
for i=1,skip do begin
  readf,1,line
endfor
;
cnta = 0L
cntb = 0L
for i=1L,lc do begin
  readf,1,line
  z = strsplit(line,' ',/extract)
  if (fix(z[7]) eq 1) then begin
    eu_alpha(0:2,cnta) = float(z[0:2])
    cnta += 1L
  end else begin
    eu_beta(0:2,cntb) = float(z[0:2])
    cntb += 1L
  end
endfor
close,1


spawn,'touch Tiequi_alpha.data'
openu,1,'Tiequi_alpha.data',/f77
writeu,1,nalpha
writeu,1,eu_alpha
close,1
print,'file Tiequi_alpha.data created with ',nalpha,' data points'


spawn,'touch Tiequi_beta.data'
openu,1,'Tiequi_beta.data',/f77
writeu,1,nbeta 
writeu,1,eu_beta
close,1
print,'file Tiequi_alpha.data created with ',nbeta,' data points'

skipit:

; and here we create a few .mrc files with and without the FZ limitation
; a value of 128 for npix is reasonable; larger files take too long, smaller doesn't
; have enough resolution
npix = 128L
np = 2L

print,'calling f90 conversion program'

cube = prep_mrc_file('Tiequi_alpha.data','cube.raw','stereogramCubicFZ.mrc',32,'DrawFZ',np,npix)
tvscl,total(cube,1)

cube = prep_mrc_file('Tiequi_alpha.data','cube.raw','stereogramHexagonalFZ.mrc',24,'DrawFZ',np,npix)
tvscl,total(cube,1)

cube = prep_mrc_file('Tiequi_alpha.data','cube.raw','stereogramTiequi_alpha.mrc',24,'FullSP',np,npix)
tvscl,total(cube,1)
cube = prep_mrc_file('Tiequi_alpha.data','cube.raw','stereogramTiequi_alphaFZ.mrc',24,'FZonly',np,npix)
tvscl,total(cube,1)
cube = prep_mrc_file('Tiequi_beta.data','cube.raw','stereogramTiequi_beta.mrc',32,'FullSP',np,npix)
tvscl,total(cube,1)
cube = prep_mrc_file('Tiequi_beta.data','cube.raw','stereogramTiequi_betaFZ.mrc',32,'FZonly',np,npix)
tvscl,total(cube,1)

end
