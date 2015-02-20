function prep_mrc_file,eulerfile,cubefile,mrc_file,pgnum,FZmode,np,npix,smooth=smooth

; creates a namelist file for the CTEMquatRC program, executes it, load
; the output and converts that to a .mrc file.  Then the result array is 
; returned to the calling program.
;
; written by MDG, 02/19/15
;
nump = 2L*(npix + np + 1)+1L

; create namelist file
openw,1,'Stereogram.nml'
printf,1,'&Stereogram'
printf,1,'eulerdatafile = '''+eulerfile+''''
printf,1,'cubefile = '''+cubefile+''''
printf,1,'pgnum = '+string(pgnum,format="(I2)")
printf,1,'FZmode = '''+FZmode+''''
printf,1,'verbose = ''y'''
printf,1,'np = '+string(np,format="(I2)")
printf,1,'npix = '+string(npix,format="(I4)")
printf,1,'/'
close,1
cmd = '/users/mdg/Files/OSU/Software/CTEMsoft2013/Build/Bin/CTEMquatRC'
print,cmd
spawn,cmd

openu,1,cubefile,/f77
cube = fltarr(nump,nump,nump)
readu,1,cube
close,1

if (keyword_set(smooth)) then cube = smooth(cube,3)

make_mrc,cube,mrc_file
print,'data stored in mrc format ',mrc_file

return,cube
end
