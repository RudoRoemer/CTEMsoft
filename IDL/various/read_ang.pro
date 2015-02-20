;pro read_ang
;+
; NAME: read_ang
;     
; PURPOSE: this routine reads a TSL ANG formatted file and extracts the requested data sets
;     
; CALLING SEQUENCE: read_ang,fname,euler=euler,pos=pos,iq=iq,ci=ci
;     
;
; OPTIONAL INPUT:
;
; OPTIONAL KEYWORD INPUT:
; 
; METHOD: simple I/O, with some trickier stuff to interpolate the data onto a square grid
;
; NOTES:
;
; PROCEDURES/FUNCTIONS USED:
;   
; PROJECT: IMTS/D3D ONR program
;
; REVISION HISTORY
;      Written, MDG, 7/13/06
;-
;On_error,2                    ;Return to caller

pro read_ang,fname,euler_angles,positions,image_quality,confidence_index,graphics=graphics,lskip=lskip,phase=p
;
; make sure the fname file exists
;

res = file_test(fname)
if (res ne 1) then begin 
  print,'Error: file does not exist'
  return
endif

if arg_present(lskip) then skip = lskip else skip=32

; then get the number of lines in the file
spawn,'wc -l '+fname+'>linecount'
openr,1,'linecount'
lc = long(0)
readf,1,lc
close,1
spawn,'rm linecount'
lc = lc - skip
e1 = fltarr(lc)
e2 = fltarr(lc)
e3 = fltarr(lc)
px = fltarr(lc)
py = fltarr(lc)
iq = fltarr(lc)
ci = fltarr(lc)
p = intarr(lc)

; the first skip lines are comments and can be skipped
line = ''
openr,1,fname
for i=1,skip do readf,1,line

; here we have the real data
a=0.0 & b=0.0 & c=0.0 & d=0.0 & e=0.0 & f=0.0 & g=0.0
ph = 0 & x=0 & y = 0.0
for i=0L,lc-1 do begin
  readf,1,a,b,c,d,e,f,g,ph,x,y
  e1[i] = a
  e2[i] = b
  e3[i] = c
  px[i] = d
  py[i] = e
  iq[i] = f
  ci[i] = g
  p[i] = ph
endfor
close,1

; reform the array to 2D form
mx = max(px) 
my = max(py)

dx = px[1]-px[0]
nx = fix(mx/dx)+1
dy = py[nx]-py[0]
ny = fix(my/dy)+1

; the ANG format puts 4pi in each Euler angle when the confidence index equals -1
; We change this here to -0.1 in each angle, since they must be positive anyway
q=where(ci eq -1.0,count)
if (count ne 0) then begin 
	e1[q] = 0.0
	e2[q] = 0.0
	e3[q] = 0.0
endif

euler_angles = fltarr(3,nx,ny)
euler_angles[0,0:*,0:*] = reform(e1,nx,ny) 
euler_angles[1,0:*,0:*] = reform(e2,nx,ny) 
euler_angles[2,0:*,0:*] = reform(e3,nx,ny) 

image_quality = reform(iq,nx,ny)
confidence_index = reform(ci,nx,ny)

if arg_present(phase) then p = reform(p,nx,ny)

positions = [[px],[py]]

euler = euler_angles

if keyword_set(graphics) then begin 
  window,10,xsi = nx,ysi = ny,retain=2
  print,min(euler[0,*,*],max=ma)/!dtor,ma/!dtor
  print,min(euler[1,*,*],max=ma)/!dtor,ma/!dtor
  print,min(euler[2,*,*],max=ma)/!dtor,ma/!dtor
  euler[0,*,*] /= (2.0*!pi)
  euler[1,*,*] /= (0.5*!pi)
  euler[2,*,*] /= (2.0*!pi)
  tvscl,bytscl(euler,min=0.0,max=1.0),true=1
endif

end

