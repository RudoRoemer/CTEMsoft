
function read_mrc_image,fname,M,i,signed=signed
;
; read a single image from a previously opened .mrc file (FEI format only)
;

openu,1,fname
; compute the offset 
offset = 1024L + M.next + long(i)*M.nx*M.ny*2L
if keyword_set(signed) then q=assoc(1,intarr(M.nx,M.ny),offset) else q=assoc(1,uintarr(M.nx,M.ny),offset)
image = q[0]
close,1

return,image
end


pro MRCinit,MRCHeader,FEIHeaders,create=create,fname = fname
;
; initialize the FEIHeaders and MRCHeader structures; these must 
; be filled by the main program, or they are read from the file if 
; fname is present.
;
; written by MDG, 12/27/11 for EIC project, based on M. Jackson's header file
;

;/**
; * @brief this was take from http://www.biochem.mpg.de/doc_tom/index.html using the
; * tom_mrcfeistack2emseries code.
; */
FEIHeaders = replicate({FEIstruct,  $
	a_tilt : float(0), $
	b_tilt : float(0), $
	x_stage : float(0), $
	y_stage : float(0), $
	z_stage : float(0), $
	x_shift : float(0), $
	y_shift : float(0), $
	defocus : float(0), $
	exp_time : float(0), $
	mean_int : float(0), $
	tiltaxis : float(0), $
	pixelsize : float(0), $
	magnification : float(0), $
	voltage : float(0), $
	unused : string(' ',format='(A72)') $
},1024)

;/**
; * @brief This spec was taken from http://bio3d.colorado.edu/imod/doc/mrc_format.txt
; * and we are going to assume an IMOD version of 2.6.20 and above:
; */
MRCHeader = {MRCstruct, $
	nx : long(0), $  ; number of columns
	ny : long(0), $  ; number of rows
	nz : long(0), $ ; number of sections
	mode : long(1), $ ; type of image pixel 
	nxstart : long(0), $ ; starting point of subimage
	nystart : long(0), $ 
	nzstart : long(0), $ 
	mx : long(0), $ ; grid size in x
	my : long(0), $ ; grid size in y
	mz : long(0), $ ; grid size in z
	xlen : float(0), $ ; cell size; pixel spacing = xlen/mx, ylen/my, zlen/mz
	ylen : float(0), $
	zlen : float(0), $
	alpha : float(90), $ ; cell angles - ignored by IMOD
	beta : float(90), $
	gamma: float(90), $
	mapc : long(1), $ ; map column  1=x,2=y,3=z
	mapr : long(2), $ ; map row     1=x,2=y,3=z
	maps : long(3), $ ; map section 1=x,2=y,3=z
	amin : float(0), $ ; minimum pixel value (needs to be set for proper scaling of data)
	amax : float(0), $ ; maximum pixel value
	amean : float(0), $ ; mean pixel value
	ispg : fix(0), $ ; space group number
	nsymbt : fix(0), $ ; NOT SURE
	next : long(131072), $ ; number of bytes in extended header (1024 * 128 for FEI)
	creatid : fix(0), $ ; used to be an ID number, is 0 as of IMOD 4.2.23
	extra_data : string(' ',format='(A30)'), $ ; not used, first two bytes should be 0
	nint : fix(0), $ ; number of bytes per section (SerialEM interpretation)
	nreal : fix(32), $ ; bit flags for short data type
	extra_data_2 : string(' ',format='(A20)'), $ ; not used
	imodStamp : long(0), $ ; 
	imodFlags : long(0), $ ;
	idtype : fix(0), $ ; ( 0 = mono, 1 = tilt, 2 = tilts, 3 = lina, 4 = lins)
	lens : fix(0), $
	nd1 : fix(0), $ ; for idtype = 1, nd1 = axis (1, 2, or 3)
	nd2 : fix(0), $ 
	vd1 : fix(0), $ ; vd1 = 100. * tilt increment
	vd2 : fix(0), $ ; vd2 = 100. * starting angle
	tiltangles : fltarr(6), $ ; 0,1,2 = original:  3,4,5 = current
	xorg : float(0), $ ; origin of image
	yorg : float(0), $
	zorg : float(0), $
	cmap : 'MAP ', $
	stamp : '    ', $ ; First two bytes have 17 and 17 for big-endian or 68 and 65 for little-endian
	rms : float(0), $ ; RMS deviation of densities from mean density
	nLabels : long(0), $ ; Number of labels with useful data
	labels : string(' ',format='(A800)') $ ; 10 labels of 80 characters each
}

; set the big-endian / little-endian variable stamp

if keyword_set(create) then begin
 MRCHeader.nlabels = long(1)
 s = 'EIC project, MDG, Copyright 2011'
 MRCHeader.labels = s + string(' ',format="(A768)")
endif

if arg_present(fname) then begin
  openu,1,fname
  readu,1,MRCHeader,FEIHeaders
  close,1
endif

end

