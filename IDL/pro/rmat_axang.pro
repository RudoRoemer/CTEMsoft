function rmat_axang,axis,angle,hand
  ;;This function returns a rotation matrix for rotating by an angle theta (in degrees) around an axis defined by a unit vector.
  ;; Defaults to right handed but if you want left handed say Rmat=rmat_axang([x,y,z],angle,'left')
  
  ;;For proper use if you had a vector [1,1,0] that you wanted to rotate by 45deg around
  ;;the z-axis ([0,0,1]), then you would enter
  ;; [1,1,0]##rmat_axang([0,0,1],45)
  ;; which will return 1.41421
  if n_elements(hand) eq 0 then hand = 'right'
  
;; Unitize the vector.
  axis=axis/norm(axis)
  
  ux=axis[0]
  uy=axis[1]
  uz=axis[2]
  angle=angle*!pi/180.0
  sth=sin(angle)
  cth=cos(angle)
  
  if hand eq 'left' then begin
    R=[[    cth + ux^2*(1-cth),  ux*uy*(1-cth)-uz*sth,  ux*uz*(1-cth)+uy*sth],$
       [  uy*ux*(1-cth)+uz*sth,      cth+uy^2*(1-cth),  uy*uz*(1-cth)-ux*sth],$
       [  uz*ux*(1-cth)-uy*sth,  uz*uy*(1-cth)+ux*sth,      cth+uz^2*(1-cth)]]
  endif
  if hand eq 'right' then begin
    R=[[    cth + ux^2*(1-cth),  ux*uy*(1-cth)+uz*sth,  ux*uz*(1-cth)-uy*sth],$
       [  uy*ux*(1-cth)-uz*sth,      cth+uy^2*(1-cth),  uy*uz*(1-cth)+ux*sth],$
       [  uz*ux*(1-cth)+uy*sth,  uz*uy*(1-cth)-ux*sth,      cth+uz^2*(1-cth)]]
  endif
     
;     stop
  return,R

end