pro lco_phot,name, phot

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'lco_phot, phot_fil, phot_str' 
    return
  endif 

 readcol, name, id, xpix, ypix, B, Bs, R, Rs, area, karea, star

 tmp = { $
         id: 0L, $
         xpix: 0., $
         ypix: 0., $
         B: 0., $
         Bs: 0., $
         R: 0., $
         Rs: 0., $
         area: 0., $
         karea: 0., $
         star: 0. $
       }

 phot = replicate(tmp, n_elements(id))

 phot.id = id
 phot.xpix = xpix
 phot.ypix = ypix
 phot.B = B
 phot.Bs = Bs
 phot.R = R
 phot.Rs = Rs
 phot.area = area
 phot.karea = karea
 phot.star = star

 return

end
