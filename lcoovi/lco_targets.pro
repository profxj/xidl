pro lco_targets, name, galstr

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'lco_targets, targ_fil, targ_str' 
    return
  endif 

 readcol, name,idg,ragal,decgal,Rgal,xgal,ygal,sgal,kgal,algal

 tmp = { $
         idg: 0L, $
         ra: 0., $
         dec: 0., $
         R: 0., $
         xpix: 0., $
         ypix: 0., $
         sa: 0., $
         ka: 0., $
         al: 0. $
       }

 galstr = replicate(tmp, n_elements(idg))

 galstr.idg = idg
 galstr.ra = ragal
 galstr.dec = decgal
 galstr.R = Rgal
 galstr.xpix = xgal
 galstr.ypix = ygal
 galstr.sa = sgal
 galstr.ka = kgal
 galstr.al = algal

 return

end
