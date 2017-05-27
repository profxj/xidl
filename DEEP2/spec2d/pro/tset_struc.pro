; Finkbeiner
; Define traceset structure for wavelength solution
function tset_struc, func, ncoeff, ntrace

   if (NOT keyword_set(ntrace)) then ntrace = 1

   tset = $      
    { func    :    func               , $
      xmin    :    0.0d               , $
      xmax    :    0.0d               , $
      coeff   :    dblarr(ncoeff, ntrace) $
    }

   return, tset
end
