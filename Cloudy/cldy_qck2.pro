;+ 
; NAME:
; cldy_qck2
;   Version 1.0
;
; PURPOSE:
;   Calculates a quick 2-phase model given a set of ions
;   and the ratio of these ions (and the error in the ratio)
;   and finally the two elements in the Cloudy grid which
;   are to give the input ratios.
;
; CALLING SEQUENCE:
;   
; cldy_qck2, grid, ratio, sig, ion1, ion2
;
; INPUTS:
;   grid  - CLOUDY grid
;   ratio -  Observed ionic ratios (array)
;   sig  -  Error on the ratios
;   ion1 - Z, i for ion1
;   ion2 - Z, i for ion2 
;   model - Two element array of the Cloudy grid
;
; RETURNS:
;   
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   NHILMT - 2-element array giving NHI min,max
;   FEHLMT - 2-element array giving FeH min,max
;
; OPTIONAL OUTPUTS:
;  CHISQ - chisq
;
; COMMENTS:
;
; EXAMPLES:
;   cldy_qck2, grid, [0.4], [0.1], [ 14,3 ], [14,2], model
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro cldy_qck2, grid, ratio, sig, ion1, ion2, model

;
  if  N_params() LT 6  then begin 
      print, 'Syntax - ' +$
        'cldy_qck2, grid, ratio, sig, ion1, ion2, model'
      return
  endif 

; Optional keywords

  npnt = n_elements(ratio)
  Aval = dblarr(npnt)
  sigA = dblarr(npnt)
  val = dblarr(npnt)


; Weight is inverse variance 

   chk = where(sig EQ 0.0, cnt)
   if cnt GT 0 then begin
       print, 'Cant have zero error!'
       return
   endif

;  Calculate the Ratios
;   lratio = 10^ratio
;   lsig = sig*alog(10.)*lratio
   lgrid = fltarr(2, npnt)
   for i=0,1 do begin
       for k=0,npnt-1 do $
         lgrid[i,k] = grid[model[i]].X[ion1[0,k],ion1[1,k]]-$
                   grid[model[i]].X[ion2[0,k],ion2[1,k]]
   endfor
                   

;  Calculate weighted answer

   for k=0,npnt-1 do begin
       denom = lgrid[0,k] - lgrid[1,k]
       if(denom EQ 0.) then continue
;          Calculate A, sigA
       Aval[k] = (ratio[k]-lgrid[1,k])/denom
       sigA[k] = sig[k]/abs(denom)
   endfor

   weight = 1./sigA^2
   ans = total(Aval*weight)/total(weight)
   sigans = sqrt(1./total(weight))

;  Chisq
   chisq = 0.d
   for k=0,npnt-1 do begin
       chisq = chisq + ((ratio[k] - ans*lgrid[0,k] - $
                           (1.-ans)*lgrid[1,k])/sig[k])^2
   endfor

;  Values
   for k=0,npnt-1 do $
     val[k] = ans*lgrid[0,k] + (1.-ans)*lgrid[1,k]

;  Output
   print, grid[model[0]].NHI, grid[model[0]].U, $
     grid[model[1]].NHI, grid[model[1]].U

   print, 'Ans: ', ans
   print, val, ratio
   print, 'Chisq: ', chisq
; Delete

   delvarx, lgrid

   return
end

