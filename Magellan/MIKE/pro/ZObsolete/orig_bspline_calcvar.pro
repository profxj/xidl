
;;
;;  Try to calculate inverse covariance matrix of pixel values
;;
;;   If one tests a model against bspline_valu with chi^2:
;;
;;  $ \chi^2 = (Flux - Model) # (invert(covariance)) # (Flux-Model) $
;;
;;  This routine will return, (invert(covariance))  [Npix, Npix]
;;

function  bspline_calcvar, action, loweraction, upperaction, alpha

   nn = n_elements(loweraction)

   npix     = (size(action))[1]

   bw       = (size(alpha))[1] 
   ncoeff   = (size(alpha))[2] 

   al = alpha[*, lindgen(ncoeff/2)*2]
   al = al[lindgen(bw/2)*2, *]

   bw       = (size(al))[1] 
   ncoeff   = (size(al))[2]  - bw


   infl = where(upperaction - loweraction GE 0, nin)
   tot_isig = fltarr(npix)
   tal = (total(al,1))[0:ncoeff-1]
   for i=0,nin-1 do begin & $
      l1 = loweraction[infl[i]] & $
      u1 = upperaction[infl[i]] & $
      tot_isig[l1:u1] = action[l1:u1,*] # tal[infl[i]:infl[i]+bw-1] & $
   endfor

;   tot_isig = full_action # ((total(al,1))[0:ncoeff-1] )
   tot_covar = 1./(tot_isig^2  + (tot_isig EQ 0)) * (tot_isig GT 0)
return, tot_covar 

;   print, '  First Invert alpha ...', format ='(a,$)'
;   inv_al = identity(ncoeff+bw)

;;; first round
;   spot = lindgen(bw-1) + 1
;   for j=0,ncoeff-1 do begin 
;     inv_al[*,j] = inv_al[*,j] / al[0,j] 
;     inv_al[*,j+spot] = inv_al[*,j+spot] - inv_al[*,j] # al[spot,j] 
;   endfor

;;   print, 'Inverse is done'
;;;; second round

;   spot = bw-1 - lindgen(bw-1)
;   for j=ncoeff-1,0,-1 do $
;     inv_al[*,j] = (inv_al[*,j] - al[spot,j] ## inv_al[*,j+spot])/al[0,j]

;   inv_al = inv_al[0:ncoeff-1, 0:ncoeff-1]
;   print, '  Then sum up variance...', format ='(a,$)'

;   data_covar    = fltarr(npix,npix)
;   data_infl    = fltarr(npix,npix)
;   full_action   = [[action], [replicate(0,npix,ncoeff-bw)]]
;   for i=0,nin-1 do begin & $
;      l1 = loweraction[infl[i]] & $
;      u1 = upperaction[infl[i]] & $
;      full_action[l1:u1, * ] = shift(full_action[l1:u1, * ], 0, infl[i])  & $
;   endfor
 
;   for i=0,nin-bw do begin
;       l1 = loweraction[infl[i]]
;       u1 = upperaction[infl[i+bw-1]]
;       p_l = infl[i]
;       p_u = (infl[i+bw-1]+bw-1) < (ncoeff-1)
;
;       a1 = full_action[l1: u1, p_l:p_u]
;
;       data_covar[l1:u1, l1:u1] = data_covar[l1:u1, l1:u1] + $
;           a1 # inv_al[p_l:p_u, p_l:p_u] # transpose(a1)
;       data_infl[l1:u1, l1:u1] = data_infl[l1:u1, l1:u1] + $
;           a1 # transpose(a1)
;   endfor

;   tot_covar = total(data_covar,1) 
;   tot_infl = total(data_infl,1) 
   
end

