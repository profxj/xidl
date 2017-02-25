;+ 
; NAME:
; mike_twisky
;     Version 1.0
;
; PURPOSE:
;    Stores slit profile and gradient along each order from twilight flats.
;        Takes TFLAT as the default
;
; CALLING SEQUENCE:
;   
;  mike_twisky, mike, setup, side
;
; INPUTS:
;   mike     -  ESI structure
;
; RETURNS:
;
; OUTPUTS:
;   nothing currently, it's just a test
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_twisky, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_profile_return
;
; REVISION HISTORY:
;   12-Aug-2003 Written by SB
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_twisky_work, tflat, tflativar, ordr_str, $
     maskimage=maskimage, chk=chk, jacobian=jacobian, residual=residual 
     

     sz        = size(tflat, /dimen)
     n_orders = n_elements(ordr_str)

;
;    Mapping out orders and gaps, this is kind of slow
;
     print, 'Calling mike_ordermask...', format='(a,$)'
     maskimage = mike_ordermask(sz[0], sz[1], ordr_str, trim=1.5)
     print, 'Done.'

;
;    Fit scattered light in order gaps
;
     print, 'Calling mike_fitgap...', format='(a,$)'
     gapfit = mike_fitgap(tflat, tflativar, maskimage)
     print, 'Done'

     tflat_sub = tflat - gapfit
     tflat_sub_ivar = tflativar

;     if keyword_set(jacobian) then begin
;        tflat_sub = tflat_sub / (jacobian + (jacobian EQ 0)) * (jacobian GT 0)
;        tflat_sub_ivar = tflat_sub_ivar * jacobian^2
;     endif 

;
;     Make profile image
;
    
     ncol= sz[0]
     nrow= sz[1]
     nord = n_orders
     profile_img = tflat*0.0

     for i=0,nord-1 do begin
       inorder = where(maskimage EQ ordr_str[i].order, nin)
;       print, 'Working on Order #', i, ordr_str[i].order, ' including ', nin,$
;              '  Pixels'
       if nin LT sz[1] then continue

       slit_length = ordr_str[i].rhedg[nrow/2] - ordr_str[i].lhedg[nrow/2]
       ystart      = inorder / ncol
       ordrcen    =  (ordr_str[i].lhedg[ystart] + ordr_str[i].rhedg[ystart])/2.0
       xstart     = 1.0*(inorder mod ncol)  
       slit_frac  = 2.0* (xstart - ordrcen) / slit_length
       y          = (2.0* ystart - nrow)/nrow

       profile =  mike_profile_return(xstart, ystart,  ordr_str[i])

       ns      = 25L
       pixel_step = 0.2/ slit_length
       f        = replicate(1,ns)
       ss      = slit_frac # f - replicate(pixel_step,nin) # findgen(ns) + $
                  (0.5*(ns-1)*pixel_step)

       profile_step = profile_ord(ss, y # f, $
                      ordr_str[i].profile0, ordr_str[i].profile1)

       cc           = total(profile_step * (tflat_sub[inorder] # f),1) / $
                      total(profile_step, 1)
       if NOT keyword_set(totalcc) then totalcc = cc else $
       totalcc      = [[totalcc],[cc]]

       profile_img[inorder] = profile

     endfor
     finalcc = djs_median(totalcc,2)
       plot,finalcc, /yno
    print, max(finalcc, pl), pl, format='(f,i,$)'
     a  = 0.5*(finalcc[pl-1] + finalcc[pl+1]) - finalcc[pl]
     b =  0.5*(finalcc[pl+1]-finalcc[pl-1])
     slit_shift = 0.2*((ns-1)/2-pl +1.0*b/(2.0*a))
     print, -1.0*b/(2.0*a) + pl, a
     print, 'Image shift in columns is ', slit_shift, " pixels"

;
;   For now, no shift and no jacobian
;

    tflat_norm = tflat_sub / (profile_img + (profile_img EQ 0)) $
                      * (profile_img GT 0)
    tflat_norm_ivar = tflat_sub_ivar * profile_img^2
    tflat_fit = tflat_sub*0.0
    residual   = tflat_sub*0.0

    for i=0,nord-1 do begin
       inorder = where(maskimage EQ ordr_str[i].order, nin)
       if nin LT sz[1] then continue

       slit_length = ordr_str[i].rhedg[nrow/2] - ordr_str[i].lhedg[nrow/2]
       ystart      = inorder / ncol
       ordrcen    =  (ordr_str[i].lhedg[ystart] + ordr_str[i].rhedg[ystart])/2.0
       xstart     = 1.0*(inorder mod ncol)  
       slit_frac  = 2.0* (xstart - ordrcen) / slit_length
       y          = (2.0* ystart - nrow)/nrow
       nbkpts = 1.5*(max(ystart) - min(ystart))
       ywave = mike_qw(xstart, ystart, ordr_str[i].arc_m, arc_slope=arc_slope)

       tempivar = tflat_norm_ivar[inorder] * (abs(slit_frac) LT 1.00)
      twi_set = bspline_iterfit(ywave, tflat_norm[inorder], $
               invvar=tempivar, yfit=yfit, $
                nbkpts=nbkpts, /groupbadpix, maxrej=5, /silent)

       tflat_fit[inorder] = yfit
       residual[inorder] = (tflat_norm[inorder] - yfit) *  $
                           sqrt(tflat_norm_ivar[inorder])
    endfor
    

    return, tflat_fit * profile_img
end

pro mike_twisky, mike, setup, side, chk=chk, $
            tflat_fil=tflat_fil, ordr_fil=ordr_fil, profile_fil=profile_fil, $
            jacob_fil=jacob_fil, mflat_fil=mflat_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_slitflat, mike, setup, [side]'
      return
  endif

   if not keyword_set( SIDE ) then side = [2L]
   if not keyword_set( chk  ) then chk=0

; Setup
   if setup GE 10 then c_s = strtrim(setup,2) else c_s = '0'+strtrim(setup,2) 

;  Loop on side
   for ii=0L,n_elements(side)-1 do begin
       qq = side[ii]
       ;; SIDE
       if qq EQ 1 then begin
           print, 'mike_twisky: Profiling BLUE trace flat'
           nm = 'B'
       endif else begin
           print, 'mike_twisky: Profiling RED trace flat'
           nm = 'R' 
       endelse

       ;; Order structure
       if NOT keyword_set(ordr_fil) then $
         ordr_fil = 'Flats/OStr_'+nm+'_'+c_s+'_P.fits'
       if x_chkfil(ordr_fil+'*',/silent) EQ 0 then begin $
           print, 'mike_twisky: Order structure does not exist', ordr_fil
           continue
       endif
       print, 'mike_twisky: Reading...', ordr_fil
       ordr_str   = xmrdfits(ordr_fil,1,/silent)

       jacobian_fil = 'Flats/Jacobian_'+nm+'_'+c_s+'.fits'
       if x_chkfil(jacobian_fil+'*',/silent) EQ 0 then begin $
           print, 'mike_twisky: Order structure does not exist', jacobian_fil
           jacobian = 0
       endif else begin
           print, 'mike_twisky: Reading...', jacobian_fil
           jacobian  = xmrdfits(jacobian_fil,0,/silent)
       endelse


       ;; TFLAT

       ;; do each TFLT one at a time
       flt = where(mike.side EQ qq AND mike.flg_anly NE 0 AND $
                  strtrim(mike.type,2) EQ 'TFLT' $
                  AND mike.setup EQ setup, nflt)

       if nflt EQ 0 then begin
          print, 'I found exactly 0 TFLTs on this side in this setup'
          continue
       endif

       mike_proc, mike, flt, setup=setup, side=qq  

       for iflt = 0,nflt-1 do begin
         tflat_fil = mike[flt[iflt]].img_final

         print, 'mike_twisky: Reading...', tflat_fil 
         image = xmrdfits(tflat_fil, 0, hdr)
         invvar = xmrdfits(tflat_fil, 1)
         
         fit = mike_twisky_work(image, invvar, ordr_str, jacobian=jacobian)
         
       endfor  

   endfor

   ;; All done
   print, 'mike_twisky: All done'
   return

end


