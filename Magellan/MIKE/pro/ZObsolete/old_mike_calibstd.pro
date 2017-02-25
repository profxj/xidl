;+ 
; NAME:
; mike_calibstd   
;   Version 1.1
;
; PURPOSE:
;    Create a sensitivity function given a standard star and its
;    appropriate calibration file.  
;
; CALLING SEQUENCE:
; mike_calibstd, mike, indx, HSTFIL=, CHKFIT=, ESOFIL=
;   
; INPUTS:
;   mike   -  MIKE structure
;   indx   -  Index of the Standard star in the MIKE structure to process
;
; RETURNS:
;
; OUTPUTS:
;  Sensitivity file in the 'Extract/' directory
;
; OPTIONAL KEYWORDS:
;  OUTFIL - Name of the output file (default:
;           'Extract/Sens_{mikeroot}.fits' )
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_calibstd, mike, 1, ESOFIL='fhr4469.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   04-Jan-2004 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_calibstd, mike, indx, HSTFIL=hstfil, CHKFIT=chkfit, $
                BSPLIN=bsplin, SWV=swv, SFX=sfx, ESOFIL=esofil, $
                BSET=bset, YFIT=yfit, SENS=sens, EVERYN=everyn, $
                OBJFIL=objfil, indx2=indx2, INTER=INTER, NCOEFF=ncoeff, $
                   SV_ALLFIT=sv_allfit

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_calibstd, mike, indx, HSTFIL=, FXFIL=, NCOEFF= [v1.1]'
    return
  endif 

;  Optional Keywords
  if n_elements(indx) NE 1 then stop
  if not keyword_set( EVERYN ) then everyn = 36/mike[indx].rowbin
  if not keyword_set( OBJFIL ) then objfil = mike[indx].obj_fil
  if not keyword_set( OUTFIL ) then begin
      slen = strlen(mike[indx].img_root)
      outfil = 'Extract/Sens_'+mike[indx].img_root
  endif
  if not keyword_set( NCOEFF ) then ncoeff= 7

; HST
  ;; Read file
  if keyword_set( HSTFIL ) then begin
     if strlen(HSTFIL) LT 30 then $
         hstfil = getenv('MIKE_DIR')+'/pro/Std/HSTFIL/'+hstfil
       readcol, hstfil, swv, sfx
  endif

  if keyword_set( ESOFIL ) then begin
      if strlen(ESOFIL) LT 30 then $
        esofil = getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+esofil
      readcol, esofil, swv, sfx
;
;     let's keep it at 1d-16
;      sfx = sfx*1d-16
  endif

  objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)

; Create summed array
  
  velpix = (mike[indx].side EQ 1 ? 1.50d : 2.10d) * mike[indx].rowbin
  loglam = alog10(1.0d + velpix / 299792.458d)
  wave0  = alog10(3000.0d)
  npix = 100000L
  wvtot = 10^(wave0 + dindgen(npix)*loglam)
  fxtot = dblarr(npix)
  
;; Loop
;;  for qq=0L,n_elements(objstr)-1 do begin
;;      gd = where(objstr[qq].var GT 0.,ngd)
;;      ii = where(abs(wvtot-objstr[qq].wave[gd[0]]) LT 0.01, nii)
;;      if nii NE 1 then stop
;;      ii = ii[0]
;;      np = gd[ngd-1]-gd[0] + 1L
;;      ;; Sum it in
;;      fxtot[ii:ii+np-1] = fxtot[ii:ii+np-1] + objstr[qq].fx[gd[0]:gd[0]+np-1]
;;  endfor
      
  ntrace = n_elements(objstr)
  nrow = long(median(objstr.nrow))

  if keyword_set(use_tflat) then begin
      stop
      tflat_spec_fil = string('Flats/Tflt_SPEC_', $
                              mike[indx].side EQ 1 ? 'B_' : 'R_', $
                              mike[indx].setup, '.fits', format='(a,a,i2.2,a)') 
      tflat_spec = mrdfits(tflat_spec_fil)			    
  endif else if n_elements(indx2) EQ 1 then begin
    obj2fil = mike[indx2].obj_fil
    obj2str = xmrdfits(obj2fil, 1, STRUCTYP='mikeobjstrct', /silent)
    tflat_spec = obj2str.box_fx[0:nrow-1]
  endif

  boxfx = objstr.box_fx[0:nrow-1]

;
;   2 iterations
;   in the first one, no guess for correction
;   in the second one, use the smoothing correction.
;
  nord = n_elements(objstr)
  tmp = {fitstrct}
  sv_allfit = replicate(tmp, nord)
  flg_sv = 0
  corr = 1.0
  for ll=1,3 do begin
    
;
;   First is to make a mask of large residuals in the ratio of the two
;     spectra
;
;      if keyword_set(tflat_spec) then begin
;
;          ratio = tflat_spec * 0.0
;          ratio_ivar = tflat_spec * 0.0
;          
;          pospix = where(tflat_spec GT 1 AND boxfx GT 1 AND $
;                         objstr.box_var[0:nrow-1] GT 0)
;          ratio[pospix] = alog(boxfx[pospix]) - alog(tflat_spec[pospix])
;          ratio_ivar[pospix] = 1.0/(1.0/boxfx[pospix] + 1.0/tflat_spec[pospix])
;          pixn = findgen(nrow)# replicate(1,ntrace)
;          mask = ratio_ivar GT 0 
;          ntop = n_elements(mask)-1
;          for jj=1,3 do begin
;              xy2traceset, pixn, ratio, rset, invvar=ratio_ivar*mask, $
;                           yfit=ratio_fit, ncoeff=7, upper=5, lower=5,$
;                           outmask=ratio_mask, maxrej=3, /silent, maxiter=100
;              ratio_res = convol((ratio-ratio_fit)*sqrt(ratio_ivar*mask),$
;                                 replicate(1.0,11)) / 11.0
;              flip_up = where(ratio_res LT 0.0 AND ratio_res[1:*] GE 0.0)
;              flip_dn = where(ratio_res GT 0.0 AND ratio_res[1:*] LE 0.0)
;              
;              pos_res = where(ratio_res GE 5 AND ratio_res[1:*] LT 5, npos)
;              for kk=0,npos-1 do begin 
;                  maxh = min(where(flip_dn GT pos_res[kk])) 
;                  minh = max(where(flip_up LT pos_res[kk])) 
;                  if (minh NE -1 AND maxh NE -1) then begin
;                      grow = ratio_res[pos_res[kk]]^2 
;                      mask[(flip_up[minh]-1-grow)>0:(flip_dn[maxh]+grow)<ntop] = 0 
;                  endif
;              endfor
;              
;              neg_res = where(ratio_res LE -3.0 AND ratio_res[1:*] GT -3.0, nneg)
;              for kk=0,nneg-1 do begin 
;                  maxh = min(where(flip_up GT neg_res[kk])) 
;                  minh = max(where(flip_dn LT neg_res[kk])) 
;                  if (minh NE -1 AND maxh NE -1) then begin 
;                      grow = ratio_res[neg_res[kk]]^2 
;                      mask[(flip_dn[minh]-1-grow)>0:(flip_up[maxh]+grow)<ntop] = 0 
;                  endif
;              endfor
;          endfor 
;      endif

      boxwv = objstr.box_wv[0:nrow-1]
      sb = sort(boxwv)
      
      
      if not keyword_set(INTER) then begin
          legendre_coeffs = fit2std(objstr, ncoeff=ncoeff, mask=mask, $
                                    wmin=wmin, wmax=wmax, corr=corr)
          wave = objstr.wave
          n = n_elements(wave[*,0])
          good = where(wave GT 2500.0 AND objstr.var GT 0)
          ord = good / n
          x = 2.0d*((wave)[good] - wmin[ord]) / (wmax[ord] - wmin[ord]) - 1.0

          sens = exp(total(flegendre(x, ncoeff) * $
                           transpose(legendre_coeffs[*, ord]),2))
      endif else begin
          legendre_coeffs = fltarr(max(sv_allfit.nord), nord)
          for qq=0L,nord-1 do begin
              if flg_sv then fitstr = sv_allfit[qq] else $
                fitstr = x_setfitstrct(FUNC='LEGEND', NORD=9L, LSIG=3., HSIG=3., $
                                       FLGREJ=1)
              wv = objstr[qq].box_wv[0:objstr[qq].npix-1]
              fx = objstr[qq].box_fx[0:objstr[qq].npix-1]
              var = objstr[qq].box_var[0:objstr[qq].npix-1]
              fit = x1dfit(wv, fx, sig=sqrt(var>0), $
                           /inter, fitstr=fitstr)
              if fitstr.nord GT 20 then stop 
              sv_allfit[qq] = fitstr
              delvarx, fitstr
              legendre_coeffs[0:sv_allfit[qq].nord-1
          endfor
          if flg_sv EQ 0 then flg_sv = 1
          ans = x_guinum(title='Are you finished, 1=Yes')
          if ans EQ 1 then break else stop
          wmin = sv_allfit.nrm[0]
          wmax = sv_allfit.nrm[1]
      endelse


;  plot, wave[good], sens
;  plot, wave[good], (objstr.fx)[good]/sens, ps=3, yr=[0.5,1.5]

      if ll EQ 3 then break
      
      bf = (objstr.fx)[good]/sens
      bi = sens^2/((objstr.var)[good])
      bi = bi * ((bf^2*bi) LT 3.0e4)
      
      bs = sort(wave[good])
      bset = bspline_iterfit(wave[good[bs]], bf[bs],everyn=50, $
                             invvar=bi[bs], yfit=bmodel, nord=3, outmask=outmask, $
                             upper=3, lower=3, /groupbadpix, maxrej=1,maxiter=100)
      
      
      bset_smooth = bset
      ss = 2.0
      hpix = long(3*ss)
      kernel = gauss_kernel(ss, hpix=hpix)
      
      bset_smooth.coeff = (convol([replicate(1,hpix), bset.coeff,  $
                                   replicate(1,hpix)], kernel)) $
                          [hpix:hpix+n_elements(bset.coeff)-1]
      
      linterp, swv, sfx, boxwv[*], standard
      bunsmooth     = bspline_valu(boxwv, bset)
      bsmooth = bspline_valu(boxwv, bset_smooth)
      
      bband = where(objstr.order EQ 45)
      if bband[0] NE -1 then bsmooth[*,bband] = bunsmooth[*,bband]
      
      
      denom = bunsmooth * standard * objstr[0].exp
      corr = bsmooth / (denom  + (denom EQ 0)) 
      corr = corr * (corr GT 0 AND denom GT 0) 
   endfor 

   svset = replicate({  order: -1L, $
                         wmin: 0.0, $
                         wmax: 0.0, $
                         lcoeff: dblarr(ncoeff) }, ntrace)

   svset.wmin   = wmin
   svset.wmax   = wmax
   svset.lcoeff = legendre_coeffs
   svset.order  = objstr.order

   if keyword_set( CHKFIT ) then begin
       linterp, boxwv[sb], corr[sb], wave[good], corr1
        
       x_splot, wave[good],(1000. * (objstr.fx)[good]*corr1/wave[good] > 0)<1, $
          ytwo=1000.0*sens / wave[good], /block, psym1=3, psym2=3
   endif
      
      ;; Save
   mwrfits, svset, outfil, /create 
   print, 'mike_calibstd:  Writing ', outfil


  print, 'mike_calibstd: All Done!'
  return
end
