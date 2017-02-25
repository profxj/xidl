;+ 
; NAME:
; mike_teststd
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

pro mike_teststd, mike, indx, HSTFIL=hstfil, CHKFIT=chkfit, $
                BSPLIN=bsplin, SWV=swv, SFX=sfx, ESOFIL=esofil, $
                BSET=bset, YFIT=yfit, SENS=sens, EVERYN=everyn, OBJFIL=objfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_calibstd, mike, indx, HSTFIL=, FXFIL= [v1.0]'
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

; HST
  ;; Read file
  if keyword_set( HSTFIL ) then begin
      hst = xmrdfits(hstfil, 1, /silent)
      swv = hst.wavelength
      sfx = hst.flux
  endif

  if keyword_set( ESOFIL ) then begin
      if strlen(ESOFIL) LT 30 then $
        esofil = getenv('MIKE_DIR')+'/pro/Std/ESOFIL/'+esofil
      readcol, esofil, swv, sfx, sjy, spix
  endif

  objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)

; Create summed array
  
  velpix = (mike[indx].side EQ 1 ? 1.50d : 2.10d) * mike[indx].rowbin
  loglam = alog10(1.0d + velpix / 299792.458d)
  wave0  = alog10(3000.0d)
  npix = 100000L
  wvtot = 10^(wave0 + dindgen(npix)*loglam)
  fxtot = dblarr(npix)
  contrib = lonarr(npix)
  frac  = objstr.fx * 0.0
  frac_fit  = objstr.fx * 0.0
  norm_fx  = objstr.fx * 0.0 
  norm_var  = objstr.var * 0.0 

  ;; Loop
  nordr= n_elements(objstr)
  for qq=0L,nordr - 1 do begin
      ii = where(abs(wvtot-objstr[qq].wave[0]) LT 0.01, nii)
      if nii NE 1 then stop
      ii = ii[0]
      np = objstr[qq].npix
      spot = lindgen(np) 
      ;; Sum it in

      pixsize = (objstr[qq].box_wv[1:*] - objstr[qq].box_wv)
      mid_wave = 0.5*(objstr[qq].box_wv[1:*] + objstr[qq].box_wv)

      good = where(objstr[qq].box_var GT 0 AND objstr[qq].box_var[1:*] GT 0)
      wave_size = abs(poly(objstr[qq].wave[spot], ladfit(mid_wave[good], pixsize[good])))

      norm_fx[spot,qq] = objstr[qq].fx[spot] / (wave_size + (wave_size EQ 0))
      norm_var[spot,qq] = objstr[qq].var[spot] / (wave_size + (wave_size EQ 0))^2 * (wave_size GT 0)

      fxtot[ii + spot] = fxtot[ii+spot] + norm_fx[spot,qq] * (norm_var[spot,qq] GT 0)
      contrib[ii + spot] = contrib[ii+spot] + (norm_var[spot,qq] GT 0)
  endfor
      
  for qq=0L,nordr - 1 do begin
      ii = where(abs(wvtot-objstr[qq].wave[0]) LT 0.01, nii)
      if nii NE 1 then stop
      ii = ii[0]
      np = objstr[qq].npix
      spot = lindgen(np) 
      keep = where(fxtot[ii+spot] GT 0 AND contrib[ii+spot] EQ 2 AND norm_var[spot,qq] GT 0, nkeep)
      ;; Sum it in
      frac[spot[keep], qq] = norm_fx[spot[keep], qq] / (fxtot[ii+spot])[keep]
      y = frac[spot[keep], qq]
      bset = bspline_iterfit(spot[keep], y, nbkpts=8, invvar=(y GT 0))
      frac_fit[spot,qq] = bspline_valu(spot, bset)
  endfor
    

  ;; Find half-flux points....

  lefthalf = fltarr(nordr)
  righthalf = fltarr(nordr)

  for qq=1,nordr-2 do begin

    np = objstr[qq].npix
    lefthalf[qq] = max(where(frac_fit[0:np/2,qq] LT 0.5))
    righthalf[qq] = min(where(frac_fit[np/2:*,qq] LT 0.5)) + np/2

  endfor

  x = [[objstr.order],[objstr.order]]
  y =  [[lefthalf],[righthalf]]
  xy2traceset, x, y, halfset, ncoeff=3, yfit=half_fit, invvar=(y GT 0), upper=20, lower=20, outmask=outmask
  print, sqrt(total((y-half_fit)^2*outmask,1) / total(outmask, 1))


  ;; Now do a 2-d bspline fit to good fractional points, with basis set by half-points.

  half = y + (half_fit - y)*(y EQ 0)
  np = max(objstr.npix)
  x2 = objstr.order ## replicate(1,np)
  yb  = frac[0:np-1,*]

  range = half[*,1] - half[*,0]
  midpoint = 0.5*(half[*,1] + half[*,0])
  xb = (findgen(np) # replicate(1,nordr)  - midpoint ## replicate(1,np)) / (range ## replicate(1,np))

  ib = (total(outmask,2) ## replicate(1,np)) * (yb NE 0) * (abs(xb) LT 1.2)
;  good = where(ib GT 0)
;  yb_log = alog(yb[good])

;  bset_final = bspline_iterfit(xb[good], yb_log, bkpt=bkpt, invvar = ib[good]*yb[good], upper=1, lower=1, $
;                                         yfit=yfit_log, x2=x2[good], npoly=2)
;  frac_final = exp(bspline_valu(xb, bset_final, x2=x2))

;  

  if mike[indx].side EQ 1 then begin
    bkpt = [-1.0, -0.8, -0.7, -0.6, -0.5, -0.35, -0.2, 0.0, 0.2, 0.35, 0.5, 0.6, 0.7, 0.8, 1.0]
  endif else begin
    bkpt = [-0.8, -0.7, -0.6, -0.5, -0.3, 0.0, 0.3, 0.5, 0.6, 0.7, 0.8]
  endelse

;  bkpt = 0.95*(findgen(15)-7)/7.
;  xb_fix = [xb,(2.2*findgen(2)-1.1) # replicate(1,nordr)]
;  yb_fix = [yb, fltarr(2,nordr)]
;  ib_fix = [ib, fltarr(2,nordr)+10.]
;  x2_fix = [x2, objstr.order ## replicate(1,2)]


  bset_final = bspline_iterfit(xb, yb, bkpt=bkpt, invvar=ib, x2=x2, npoly=3, maxdev=0.3)
  frac_final = bspline_valu(xb, bset_final, x2=x2)
;  bset_final = bspline_iterfit(xb_fix, yb_fix, bkpt=bkpt, invvar=ib_fix, x2=x2_fix, npoly=2, maxdev=0.3)
;  frac_final = frac_final * (abs(xb) LT 1.0)


 bset_1 = bspline_iterfit(xb, yb, nbkpt=15, invvar=ib, yfit=frac_1, maxdev=0.3)

  xs=sort(xb)
  bset_action = bspline_action(xb[xs], bset_1, lower=lower, upper=upper)

  bset_action = [[bset_action], [fltarr(n_elements(xs),n_elements(upper)-1)]]
  for i=1,n_elements(lower)-1 do bset_action[lower[i]:upper[i],*] = shift(bset_action[lower[i]:upper[i],*],0,i)
stop 

  
; Loop on orders

  for qq=0L,n_elements(objstr)-1 do begin
      ;; 
      npix = objstr[qq].npix
      wave = objstr[qq].wave[0:npix-1]
      flux = objstr[qq].fx[0:npix-1]
      sig = sqrt(objstr[qq].var[0:npix-1])

      ;; Dwv
      npix = objstr[qq].npix
      if wave[1] GT wave[0] then begin
          dwv = (shift(wave,-1)-shift(wave, 1))/2.
          dwv[0] = wave[1]-wave[0]
          dwv[npix-1] = wave[npix-1]-wave[npix-2]
      endif else begin
          dwv = (shift(wave,1)-shift(wave, -1))/2.
          dwv[0] = wave[0]-wave[1]
          dwv[npix-1] = wave[npix-2]-wave[npix-1]
      endelse

      ;; BSpline
      if keyword_set( BSPLIN ) then begin
          bset = bspline_iterfit(swv, sfx, yfit=yfit, everyn=everyn, /silent)
          if keyword_set( CHKFIT ) then begin
              x_splot, swv, sfx, ytwo=yfit, /block
              stop
          endif
          ;; Calculate at wavelength
          sens = bspline_valu(wave, bset) / (flux / mike[indx].exp / dwv)
          
      endif else begin ;; SPLINE
          ;; Calculate at wavelength
          splin = spl_init(swv, sfx, /double)
          std = spl_interp(swv, sfx, splin, wave)
          
          ;; Sens function
          inv_sens = (flux / mike[indx].exp / dwv) / std
          inv_sens_ivar = (std * mike[indx].exp * dwv)^2/(sig^2 + (sig EQ 0)) * (sig GT 0)
          objstr[qq].flux[0:npix-1] = inv_sens
          objstr[qq].sig[0:npix-1] = inv_sens_ivar
if qq EQ 5 then stop
      endelse

      continue
stop

      ;; Mask
      good = where( sig GT 0., ngood)
      if ngood EQ 0 then stop
      srt = sort(wave)
      wave = wave[srt]
      sens = sens[srt]
      inv_sens = inv_sens[srt]
      inv_sens_ivar = inv_sens_ivar[srt]
  
      ;; BSpline
      ivar = fltarr(n_elements(sens))
      ivar[*] = 1.

      ;; Zero out neg values (for bluest orders)
      a = where(sens LE 0., na)
      if na NE 0 then ivar[a] = 0.

      ;; Abs
;      bd = where((wave GT 6864. AND wave LT 6920.) OR $
;                 (wave GT 7585. AND wave LT 7684.), nbd)
;      if nbd NE 0 then ivar[bd] = -1.
      
      bset = bspline_iterfit(wave, inv_sens, yfit=yfit, everyn=everyn, /groupbadpix, maxrej=5, $
                         lower=5., upper=5., grow=21, invvar=inv_sens_ivar, /silent, outmask=outmask)

      next_mask = smooth(outmask, 61)
      bset = bspline_iterfit(wave, inv_sens, yfit=yfit2, everyn=everyn, /groupbadpix, maxrej=5, $
                             lower=2., upper=2., invvar=inv_sens_ivar*next_mask, /silent, outmask=outmask2)

      bset = bspline_iterfit(wave, sens, yfit=yfit, everyn=everyn, maxiter=1, $
                             lower=2., upper=2., invvar=ivar, /silent)

      
      svset = { $
                fullbkpt: bset.fullbkpt, $
                bkmask: bset.bkmask, $
                nord: bset.nord, $
                coeff: bset.coeff, $
                icoeff: bset.icoeff, $
                ordr: objstr[qq].order $
              }

      if keyword_set( CHKFIT ) then x_splot, wave, sens, ytwo=yfit, /block
      
      ;; Save
      if qq EQ 0 then mwrfits, svset, outfil, /create else $
        mwrfits, svset, outfil
  endfor

stop

  print, 'mike_calibstd: All Done!'
  return
end
