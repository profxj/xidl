;+ 
; NAME:
; ism_lowkin
;  V1.1
;
; PURPOSE:
;    Calculate kinematic characteristics for low-ion transitions
;    following the prescriptions described in Prochaska & Wolfe 1997.
;    The values are written to the ISM or DLA structure.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   dla -- DLA structure
;   nn  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  KBIN = Binning in km/s  [default: 22]
;  PER  = Percentage to toss out [default: 5% per side]
;  /CHK  --  Plot to the screen
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   ism_allabd, ism
;
;
; PROCEDURES CALLED:
;  upd_fd
;
; REVISION HISTORY:
;   14-Apr-2006 Written by JXP
;   14-Feb-2012 Added a fix for bad pixels, the value is taken as the
;               average value of the nearest correct pixels, only works 
;               if the bad pixel range is smaller than 5 pixels (MN)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro ism_lowkin, dla, nn, KBIN=kbin, PER=per, CHK=chk, SPEC_ROOT=spec_root 

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'ism_lowkin, dla, nn  [v1.0]'
    return
 endif 

  if not keyword_set(SPEC_ROOT) then spec_root = ''

  if not keyword_set(KBIN) then kbin = 22.  ; km/s
  if not keyword_set(PER) then per = 0.05  ; Percentage of total tau

  if dla[nn].flglw LE 0 then return
  c = x_constants()

  ;; Read in data
  if dla[nn].flglw GT 1 then flg = dla[nn].flglw else flg = 0
  dat = x_readspec(SPEC_ROOT+dla[nn].lwfil, /struct, /autofsig, INFLG=flg)

  
  ;; Create the velo array
  x_pixminmax, dat.wv, dla[nn].lwwav, dla[nn].zabs, dla[nn].lwvmn,  $
               dla[nn].lwvmx,  PIXMIN=pmn, PIXMAX=pmx, VELO=velo


  ;; Pad for the smoothing
  pix = pmn + lindgen(pmx-pmn+1)
  npix = n_elements(pix)

  dv = abs(median( (dat.wv[pix]-shift(dat.wv[pix],1))/dat.wv[pix] $
                   * c.c /1e5) )

  ;; Fix for small set of bad pixels
  badzero=where(dat.fx[pix] eq 0 and dat.sig[pix] LE 0)
  if max(badzero)-min(badzero) ge 5 then begin
     print, 'ism_lowkin: too many or too large sections of bad data'
     stop
  endif
  dat.fx(pix(badzero))=mean([dat.fx(pix(min(badzero)-1)),dat.fx(pix(max(badzero)+1))])
  dat.sig(pix(badzero))=mean([dat.sig(pix(min(badzero)-1)),dat.sig(pix(max(badzero)+1))])

;  nbin = 7L
  nbin = round(kbin/dv)
  if keyword_set(CHK) then print, 'ism_lowkin: Binning with ', nbin, ' pixels'

  ;; Pad for smoothing
  ppix = (pmn-nbin) + lindgen(npix+2*nbin)
  tau = fltarr(npix+2*nbin)

  ;; Create the optical depth array
  gd = where(dat.fx[ppix] GT dat.sig[ppix]/2. AND $
             dat.sig[ppix] GT 0.,ngd, $
             complement=sat,ncomplement=nsat)
  if ngd EQ 0 then begin
      print, 'ism_lowkin: Profile too saturated'
      stop
  endif
  tau[gd] = alog(1.d/dat.fx[ppix[gd]])
  ;; Set saturated pix to the median S/N ratio
  if nsat NE 0 then tau[sat] = alog(median(dat.sig[ppix]/dat.fx[ppix]))

  ;; Smooth
  stau = smooth(tau,nbin)
  ftau = stau[nbin:nbin+npix-1] 


  ;; Velocity width
  tottau = total( ftau )
  
  ;; New edges
  cumtau = total(ftau, /cumulative) / tottau
  lft = (where(cumtau GT per))[0]
  rgt = (where(cumtau GT (1.-per)))[0] - 1   ;; To match original definition


  ;; Evaluate
  dla[nn].lwfvel = float(round(abs(velo[pix[rgt]]-velo[pix[lft]])))
  if keyword_set(CHK) then $
    print, 'ism_lowkin: Velocity interval = ', dla[nn].lwfvel

  ;;;;;
  ;; Mean/Median
  vcen = (velo[pix[rgt]]+velo[pix[lft]])/2.
  mean = dla[nn].lwfvel/2.
  mn = min(abs(cumtau-0.5), imn)

  dla[nn].lwfmm = abs( (velo[pix[imn]]-vcen)/mean )
  if keyword_set(CHK) then print, 'ism_lowkin: Mean median = ', dla[nn].lwfmm

  ;; ;;;;;;;;;;;;;;;;;;;
  ;; Tau-weighted mean
  if tag_exist(dla[nn], 'tau_mean') then begin
     dla[nn].tau_mean = $ 
        total( velo[pix] * ftau ) / total( ftau )  ;; Using the smoothed tau
  endif
  
  ;; Peak
  mx = max(ftau,imx)
;  mean = (velo[pix[rgt]]+velo[pix[lft]])/2.
  dla[nn].lwfedg = abs( (velo[pix[imx]]-vcen) / mean )
  if keyword_set(CHK) then print, 'ism_lowkin: Edge test = ', dla[nn].lwfedg

  ;; CHK
  if keyword_set(CHK) then begin
      clr = getcolor(/load)
      csz = 1.5
      plot, velo[pix], ftau, color=clr.black, backgr=clr.white, psym=10
      xyouts, velo[pix[lft]]+dv, 0.1, dla[nn].qso+', z='+ $
              string(dla[nn].zabs,format='(f5.3)'), color=clr.black, charsize=csz
      ;; Limits
      oplot, replicate(velo[pix[lft]],2), [-1e9,1e9], color=clr.green
      oplot, replicate(velo[pix[rgt]],2), [-1e9,1e9], color=clr.green
      xyouts, velo[pix[lft]]+dv, 0.05, $
              'Dv = '+string(dla[nn].lwfvel,format='(f5.1)')+' km/s',$
              color=clr.black, charsiz=csz
      ;; Mean/median
      oplot, replicate(velo[pix[imn]],2), [-1e9,1e9], color=clr.blue, linest=1
      ;; Edge
      oplot, replicate(velo[pix[imx]],2), [-1e9,1e9], color=clr.red, linest=2
      ;; Tau-weighted mean
      if tag_exist(dla[nn], 'tau_mean') then $
         oplot, replicate(dla[nn].tau_mean, 2), [-1e9,1e9], color=clr.brown, linest=2
      wait, 3
  endif
  
  return
end
