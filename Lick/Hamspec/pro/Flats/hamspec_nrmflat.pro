;+ 
; NAME:
; hamspec_nrmflat
;     Version 1.1
;
; PURPOSE:
;  To normalize the flat (remove the blaze) to create an image useful
;  for pixel-to-pixel corrections.  Ideally, I recommend you use one
;  of the archived flat frames as the flats created here have the
;  following problems:
;  1.  They fail at the slit edges
;  2.  They contain scattered light which cannot be removed
;  3.  Large chip defects cause more trouble
;
; CALLING SEQUENCE:
;   
;  hamspec_nrmflat, hamspec, setup, [side], /CHK, /CLOBBER, /DEBUG, SCATT=, XSEP=
;
; INPUTS:
;   hamspec     -  MIKE structure
;   setup    -  Setup identifier 
;   [chip]   -  Blue (1), Green (2), Red (3), or multiple (Default:
;              [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  A normalized flat image
;
; OPTIONAL KEYWORDS:
;  XSEP=       --  [Default: -1]
;
; OPTIONAL OUTPUTS:
;  SCATT_IMG=  -- Scattered light image
;
; COMMENTS:
;
; EXAMPLES:
;   hamspec_nrmflat, hamspec, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  
pro hamspec_nrmflat, hamspec, setup, CHK=chk, CLOBBER=clobber, DEBUG=debug, $
                  SCATT=scatt_img, XSEP=xsep

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hamspec_nrmflat, hamspec, setup, [chip], /CLOBBER,'
      print, '           /CHK, /DEBUG, XSEP= [v1.1]'
      return
  endif 

  ;; Optional Keywords
  if not keyword_set( NRMMED ) then NRMMED = 5L
  if not keyword_set( XSEP ) then xsep = 1L

;  x_psclose
;  !p.multi=[0,1,1]

; Loop on chip

  idx = where(hamspec.setup EQ setup and hamspec.flg_anly NE 0 and $
              hamspec.type EQ 'PFLT',nidx)
  if nidx EQ 0 then stop
  rbin = hamspec[idx[0]].rowbin
  ;; Vignetting due to CCD cover
  ;upedg = round(3990./rbin) 
  nxbkpt=5L 
  nybkpt=5L

  ;; Check output file
  ;out_fil = hamspec_getfil('nqtz_fil', setup, /name, CHKFIL=chkf)
  out_fil = hamspec_getfil('npixflt_fil', setup, /name, CHKFIL=chkf)
  if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
     print, 'hamspec_nrmflat: Normalized flat exists, moving on..'
     return
  endif

  ;; Read in flat 
  ;flat_fil = hamspec_getfil('qtz_fil', setup, /name)
  flat_fil = hamspec_getfil('pixflt_fil', setup, /name)
  flat = xmrdfits(flat_fil, 0, /silent)
  flat_ivar = xmrdfits(flat_fil, 1, /silent)

  sz = size(flat, /dimensions)
  nrm_flat = fltarr(sz[0],sz[1])
  nrm_ivar = fltarr(sz[0],sz[1])
  img_nrm = fltarr(sz[0],sz[1])

  ;; Read in Order structure
  ordr_str = hamspec_getfil('ordr_str', setup)
  
  ;; Fit parameters
  slit_cen = round((ordr_str.lhedg + ordr_str.rhedg)/2.)
  med_img = fltarr(sz[1], NRMMED*2+1L) 
  med_img[*] = 1.

  ;; Subtract scattered light
  msk = x_ordermask(sz[0], sz[1], ordr_str, trim=2)
  ;scatt_img = x_fitgap(flat, flat_ivar, msk, $
  ;                     nxbkpt=nxbkpt, nybkpt=nybkpt)
  ;flat = flat - scatt_img

  ;; Loop on Orders
  nordr= n_elements(ordr_str)
  for q=1L,nordr-1 do begin
;          if q GT nordr-5 then stop
     ;; Take median along slit centers
     gd = where(slit_cen[*,q] GT NRMMED and $
                (slit_cen[*,q]+NRMMED) LT (sz[0]-1),ngd)
     frstj = gd[0]
     lstj = gd[ngd-1] ;< upedg
     for j=frstj,lstj do begin
        p0 = (slit_cen[j,q]-NRMMED)>0
        p1 = (slit_cen[j,q]+NRMMED)<(sz[0]-1)
        np = p1-p0+1
        if np GT 0 then med_img[j,0:np-1] = flat[p0:p1,j]
     endfor
     med_nrm = djs_median(med_img, 2)
     ;if q EQ 10 then stop
     ;; FIT
     bset = bspline_iterfit(frstj+findgen(lstj-frstj+1),med_nrm[frstj:lstj], $
                            yfit=fit, everyn=round(15L/float(rbin)))
;          x_splot, frstj+findgen(lstj-frstj+1),med_nrm[frstj:lstj], $
;            ytwo=fit, /blo
;          stop
     ;; LEFT and RIGHT
     if q EQ 0 then lft = (ordr_str[q].lhedg[*]-XSEP) > 0L $
     else lft = (ordr_str[q-1].rhedg[*]+XSEP) > $
                (ordr_str[q].lhedg[*]-XSEP) < (sz[0]-1)
     if q EQ (nordr-1) then rgt = lft > (ordr_str[q].rhedg[*]+XSEP) $
                                  < (sz[0]-1) $
     else rgt = lft > (ordr_str[q].rhedg[*]+xsep) < $
                (ordr_str[q+1].lhedg[*]-XSEP)
     ;; Cuts
     rgt = round(rgt) < (sz[0]-1)
     lft = round(lft) > 0L
     ;; MAP
     for j=frstj,lstj do begin
        if lft[j] LE rgt[j] then img_nrm[lft[j]:rgt[j],j] = fit[j-frstj]
     endfor
  endfor
  a = where(img_nrm NE 0)
  nrm_flat[a] = flat[a]/img_nrm[a]
  nrm_ivar[a] = flat_ivar[a] * (img_nrm[a])^2
  
  ;; CHK
  if keyword_set(CHK) then xatv, nrm_flat, min=0.9, max=1.1, /block
  ;; Write 
  print, 'hamspec_nrmflat: Writing ', out_fil
  mwrfits, nrm_flat, out_fil, /create
  mwrfits, nrm_ivar, out_fil
  mwrfits, scatt_img, out_fil
  spawn, 'gzip -f '+out_fil
  
  ;; 
  print, 'hamspec_nrmflat: All done!'

  return
end
