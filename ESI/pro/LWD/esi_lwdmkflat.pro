;+ 
; NAME:
; esi_lwdmkflat   
;     Version 1.0
;
; PURPOSE:
;    Create a bias subtracted, median FLAT
;      Defaults to Dome Flat
;
; CALLING SEQUENCE:
;   
;  esi_lwdmkflat, esi, /INTERNAL
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   esi   -  Creates a combined ZRO frame for data reduction
;
; OPTIONAL KEYWORDS:
;  FLAT_FIL= -- 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdmkflat, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdmkflat, esi, slit, CHK=chk, FLAT_FIL=flat_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdmkflat, esi, slit, /CHK, FLAT_FIL= [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
  flt = 'DFLT'

  ;; Slit
  c_s = esi_slitnm( slit )
; Check for Arc IMG
  arcfil = 'Arcs/AIMG_LWD'+c_s+'.fits'
  a = findfile(arcfil, count=na)
  if na EQ 0 then begin
      print, 'esi_mklwdflat: No Arc Image!  Make this first...'
      return
  endif
  
  ;;  Grab the files
  flat = where(esi.mode EQ 1 AND esi.flg_anly NE 0 AND $
               esi.type EQ flt AND esi.slit EQ slit, nflat)
  if nflat EQ 0 then begin
      print, 'esi_lwdmkflat:  No Flats to Process!!'
      stop
      return
  endif

  if not keyword_set( flat_fil ) then begin
      ;; BIAS subtract
      if not keyword_set( REDOOV ) then $
        bias = where(esi[flat].flg_ov EQ 0, nbias) $
      else begin
          nbias = n_elements(flat)
          bias = lindgen(nbias)
      endelse
      if nbias NE 0 then esi_subbias, esi, flat[bias]
      
      ;; Median Combine
      xcombine, 'OV/ov_'+esi[flat].img_root, img_flat, head, $
        FCOMB=2, SCALE='MED', GAIN=esi[flat[0]].gain, RN=esi[flat[0]].readno
    
  endif else img_flat = xmrdfits(flat_fil,/silent)
  sz_flat = size(img_flat, /dimensions)
      
  ;; Read Arc Image
  img_arc = mrdfits(arcfil, /silent)
  sz_arc = size(img_arc, /dimensions)
  if sz_arc[0] NE sz_flat[0] OR sz_arc[1] NE sz_flat[1] then begin
      print, 'esi_mklwdflat: Images are wrong size!'
      stop
      return
  endif

  ;; Normalize
  wvmn = 4000.
  gd = where(img_arc GT wvmn AND img_flat GT 100., ngd)
  srt = sort(img_arc[gd])
  gd = gd[srt]

  ;; Median
  nstp = 500L
  med_wv = fltarr(nstp)
  med_fx = fltarr(nstp)

  nimg = ngd / nstp
  for qq=0L,nstp-1 do begin
      med_wv[qq] = median(img_arc[gd[qq*nimg:(qq+1)*nimg-1]])
      med_fx[qq] = median(img_flat[gd[qq*nimg:(qq+1)*nimg-1]])
  endfor

  print, 'esi_mklwdflat: Fitting...'  ; Might consider a straight spline
  bset = bspline_iterfit(med_wv, med_fx, everyn=3, maxiter=5L,$
                        upper=5, lower=5, nordr=3, yfit=yfit)

  ;; Calculate everywhere
  fit_flat = img_flat*0.
  a = where(img_arc GT wvmn)
  fit_flat[a] = bspline_valu(img_arc[a], bset)

  ;; Divide
  nrm_flat = img_flat*0.
  nrm_flat[a] = img_flat[a] / fit_flat[a]

  if keyword_set( CHK ) then xatv, nrm_flat, /block

; Output
  outfil = 'Flats/FlatLWD_'+c_s+'.fits'
  mwrfits, nrm_flat, outfil, head, /create, /silent
  spawn, 'gzip -f '+outfil
  print, 'esi_lwdmkflat: Flat created ', outfil

; Del OV
  if not keyword_set( SVOV ) then esi_delov, esi, flat

; Set Flat_fil in strcture
  obj = where(esi.mode EQ 1 and esi.flg_anly NE 0 AND $
              esi.type EQ 'OBJ' AND esi.slit EQ slit, nobj)
  if nobj NE 0 then esi[obj].flat_fil = outfil

  print, 'esi_lwdmkflat: All done!'
  return
end
