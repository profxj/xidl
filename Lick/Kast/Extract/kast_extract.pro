;+ 
; NAME:
; kast_extract   
;   Version 1.1
;
; PURPOSE:
;    Plots any array interactively
;
; CALLING SEQUENCE:
;   
;   spec = x_apall(ydat, [head])
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   wave       - wavelength array
;   DISPLAY    - Display the sky subtracted image with xatv
;   OVR        - String array for ov region:  '[2050:2100, *]'
;   ERROR      - Variance array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   spec = kast_extract('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_extract, kast, setup, side, obj_id, exp, $
                  APER=aper, SKYNORD=skynord, STD=std, SILENT=silent, $
                  YMODEL=ymodel, BOXCAR=boxcar


;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_extract, kast, setup, side, obj_id, [exp]  [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( SKYNORD ) then skynord = 2
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0
  if not keyword_set( AOFF ) then aoff = 5.

; Find objects
  if not keyword_set( STD ) then begin
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.obj_id EQ obj_id AND kast.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'kast_extract: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'kast_extract, kast, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
      indx = indx[obj_id]
      nindx = 1L
      ;; Set boxcar
      BOXCAR=1
      aoff = 15.
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Loop
  for q=0L,n_elements(exp)-1 do begin

      ;; Open objfil
      objfil = kast[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'kast_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+kast[indx[exp[q]]].img_root
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'kast_extract: No Final file ', imgfil
          stop
      endif
      img = xmrdfits(imgfil, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz = size(img, /dimensions)

      ;; Set Sigma
      if not keyword_set( APER ) then begin
          sigma = fltarr(nobj)
          if not keyword_set( SILENT) then $
            print, 'kast_extract: Define the Aperture'
          for kk=0L,nobj-1 do begin
              cline = objstr[kk].xcen
              spec = djs_median(img[cline-15L:cline+15L, *],1)
              aper = x_setaper(spec, objstr[kk].ycen, 0.05, RADIUS=20L)
              ;; STD
              if keyword_set(STD) then aper = aper + 3.
              objstr[kk].aper = aper
              sigma[kk] = total(aper)/4.
          endfor
      endif

      if not keyword_set( BOXCAR ) then begin
          ;; TRIM
          img = img[*,20:160]
          var = var[*,20:160]

          ;; IVAR
          ivar = 1. / transpose(var)

          ;; Setup the trace
          xcen = objstr.trace[0:sz[0]-1] - 20.  ; Trim offset

          ;; Kludging!
          if nobj EQ 1 then xcen = [[xcen], [xcen*0.+9999.]]
          
          ;; Optimal
          extract_image, transpose(img), ivar, xcen, $
            sigma, flux, finv, ymodel=ymodel, nPoly=skynord, wfixed=[1,1,1]

          if keyword_set( CHK ) then xatv, ymodel, /block

          ;; Fill up objstr
          for kk=0L,nobj-1 do begin
              objstr[kk].fx[0:sz[0]-1] = flux[*,kk]
              objstr[kk].var[0:sz[0]-1] = 1./finv[*,kk]
              objstr[kk].trace[0:sz[0]-1] = xcen[*,kk]
              objstr[kk].npix = sz[0]-1
          endfor

      endif else begin
          ;; SKY SUBTRACTION
          if not keyword_set( SKYFSTR ) then begin
              skyfstr = { fitstrct }
              if not keyword_set( SKYFUNC ) then skyfstr.func = 'POLY'
              if not keyword_set( SKYNORD ) then skyfstr.nord = 2
              if not keyword_set( SKYLOW ) then skyfstr.lsig = 2.
              if not keyword_set( SKYHIGH ) then skyfstr.hsig = 2.
              skyfstr.niter = 2
              skyfstr.maxrej = 5
          endif
          ;; Key on primary science object only
          skyreg = [ [-30.-objstr[0].aper[0]-aoff, objstr[0].aper[1]+aoff], $
                     [-objstr[0].aper[0]-aoff, objstr[0].aper[1]+aoff+30.]] 
          skyreg[*] = skyreg[*] + objstr[0].ycen
          fimg = x_skysub(img, sky, REG=skyreg, CREG=objstr[0].xcen, $
                          TRACE=objstr[0].trace[0:sz[0]-1], $
                          NOREJ=norej, SKYFSTR=skyfstr, SKYRMS=skyrms) 
          ;; BOXCAR
          flux = x_extract(fimg, [objstr[0].ycen-objstr[0].aper[0], $
                                  objstr[0].ycen+objstr[0].aper[1]], $
                           objstr[0].trace[0:sz[0]-1], $
                           fvar, sky, $
                           CAPER=objstr[0].xcen, $
                           SKYRMS=skyrms, RN=kast[indx[exp[q]]].readno, $
                           GAIN=1.)
          ;; Fill it up
          objstr[0].fx[0:sz[0]-1] = flux
          objstr[0].var[0:sz[0]-1] = fvar
          objstr[0].npix = sz[0]-1
      endelse
          
      ;; Write
      print, 'kast_extract: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil
                     
  endfor


  print, 'kast_extract: All Done!'
  return
end
