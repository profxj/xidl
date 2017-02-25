;+ 
; NAME:
; wfccd_extobjbox
;    Version 1.0
;
; PURPOSE:
;
; CALLING SEQUENCE:
;   
;   wfccd_extobjbox, wfccd, WFARC=
;
; INPUTS:
;   wfstrct     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   wfarc      -  WFCCD arc structure (fits file)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_extobjbox, wfstrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_extobjbox, wfccd, mask_id, exp_id, SLIT_ID=slit_id, DEBUG=debug

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_extobjbox, wfccd, mask_id, [exp_id] SLIT_ID=, /DEBUG, '
    print, '          [v1.0]'
    return
  endif 

;  Optional Keywords

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'wfccd_extobjbox: Loading up the files...'

; Open Slit file
  wfslit = xmrdfits(wfccd[exp].slit_fil, 1, STRUCTYP='mslitstrct', /silent) 

; Open Obj file
  wfobj = xmrdfits(wfccd[exp].obj_fil, 1, STRUCTYP='specobjstrct', /silent)

; Open Flux, Wave, VAR, RMS
  var = xmrdfits(wfccd[exp].img_final, 1, /silent)
  wave = xmrdfits(wfccd[exp].img_final, 2, /silent)
  fx = xmrdfits(wfccd[exp].img_final, 3, /silent) ; Sky sub
  skyrms = xmrdfits(wfccd[exp].img_final, 4, /silent)

  sz = size(fx, /dimensions)

;  Mask
  msk = bytarr(sz[0],sz[1])

; Guess
  xyguess = fltarr(2)

;;;;;;;;;;;;
;  Loop on Slits

  ; DEBUG
  if keyword_set(SLIT_ID) then gdslit = where(wfslit.id EQ slit_id, nslit) $
  else gdslit = where(wfslit.flg_anly NE 0, nslit)

  for qq=0L,nslit-1 do begin
      if not keyword_set( SILENT ) then $
        print, 'wfccd_extobjbox: Extracting slit '+strtrim(qq,2)+$
        ' of '+strtrim(nslit-1,2)+'  Obj: ', strtrim(wfslit[gdslit[qq]].id,2)
      ; Reset mask
      msk = msk*0B
      ;; OBJ
      objslit = where(wfobj.slit_id EQ wfslit[gdslit[qq]].id $
                     AND wfobj.flg_anly NE 0, nobj)
      if nobj EQ 0 then begin
          print, 'wfccd_extobjbox:  No good obj in this slit! ', $
            wfslit[gdslit[qq]].id
          continue
      endif
      ; Grab the frame
      ymn = long(min(wfslit[gdslit[qq]].yedg_orig[0:sz[0]-1,0]))
      ymx = (long(max(wfslit[gdslit[qq]].yedg_orig[0:sz[0]-1,1]))+1) $
        < (sz[1]-1)
      subfx = fx[*,ymn:ymx]
      subwv = wave[*,ymn:ymx]
      subvar = var[*,ymn:ymx]

      ; Add SKY RMS to VAR 
      skyvar = skyrms[*,gdslit[qq]]^2
      szsub = size(subvar, /dimensions)
      imgsky = skyvar # replicate(1.,szsub[1])
      ; DO NOT ADD TO REJ PIX
      a = where(subvar GT 0.)
      subvar[a] = subvar[a] + imgsky[a]

      ; Make the Mask
      for i=0L, sz[0]-1 do begin
          yedg = round(wfslit[gdslit[qq]].yedg_orig[i,*])
          msk[i,yedg[0]:yedg[1]] = 1
      endfor
      submsk = msk[*,ymn:ymx]

      ; Zero out other slits
      badmsk = where(submsk EQ 0, nbad)
      if nbad NE 0 then begin
          subfx[badmsk] = 0.
          subwv[badmsk] = 0.
          subvar[badmsk] = 0.
      endif
      
      ;; Name
      allnm = wfobj[objslit].obj_id

      ;; Loop on Objects
      for jj=0L,nobj-1 do begin
          ; Find Obj
          nm = strtrim(allnm[jj],2)
          obj = (where(wfobj[objslit].obj_id EQ nm))[0]
          ; Guess
          xyguess[0] = wfobj[objslit[obj]].xcen
          xyguess[1] = wfobj[objslit[obj]].ycen - ymn

          ;; Boxcar
          x_extobjbox, subfx, subwv, xyguess, spec, MSK=submsk, $
            WVMNX=[3200.,11000], PIX=pix, FRAC=frac, VAR=subvar, $
            CRVAL1=alog10(3400.d), CDELT=1.448662d-4,$ ; 100 km/s pixels
            NPIX=3500L, COLLMNX=[3700., 8000.], DEBUG=debug, /REDBLUE, $
            /REJ_CR, /REBINC, BKAPER=[3.,3.], SIG_COLL=0.5, /SILENT

          ; Write to structure
          if spec.npix NE 0 then begin
              wfobj[objslit[obj]].npix = spec.npix
              wfobj[objslit[obj]].wave[0:spec.npix-1] = spec.wv
              wfobj[objslit[obj]].fx[0:spec.npix-1] = spec.fx
              wfobj[objslit[obj]].var[0:spec.npix-1] = spec.var
              ; Trace
              ntrc = n_elements(spec.trc)
              wfobj[objslit[obj]].trace[0:ntrc-1] = $
                spec.trc+replicate(ymn, ntrc) ; Offset!
              ; Aper
              wfobj[objslit[obj]].aper = spec.aper
          endif else wfobj[objslit[obj]].flg_anly = 0
      endfor
  endfor

  ; Finish up
  if not keyword_set( SILENT ) then print, 'wfccd_extobjbox: All done!'
  mwrfits, wfobj, wfccd[exp].obj_fil, /create

  return
end
          
  

