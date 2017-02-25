;+ 
; NAME:
; kast_extract   
;   Version 1.1
;
; PURPOSE:
;    Extract a 1D spectrum from a Kast 2D image.
;
; CALLING SEQUENCE:
;  kast_extract, kast, setup, side, obj_id, [exp], $
;                 APER=, SKYNORD=, /STD, YMODEL=, /BOXCAR
;
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;   side  --  Specific camera [blue (1) vs. red (2)]
; obj_id  --  Object value
;  [exp]  --  Exposure indices
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /BOXCAR -- Extract with a boxcar [default: optimal]
;  APER= -- Aperture (for boxcar primarily)
;
; OPTIONAL OUTPUTS:
;  YMODEL=  -- Image fit created in optimal extraction
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_extract, kast, setup, side, obj_id, exp, $
                  APER=aper, SKYNORD=skynord, STD=std, $
                  YMODEL=ymodel, BOXCAR=boxcar, SKYBND=skybnd, CHK=chk
                  


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
  if not keyword_set( SKYBND ) then skybnd = 30.

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
      if strlen(strtrim(objfil,2)) eq 0 or x_chkfil(objfil+'*') EQ 0 then begin
          print, 'kast_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+kast[indx[exp[q]]].img_root
      if strlen(strtrim(imgfil,2)) eq 0 or x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'kast_extract: No Final file ', imgfil
          stop
      endif
      img = xmrdfits(imgfil, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz = size(img, /dimensions)

      ;; Set Sigma
      if not keyword_set( APER ) then begin
          sigma = fltarr(nobj)
          for kk=0L,nobj-1 do begin
              cline = objstr[kk].xcen
              spec = djs_median(img[cline-15L:cline+15L, *],1)
              aper = x_setaper(spec, objstr[kk].ycen, 0.05, RADIUS=20L)
              ;stop
              ;; STD
              if keyword_set(STD) then aper = aper + 3.
              objstr[kk].aper = aper
              sigma[kk] = total(aper)/4.
          endfor
      endif
      
      ;; Skyregion
      msk = bytarr(2*round(SKYBND) + round(total(objstr[0].aper)) + $
                   2*aoff)
      nmsk = n_elements(msk)
      xmsk = objstr[0].ycen-objstr[0].aper[0]-SKYBND + findgen(nmsk)
      ;; Sciobj
      scim = where(xmsk GT (objstr[0].ycen-objstr[0].aper[0]-aoff) AND $
                   xmsk LT (objstr[0].ycen+objstr[0].aper[1]+aoff) )
      msk[scim] = 1B

      ;; Serendip
      for ii=1L,n_elements(objstr)-1 do begin
          serm = where(xmsk GT (objstr[ii].ycen-objstr[ii].aper[0]) AND $
                       xmsk LT (objstr[ii].ycen+objstr[ii].aper[1]),nserm)
          if nserm NE 0 then msk[serm] = 1B
      endfor

      ;; Make the Regions
      subms = xmsk[where(msk EQ 0B,nsub)]
      shfsub = subms - shift(subms,1)
      bound = where(shfsub GT 1., nbound)
      lhs = [subms[0]]
      rhs = [subms[nsub-1]]
      for ii=0L,nbound-1 do begin
          lhs = [lhs, subms[bound[ii]]]
          rhs = [subms[bound[ii]-1],rhs]
      endfor
      skyreg = fltarr(n_elements(lhs),2)
      for ii=0L,n_elements(lhs)-1 do begin
          skyreg[ii,0] = lhs[ii]
          skyreg[ii,1] = rhs[ii]
      endfor

      ;; EXTRACT
      if not keyword_set( BOXCAR ) then begin

          ;; IVAR
          ivar = 1. / transpose(var)
          sivar = ivar

          ;; Mask Out Serendip objects (just use trace of Sciobj)
          ymsk = findgen(sz[1])
          for ii=1L,n_elements(objstr)-1 do begin
              yoff = objstr[ii].ycen - objstr[0].ycen
              ;; Mask
              msku = (objstr[0].trace + yoff + objstr[ii].aper[1]) < $
                float(sz[1]-1)
              mskd = (objstr[0].trace + yoff - objstr[ii].aper[0]) > 0.
              for jj=0L,sz[0]-1 do begin
                  bad = where(ymsk GE mskd[jj] AND ymsk LE msku[jj],nbd)
                  if nbd NE 0 then sivar[bad,jj] = -1.
              endfor
          endfor

          ;; TRIM
          img = img[*,20:160]
          var = var[*,20:160]
          ivar = ivar[20:160,*]
          sivar = sivar[20:160,*]

          ;; Setup the trace
          xcen = objstr[0].trace[0:sz[0]-1] - 20.  ; Trim offset
;          xcen = objstr.trace[0:sz[0]-1] - 20.  ; Trim offset

          ;; Kludging!
          xcen = [[xcen], [xcen*0.+9999.]]

;          if flg_nobj EQ 1 then xcen = [[xcen], [xcen*0.+9999.]] $
;          else begin
;              ;;added 10/24/04 by SLM: sort xcen into ascending order
;              resort=sort(xcen[0,*])
;              sxcen=xcen
;;              ssigma=sigma
;              nloop=n_elements(xcen[0,*])
;              
;              for i=0,nloop-1 do begin
;                  sxcen[*,i]=xcen[*,resort[i]]
;                  ssigma[i]=sigma[resort[i]]
;              endfor
;          
;              xcen=sxcen      
;              sigma=ssigma 
;              ;;end SLM section
;          endelse


          ;; Optimal
          szi = size(ivar,/dimensions)
          mask = bytarr(szi[0],szi[1])
          mask[*] = 1B
          bad = where( sivar LT 0., nbad)
          if nbad NE 0 then mask[bad] = 0B
          extract_image, transpose(img), ivar, xcen, mask=mask, $
            sigma, flux, finv, ymodel=ymodel, nPoly=skynord, wfixed=[1,1,1]

          if keyword_set( CHK ) then xatv, transpose(img)-ymodel, /block

          ;; Fill up objstr
          for kk=0L,nobj-1 do begin
;          for kk=0L,1L do begin
              objstr[kk].fx[0:sz[0]-1] = flux[*,kk]
              objstr[kk].var[0:sz[0]-1] = 1./finv[*,kk]
              objstr[kk].trace[0:sz[0]-1] = xcen[*,kk]
              objstr[kk].npix = sz[0]-1
          endfor

      endif else begin  ;; BOXCAR
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

          ;; Subtract
          fimg = x_skysub(img, sky, REG=skyreg, CREG=objstr[0].xcen, $
                          TRACE=objstr[0].trace[0:sz[0]-1], $
                          NOREJ=norej, SKYFSTR=skyfstr, SKYRMS=skyrms) 
          if keyword_set(CHK) then xatv, fimg, /blo

          ;; BOXCAR
          flux = x_extract(fimg, [objstr[0].ycen-aper[0], $
                                  objstr[0].ycen+aper[1]], $
                           objstr[0].trace[0:sz[0]-1], $
                           fvar, sky, $
                           CAPER=objstr[0].xcen, $
                           SKYRMS=skyrms, RN=kast[indx[exp[q]]].readno, $
                           GAIN=1., VAR=var)
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
