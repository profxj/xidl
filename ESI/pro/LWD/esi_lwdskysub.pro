;+ 
; NAME:
; esi_lwdskysub   
;     Version 1.0
;
; PURPOSE:
;    Sky Subtract image
;
; CALLING SEQUENCE:
;   
;  esi_lwdskysub, esi, obj_id
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_lwdskysub, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_lwdskysub, esi, obj_id, expsr, SKYREG=skyreg, REFWV=refwv, CHK=chk, $
                   LOWORD=loword, STD=std

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_lwdskysub, esi, obj_id, [exspr]  [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( SKYREG ) then skyreg = 300L
  if not keyword_set(REFWV) then refwv = [6550., 6650.]
  if not keyword_set(LOWORD) then loword = 2L

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 1 AND $
                   esi.obj_id EQ obj_id AND esi.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_lwdskysub: No images to find obj for!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

  skylin = fltarr(50,3)
  skylin[0,*] = [5563, 5588, 3]
  skylin[1,*] = [7700, 8000, 1]
  skylin[2,*] = [8250, 9100, 1]
  skylin[3,*] = [9350, 10100, 1]

  if keyword_set( STD ) then skylin[*,2] = skylin[*,2] - 1

;  Loop

  for q=0L,n_elements(exp)-1 do begin
      print, 'esi_lwdskysub: Reading files...'
      ;; Open Obj file
      objfil = esi[indx[q]].obj_fil
      a = findfile(objfil+'*', count=na)
      if na EQ 0 then begin
          print, 'esi_lwdskysub: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
      nobj = n_elements(objstr)

      ;; IMG+VAR Fil 
      imgfil = objstr[0].spec2d_fil
      img = xmrdfits(imgfil, 0, head, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      bdvar = where(var LE 0.)
      var[bdvar] = -1.
      ivar = 1./var
      sz_img = size(img, /dimensions)

      ;; Read ARC
      img_arc = xmrdfits(strtrim(esi[indx[exp[q]]].arc_fil,2), /silent) 


      ;; Include around sci obj following the trace
      msk = bytarr(sz_img[0], sz_img[1])
      sci = where(objstr.obj_id EQ 'a')
      for ii=0L,sz_img[0]-1 do begin
          jmn = 0 > (objstr[sci].trace[ii] - skyreg)
          jmx = (sz_img[1]-1) < (objstr[sci].trace[ii] + skyreg)
          msk[ii,jmn:jmx] = 1
      endfor

      ;; Mask out Obj
      print, 'esi_lwdskysub: Masking Obj'
      for jj=0L,nobj-1 do begin
          for ii=0L,sz_img[0]-1 do begin
              jmn = 0 > round(objstr[jj].trace[ii] - objstr[jj].aper[0]-5.)
              jmx = (sz_img[1]-1) < round(objstr[jj].trace[ii] $
                                     + objstr[jj].aper[1]+5.)
              msk[ii,jmn:jmx] = 2
          endfor
      endfor
      
      ;; Mask only good wave region
      bdwv = where(img_arc LT 4000. OR img_arc GT 10000.)
      msk[bdwv] = 0

      ;; Sky Shape
      msk_mul = bytarr(sz_img[0], sz_img[1])
      sva = where(img_arc GT refwv[0] AND img_arc LT refwv[1] AND msk EQ 1)
      msk_mul[sva] = 1
      smsh = fltarr(sz_img[1])
      print, 'esi_lwdfndobj: Fitting Sky Shape'
      for jj=0L,sz_img[1]-1 do begin
          b = where(msk_mul[*,jj] EQ 1, nb)
          if nb NE 0 then smsh[jj] = median(img[b,jj])
      endfor
      msk_smsh = bytarr(sz_img[1])
      msk_smsh[where(smsh > 0.)] = 1
      mult_fitstr = {fitstrct}
      mult_fitstr.nord = loword
      mult_fitstr.func = 'POLY'
      mult_fitstr.hsig = 2.5
      mult_fitstr.lsig = 2.5
      mult_fitstr.niter = 2
      mult_fit = x_fitrej(findgen(sz_img[1]), smsh, MSK=msk_smsh, $
                          FITSTR=mult_fitstr)
      ;; Normalize
      nrm_mult = median(mult_fit[where(msk_smsh EQ 1)])
      if keyword_set( CHK ) then x_splot, smsh, YTWO=mult_fit, /block
      mult_fit = mult_fit/nrm_mult


      ;; Convert to 1D
      print, 'esi_lwdskysub: Grabbing sky pixels'
      skypix = where(msk EQ 1, nsky)
      srt = sort(img_arc[skypix])
      sky_wv = img_arc[skypix[srt]]

      ;; Divide Out Shape of sky
      print, 'esi_lwdskysub: Dividing out sky shape'
      ;; Find trace intersection
      msk_trc = bytarr(sz_img[0])
      for ii=0L,sz_img[0]-1 do begin
          jj = round(objstr[sci].trace[ii])
          if img_arc[ii,jj] GT refwv[0] AND img_arc[ii,jj] LT refwv[1] then $
            msk_trc[ii] = 1
      endfor
      shp_col = round(median( where(msk_trc EQ 1)))

      ;; Create Shape image
      img_shape = fltarr(sz_img[0], sz_img[1])
      for ii=0L,sz_img[0]-1 do begin
          img_shape[ii,*] = shift(mult_fit, round(objstr[sci].trace[ii]- $
                                                  objstr[sci].trace[shp_col]))
      endfor

      ;; Normalize the Sky and grab the pix
      nrm_sky = img/img_shape
      sky_fx = nrm_sky[skypix[srt]]
      sky_ivar = ivar[skypix[srt]]

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; SETUP BKPTS
      print, 'esi_lwdskysub: Setting Break points'
      wvoned = fltarr(sz_img[0])
      for i=0L,sz_img[0]-1 do wvoned[i] = img_arc[i,round(objstr[0].trace[i])]
      bkpts = wvoned
      ;; Add half pixel
      hlf = (wvoned+shift(wvoned,1))/2.
      bkpts = [bkpts,hlf[1:sz_img[0]-2]]
      ;; Good ones
      gdbk = where(bkpts GT 3990 AND bkpts LT 10050)
      bkpts = bkpts[gdbk]
      ;; Improve sky lines
      sz_slin = size(skylin, /dimensions)
      for k=0L,sz_slin[0]-1 do begin
          case round(skylin[k,2]) of
              -1: 
              0: 
              1: begin ;; Every 5 Ang
                  bkpts = [bkpts, 5.*findgen(round(skylin[k,1]-skylin[k,0])/5)+ $
                           skylin[k,0]]
              end
              2: begin ;; Every ang
                  bkpts = [bkpts, findgen(round(skylin[k,1]-skylin[k,0])+1)+ $
                           skylin[k,0]]
              end
              3: begin ;; Every 1/2 ang
                  bkpts = [bkpts, 0.5*findgen(round(skylin[k,1]-skylin[k,0])*2+2)+ $
                           skylin[k,0]]
              end
              else: stop
          endcase
      endfor
      ;; Sort
      srt = sort(bkpts)
      bkpts = bkpts[srt]
      ;;;;;;;;;;;;;;
      ;; FIT with Bspline
      print, 'esi_lwdskysub: Fitting ', nsky, ' sky pixels with Bspline'
      bset = bspline_iterfit(sky_wv, sky_fx, bkpt=bkpts, nord=3, $
                             upper=2.5, lower=3.0,INVVAR=sky_ivar, maxiter=5L)

      ;;;;;;;;;;
      if keyword_set(CHK) then begin
          nfit = 5000L
          x0 = sky_wv[0]
          xN = sky_wv[n_elements(sky_wv)-1]
          xfit = fltarr(nfit)
          for i=0L,nfit-1 do xfit[i] = x0 + $
            float(i)*(xN-x0)/float(nfit)
          yfit = bspline_valu(xfit, bset)
          x_splot, sky_wv, sky_fx, PSYM1=3, /block, XTWO=xfit, YTWO=yfit
      endif

      ;; Create sky image
      print, 'esi_lwdskysub: Creating a sky image'
      img_sky = img*0.
      gdreg = where(msk NE 0 AND img_arc GT 0.)
      img_sky[gdreg] = bspline_valu(img_arc[gdreg], bset)

      ;; Subtract from image
      print, 'esi_lwdskysub: Subtracting sky'
      tmp = img_sky*img_shape
      img_new = img - tmp
      badpix = where(msk EQ 0 OR img_arc LE 0.)
      img_new[badpix] = 0.

      ;; Display
      if keyword_set( CHK ) then begin
          xatv, img_new, min=-20., max=200., /block
      endif
      
      print, 'esi_lwdskysub: Writing sky image to ', imgfil
      ;; Ouptut New Image
      mwrfits, img, imgfil, head, /create, /silent
      mwrfits, var, imgfil, /silent
      mwrfits, img_new, imgfil, /silent
      spawn, 'gzip -f '+imgfil
  endfor
  
;  DONE
  print, 'esi_lwdskysub: All done! '
  return
end
