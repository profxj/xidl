;+ 
; NAME:
; esi_echpskysub   
;     Version 1.0
;
; PURPOSE:
;    Sky Subtract image
;
; CALLING SEQUENCE:
;   
;  esi_echpskysub, esi, obj_id
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
;   esi_echpskysub, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Oct-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echpskysub, esi, obj_id, exp, CHK=chk, STD=std, NORD=nord, ORDR=ordr, $
                   NOVAC=novac, SKLFIL=sklfil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echpskysub, esi, obj_id, [exspr], /CHK, /STD, /NOVAC, ORDR=, NORD= [v1.0]'
      return
  endif 
  
;  Optional Keywords
;  if not keyword_set(REFWV) then refwv = [6550., 6650.]
;  if not keyword_set(LOWORD) then loword = 2L
  if not keyword_set(NORD) then nord = 3L
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( ORDR ) then begin
      ordr = [0L,9L] 
      flg_ordr = 0
  endif else begin
      flg_ordr = 1
      if n_elements(ordr) NE 2 then begin
          print, 'esi_echpskysub: ORDR must be a 2 element array'
          return
      endif
  endelse

; SKYLINES
  skylin = dblarr(10,75,3)
  if not keyword_set( SKLFIL ) then begin
      skylin[4,0,*] = [5576.d, 5579, 2]
      skylin[4,1,*] = [5889.d, 5891, 2]
      skylin[4,2,*] = [5895.d, 5897, 2]
  ;;
      skylin[5,0,*] = [6299.4, 6301.4, 2]
      skylin[5,1,*] = [6306.0, 6308.0, 1]
      skylin[5,2,*] = [6362.7, 6365.0, 2]
      ;;
      skylin[6,0,*] = [6826.4, 6830.8, 1]
      skylin[6,1,*] = [6833.2, 6835.5, 1]
      skylin[6,2,*] = [6862.6, 6865.2, 1]
      skylin[6,3,*] = [6922.4, 6924.2, 1]
      skylin[6,4,*] = [6947.8, 6950.1, 1]
      skylin[6,5,*] = [6977.2, 6979.5, 1]
      skylin[6,6,*] = [7238.0, 7242.6, 1]
      skylin[6,7,*] = [7244.0, 7246.0, 2]
      skylin[6,8,*] = [7275.6, 7277.3, 2]
      skylin[6,9,*] = [7283.7, 7285.2, 1]
      
      skylin[7,0,*] = [7238.0, 7242.6, 1]
      skylin[7,1,*] = [7244.0, 7246.0, 1]
      skylin[7,2,*] = [7275.1, 7277.3, 1]
      skylin[7,3,*] = [7283.2, 7285.2, 1]
      skylin[7,4,*] = [7302.6, 7304.8, 1]
      skylin[7,5,*] = [7314.8, 7317.6, 1]
      skylin[7,6,*] = [7328.0, 7330.4, 1]
      skylin[7,7,*] = [7339.6, 7342.3, 2]
      skylin[7,8,*] = [7357.5, 7359.8, 1]
      skylin[7,9,*] = [7367.9, 7370.9, 2]
      skylin[7,10,*] = [7390.8, 7393.5, 1]
      skylin[7,11,*] = [7400.3, 7403.0, 1]
      skylin[7,12,*] = [7522.8, 7525.5, 1]
      skylin[7,13,*] = [7570.2, 7573.2, 1]
      skylin[7,14,*] = [7711.1, 7714.1, 2]
      skylin[7,15,*] = [7714.15, 7719.3, 1]
      skylin[7,16,*] = [7749.25, 7752.1, 2]
      skylin[7,17,*] = [7758.75, 7761.4, 2]
      skylin[7,18,*] = [7779.0, 7781.6, 1]
      skylin[7,19,*] = [7792.9, 7795.8, 2]
      skylin[7,20,*] = [7807.0, 7809.6, 2]
      skylin[7,21,*] = [7820.3, 7822.9, 2]
      skylin[7,22,*] = [7840.15, 7842.5, 2]
      skylin[7,23,*] = [7852.25, 7854.8, 2]
      skylin[7,24,*] = [7859.4, 7862.0, 1]
      skylin[7,25,*] = [7866.3, 7868.9, 1]
      skylin[7,26,*] = [7869.5, 7871.8, 1]
      skylin[7,27,*] = [7877.5, 7880.1, 1]
      skylin[7,28,*] = [7888.4, 7891.3, 2]
      skylin[7,29,*] = [7912.3, 7914.8, 2]
      skylin[7,30,*] = [7919.9, 7922.2, 1]
      skylin[7,31,*] = [7948.2, 7950.2, 1]
      skylin[7,32,*] = [7963.5, 7965.8, 2]
      skylin[7,33,*] = [7992.0, 7994.3, 2]
      skylin[7,34,*] = [8013.0, 8015.3, 1]
      skylin[7,35,*] = [8024.6, 8026.9, 2]
      skylin[7,36,*] = [8060.8, 8063.6, 1]
      
      skylin[8,0,*] = [8276.9, 8284.0, 1]
      skylin[8,1,*] = [8287.2, 8289.9, 1]
      skylin[8,2,*] = [8297.5, 8300.2, 2]
      skylin[8,3,*] = [8309.4, 8311.8, 1]
      skylin[8,4,*] = [8343.2, 8346.3, 2]
      skylin[8,5,*] = [8351.4, 8354.5, 2]
      skylin[8,6,*] = [8381.1, 8383.8, 1]
      skylin[8,7,*] = [8397.8, 8400.5, 2]
      skylin[8,8,*] = [8413.8, 8416.6, 1]
      skylin[8,9,*] = [8428.8, 8431.5, 2]
      skylin[8,10,*] = [8450.9, 8453.6, 1]
      skylin[8,11,*] = [8463.8, 8466.5, 2]
      skylin[8,12,*] = [8491.9, 8495.0, 1]
      skylin[8,13,*] = [8503.1, 8506.5, 2]
      skylin[8,14,*] = [8537.3, 8540.0, 1]
      skylin[8,15,*] = [8547.1, 8550.2, 1]
      skylin[8,16,*] = [8611.9, 8640.2, 1]
      skylin[8,17,*] = [8645.6, 8647.9, 1]
      skylin[8,18,*] = [8647.91, 8652.6, 2]
      skylin[8,19,*] = [8653.0, 8657.0, 2]
      skylin[8,20,*] = [8658.0, 8662.0, 1]
      skylin[8,21,*] = [8662.7, 8667.1, 2]
      skylin[8,22,*] = [8667.7, 8672.1, 2]
      skylin[8,23,*] = [8673.1, 8677.5, 2]
      skylin[8,24,*] = [8678.8, 8683.2, 1]
      skylin[8,25,*] = [8696.5, 8709.4, 1]
      skylin[8,26,*] = [8757.5, 8762.8, 2]
      skylin[8,27,*] = [8765.2, 8769.2, 2]
      skylin[8,28,*] = [8774.8, 8779.5, 2]
      skylin[8,29,*] = [8789.5, 8792.5, 2]
      skylin[8,30,*] = [8824.1, 8825.9, 1]
      skylin[8,31,*] = [8826.0, 8828.7, 2]
      skylin[8,32,*] = [8835.0, 8837.9, 2]
      skylin[8,33,*] = [8848.3, 8851.6, 1]
      skylin[8,34,*] = [8866.2, 8868.9, 2]
      skylin[8,35,*] = [8884.5, 8887.4, 2]
      skylin[8,36,*] = [8902.0, 8904.3, 2]
      skylin[8,37,*] = [8918.2, 8920.8, 2]
      skylin[8,38,*] = [8942.2, 8944.5, 2]
      skylin[8,39,*] = [8956.7, 8959.3, 2]
      skylin[8,40,*] = [8987.3, 8989.6, 1]
      skylin[8,41,*] = [8999.7, 9003.0, 1]
      skylin[8,42,*] = [9037.1, 9039.4, 1]
      skylin[8,43,*] = [9047.9, 9051.2, 1]
      skylin[8,44,*] = [9305.0, 9312.1, 1]
      skylin[8,45,*] = [9312.11,9314.7, 2]
      skylin[8,46,*] = [9322.3, 9325.2, 1]
      skylin[8,47,*] = [9337.6, 9339.2, 1]
      skylin[8,48,*] = [8363.2, 8366.2, 1]
      
      skylin[9,0,*] = [9373.7, 9377.2, 2]
      skylin[9,1,*] = [9418.4, 9420.8, 1]
      skylin[9,2,*] = [9438.4, 9440.8, 2]
      skylin[9,3,*] = [9457.3, 9459.7, 1]
      skylin[9,4,*] = [9475.3, 9478.1, 2]
      skylin[9,5,*] = [9501.3, 9503.8, 1]
      skylin[9,6,*] = [9517.8, 9520.6, 1]
      skylin[9,7,*] = [9551.1, 9553.5, 1]
      skylin[9,8,*] = [9565.5, 9568.7, 2]
      skylin[9,9,*] = [9606.0, 9608.8, 1]
      skylin[9,10,*] = [9619.2, 9622.1, 1]
      skylin[9,11,*] = [9667.2, 9670.1, 1]
      skylin[9,12,*] = [9678.8, 9683.7, 1]
      skylin[9,13,*] = [9698.0, 9700.4, 1]
      skylin[9,14,*] = [9718.4, 9721.2, 1]
      skylin[9,15,*] = [9789.7, 9792.5, 2]
      skylin[9,16,*] = [9798.5, 9801.3, 1]
      skylin[9,17,*] = [9846.7, 9849.8, 1]
      skylin[9,18,*] = [9870.3, 9873.6, 2]
      skylin[9,19,*] = [9893.3, 9896.2, 1]
      skylin[9,20,*] = [9913.2, 9916.4, 2]
      skylin[9,21,*] = [9940.0, 9949.2, 2]
      skylin[9,22,*] = [9955.0, 9957.9, 2]
      skylin[9,23,*] = [9961.0, 9963.4, 1]
      skylin[9,24,*] = [10011.5, 10015.1, 2]
      skylin[9,25,*] = [10024.0, 10026.5, 2]
      skylin[9,26,*] = [10059.0, 10062.2, 1]
      skylin[9,27,*] = [10081.0, 10083.8, 2]
      skylin[9,28,*] = [10102.0, 10105.0, 1]
      skylin[9,29,*] = [10122.0, 10125.7, 2]
      skylin[9,30,*] = [10170.0, 10173.2, 1]
      skylin[9,31,*] = [10188.0, 10190.7, 1]
      skylin[9,32,*] = [9736.4, 9439.6, 1]
  endif else begin
      readcol, sklfil, sl_i, sl_j, sl_w1, sl_w2, sl_f, FORMAT='L,L,F,F,F'
      for i=0L,n_elements(sl_i)-1 do $
        skylin[sl_i[i],sl_j[i],*] = [sl_w1[i], sl_w2[i], sl_f[i]]
  endelse
          
      
  if not keyword_set( NOVAC ) then begin
      dumwv = skylin[*,*,0:1]
      dumwv = reform(dumwv, 1500)
      tmp = dumwv
      airtovac, dumwv
      diff = dumwv - tmp
      skylin[*,*,0:1] = reform(dumwv, [10,75,2])
  endif

  if keyword_set( STD ) then skylin[*,*,2] = skylin[*,*,2] - 1

  all_mnxwv = dblarr(10,2)
  all_mnxwv[0,*] = [3900., 4380.]
  all_mnxwv[1,*] = [4000., 4720.]
  all_mnxwv[2,*] = [4320., 5070.]
  all_mnxwv[3,*] = [4660., 5500.]
  all_mnxwv[4,*] = [5055., 5980.]
  all_mnxwv[5,*] = [5618., 6572.]
  all_mnxwv[6,*] = [6230., 7300.]
  all_mnxwv[7,*] = [7000., 8210.]
  all_mnxwv[8,*] = [8000., 9380.]
  all_mnxwv[9,*] = [9330., 10200.]

; 

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
                   esi.obj_id EQ obj_id AND $
                   strtrim(esi.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_echpskysub: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = obj_id[0]
      nindx = 1L
  endelse
      

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Read Arc
  print, 'esi_echpskysub: Reading Arc file ', $
    strtrim(esi[indx[0]].arc_fil,2)
  img_arc = xmrdfits(strtrim(esi[indx[0]].arc_fil,2), /silent) 
  sz_arc = size(img_arc, /dimensions)

; Open Slit file
  c_s = esi_slitnm(esi[indx[0]].slit)
  sedg_fil = 'Flats/SEdg_ECH'+c_s+'.fits'
  print, 'esi_echpskysub: Grabbing slit edges from: ', sedg_fil
  slit_edg = xmrdfits(sedg_fil, /silent)
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = round(slit_edg)
  ;; Hole Trace
  restore, fitfil

;  Loop

  for q=0L,n_elements(exp)-1 do begin
      ;; Open Obj file
      objfil = esi[indx[exp[q]]].obj_fil
      a = findfile(objfil+'*', count=na)
      if na EQ 0 then begin
          print, 'esi_echpskysub: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
      nobj = n_elements(objstr)/10

      ;; IMG+VAR Fil 
      imgfil = objstr[0].spec2d_fil
      skyfil = 'Sky/sky_'+esi[indx[exp[q]]].img_root
      print, 'esi_echpskysub: Reading Image files... ', imgfil
      img = xmrdfits(imgfil, 0, head, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Kludge for Sky Sub only!
;      gdvar = where(var GT 0.) 
;      ivar = fltarr(sz_img[0], sz_img[1]) 
;      ivar[gdvar] = 1. / (var[gdvar] + 25.)  ; Essential kludge for bluest data
;      bdpix = where(var EQ -1., nbd) 
;      if nbd NE 0 then ivar[bdpix] = -1.  ; CR rejection

      ;; Mask
      msk = lonarr(sz_img[0],sz_img[1])
      ;; Final image
      if flg_ordr EQ 1 then img_new = xmrdfits(imgfil, 2, /silent) $ 
      else img_new = fltarr(sz_img[0],sz_img[1])
      ;; Sky image
;      img_sky = fltarr(sz_img[0],sz_img[1])

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Loop on Orders
      for qq=ordr[0],ordr[1] do begin
          ;; 
          print, 'esi_echpskysub: Subtracting order ', $
            string(15L-qq, FORMAT='(i3)')

          ;; Create Mask
          msk[*] = 0
          lhs = (rnd_edg[*,qq,0]+17L) > 0L  ;; LHS has funny edge (0.5" only?)
          rhs = (rnd_edg[*,qq,1]-9L) < (sz_img[0]-1)

          gdmsk = where(lhs LE rhs, ngdmsk)
          for ii=0L,ngdmsk-1 do begin
              j = gdmsk[ii]
              msk[lhs[j]:rhs[j],j] = 1
          endfor


          ;; Mask out Obj
          print, 'esi_echpskysub: Masking... '
          mskobj = where(objstr.slit_id EQ qq, nmsk)
          for ii=0L,nmsk-1 do begin
              for jj=0L,sz_img[1]-1 do begin
                  imn = 0L > round(objstr[mskobj[ii]].trace[jj] $
                                  - objstr[mskobj[ii]].aper[0])
                  imx = (sz_img[0]-1) < round(objstr[mskobj[ii]].trace[jj] $
                                              + objstr[mskobj[ii]].aper[1]) > 0L
                  msk[imn:imx,jj] = 2
              endfor
          endfor

          ;; Bad Arc pix
          bdwv = where(img_arc LT all_mnxwv[qq,0] OR img_arc GT all_mnxwv[qq,1], $
                      nbdwv)
          if nbdwv NE 0 then msk[bdwv] = 0
          delvarx, bdwv

          ;; Bad Flux pix
          bdfx = where(var LE 0., nbdfx)
          if nbdfx NE 0 then msk[bdfx] = 0
          delvarx, bdfx


          ;; Special for Orders
          case qq of
              0: msk[*,0:1499] = 0L
              9: msk[*,2500:sz_img[1]-1] = 0L
              else:
          endcase

          ;; VARIANCE and inverse variance
          print, 'esi_echpskysub: Variance'
          newvar = fltarr(sz_img[0], sz_img[1])
          for ii=0L,ngdmsk-1 do begin
              j = gdmsk[ii]
              b = where(msk[lhs[j]:rhs[j],j] EQ 1, nb)
              if nb NE 0 then $
                newvar[lhs[j]:rhs[j],j] = median(img[lhs[j]+b,j]) > 0.
          endfor
          ;; RN  
          newvar = newvar + esi[indx[0]].readno^2
          ivar = fltarr(sz_img[0], sz_img[1]) 
          ivar = 1. / (newvar)
          bdpix = where(var EQ -1., nbd) 
          if nbd NE 0 then ivar[bdpix] = -1. ; CR rejection

          ;; CR
          bd = where(var GT 100. AND var GT 2*newvar AND msk EQ 1, nbd)
          delvarx, newvar
          if nbd NE 0 then begin
              msk[bd] = 0
              ivar[bd] = -1.
          endif

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; FIT POLY
          print, 'esi_echpskysub: Fitting... '
          if qq LT 4 then nord = 1 else nord = 2
          pfitstr = x_setfitstrct(NORD=nord, HSIG=5., LSIG=5., NITER=2L, FLGREJ=1)
          for j=0L,sz_img[1]-1 do begin
              a = where(msk[*,j] EQ 1, na)
              b = where(msk[*,j] NE 0, nb)
              if na GT 2 then begin
                  fit = x_fitrej(img_arc[a,j], img[a,j], FITSTR=pfitstr, $
                                 IVAR=ivar[a,j])
                  if fit[0] NE -1 then begin
                      img_new[b,j] = img[b,j] - $
                        x_calcfit(img_arc[b,j],FITSTR=pfitstr)
                  endif else var[b,j] = -1
              endif else begin
                  if nb NE 0 then var[b,j] = -1.
              endelse
          endfor

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
          ;; Create sky image
;          print, 'esi_echpskysub: Creating a sky image'
;          gdreg = where(msk NE 0 AND img_arc GT 0.)
;          img_sky[gdreg] = bspline_valu(img_arc[gdreg], bset)

          ;; Subtract from image
;          print, 'esi_echpskysub: Subtracting sky'
;          img_new[gdreg] = img[gdreg] - img_sky[gdreg]

      endfor
      ;; Ouptut New Image
      print, 'esi_echpskysub: Writing output to: ', imgfil
      if keyword_set( CHK ) then xatv, img_new, WVIMG=img_arc, /block
      mwrfits, img, imgfil, head, /create, /silent
      mwrfits, var, imgfil, /silent
      mwrfits, img_new, imgfil, /silent
      ;; COMPRESS
      spawn, 'gzip -f '+imgfil
  endfor
  
;  DONE
  print, 'esi_echpskysub: All done! '
  return
end
