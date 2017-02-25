;+ 
; NAME:
; esi_echold_skysub   
;     Version 1.1
;
; PURPOSE:
;    Sky Subtract image and add sky subtracted frame to the final
;    image (e.g. Final/f_esi0020.fits).  The program also outputs
;    information related to the sky fit in the Sky directory.  The
;    user has control over which order to subtract (ORDR=), the
;    wavelength image to use (AIMG=) and can input a break point file
;    for detailed sky line subtraction (SKLFIL=).
;
; CALLING SEQUENCE:
;   
;  esi_echold_skysub, esi, obj_id, [exp], /CHK, /STD, ORDR=, /NOVAC,
;  SKLFIL=, BORDR=, AIMG=
;
; INPUTS:
;   esi   -  ESI structure
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STD     - Sky subtract standard star
;  /CHK     - Show steps along the way
;  AIMG=    - Use alternate Wavelength map for sky subtraction
;             (string filename)
;  ORDR=    - Sky subtract a single order only (e.g.  [5L,5L])
;  /NOVAC   - Do not perform vacuum wavelength correction
;  BORDR=   - Order to begin bspline fitting (default: 5L)
;  SKLFIL=  - ASCII file setting breakpoints around sky lines (string)
;  /USEOLD  - Overwrite only the new orders into the old sky sub image
;  /FCHK
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echold_skysub, esi, 1L, [0L], /CHK, ORDR=7L   
;           {Sky sub exposure 0 and order 7 only}
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Jul-2002 Written by JXP
;   10-Oct-2002 Revised: Saving bset + allows single order redo
;   04-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echold_skysub, esi, obj_id, exp, CHK=chk, STD=std, ORDR=ordr $
                   , NOVAC = novac, SKLFIL = sklfil, BORDR = bordr $
                   , AIMG = aimg, USEOLD = useold, SEDG_FIL = sedg_fil $
                   , FITFIL = fitfil, CBIN = cbin, RBIN = rbin, FCHK = fchk $
                   , LSLITE = LSLITE, RSLITE = RSLITE

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echold_skysub, esi, obj_id, [exspr], /CHK, /STD, /NOVAC, ORDR='
      print, '          BORDR=, AIMG=, /USEOLD [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set(BORDR) then bordr = 5L  ;; Order to begin bspline
  if not keyword_set( FITFIL ) then fitfil = 'Maps/hole_fit.idl'
  if not keyword_set( LSLITE ) then lslite = round(22./cbin)
  if not keyword_set( RSLITE ) then rslite = round(14./cbin)
  if not keyword_set( ORDR ) then begin  ;; Orders to sky subtract
      ordr = [0L,9L]
      flg_ordr = 0
  endif else begin
      if keyword_set( USEOLD ) then flg_ordr = 1 else flg_ordr = 0
      if n_elements(ordr) NE 2 then begin
          print, 'esi_echold_skysub: ORDR must be a 2 element array'
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
      skylin[7,37,*] = [7771.0, 7777.0, 2]
      
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
      if x_chkfil(sklfil+'*') EQ 0 then begin
          print, 'esi_echold_skysub: Break point file doesnt exist! Returning...'
          return
      endif
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
                   esi.rbin EQ rbin AND esi.cbin EQ cbin AND $
                   esi.obj_id EQ obj_id AND $
                   strtrim(esi.type,2) EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'esi_echold_skysub: No images to sky subtract!', obj_id
          return
      endif
  endif else begin
      indx = [obj_id[0]]
      nindx = 1L
  endelse
      

;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

;  Read Arc
  ;; arc_fil
  if not keyword_set(AIMG) then $
    arc_fil = strtrim(esi[indx[0]].arc_fil,2) $
  else arc_fil = aimg
  ;; READ
  if x_chkfil(arc_fil+'*') EQ 0 then begin
      print, 'esi_echold_skysub: Arc file doesnt exist!', arc_fil
      return
  endif
  print, 'esi_echold_skysub: Reading Arc file ', arc_fil
  img_arc = xmrdfits(arc_fil, /silent) 
  sz_arc = size(img_arc, /dimensions)

; Open Slit file
  if not keyword_set( SEDG_FIL ) then begin
      slit_edg = esi_getfil('sedg_fil', SLIT=esi[indx[0]].slit, $
                            cbin=cbin, rbin=rbin)
  endif else begin
      slit_edg = xmrdfits(strtrim(SEDG_FIL,2), 0, /silent)
  endelse
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
          print, 'esi_echold_skysub: No Obj file! ', objfil, ' Skipping...'
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
      nobj = n_elements(objstr)/10

      ;; IMG+VAR Fil 
      imgfil = esi_getfil('fin_fil', subfil=esi[indx[exp[q]]].img_root,/name)
      skyfil = 'Sky/sky_'+esi[indx[exp[q]]].img_root
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'esi_echold_skysub: Image file doesnt exist!', imgfil
          return
      endif
      print, 'esi_echold_skysub: Reading Image files... ', imgfil
      img = xmrdfits(imgfil, 0, head, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Mask
      msk = lonarr(sz_img[0],sz_img[1])
      ;; Final image
      if flg_ordr EQ 1 then img_new = xmrdfits(imgfil, 2, /silent) $ 
      else img_new = fltarr(sz_img[0],sz_img[1])

      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
      ;; Loop on Orders
      for qq=ordr[0],ordr[1] do begin
          ;; 
          print, 'esi_echold_skysub: Subtracting order ', $
            string(15L-qq, FORMAT='(i3)')

          ;; Create Mask
          msk[*] = 0
          lhs = (rnd_edg[*,qq,0]+LSLITE) > 0L ;; LHS has funny edge (0.5" only?)
          rhs = (rnd_edg[*,qq,1]-RSLITE) < (sz_img[0]-1)

          gdmsk = where(lhs LE rhs, ngdmsk)
          for ii=0L,ngdmsk-1 do begin
              j = gdmsk[ii]
              msk[lhs[j]:rhs[j],j] = 1
          endfor
          ;; Mask out Obj
          print, 'esi_echold_skysub: Masking... '
          mskobj = where(objstr.order EQ qq, nmsk)
          for ii=0L,nmsk-1 do begin
              case qq of 
                 0: begin
                    jjmn = 1500L
                    jjmx = sz_img[1]-1L
                    ;;jjmx = 3800L ;; changed by JFH 04-29-2008
                 end
                 9: begin
                    jjmn = 0L
                    jjmx = 2170L
                 end
                 else: begin
                   jjmn = 0L
                   jjmx = sz_img[1]-1
                 end
              endcase
              for jj=jjmn,jjmx do begin
                  imn = 0L > round(objstr[mskobj[ii]].trace[jj] $
                                  - objstr[mskobj[ii]].aper[0])
                  imx = (sz_img[0]-1) < round(objstr[mskobj[ii]].trace[jj] $
                                              + objstr[mskobj[ii]].aper[1]) > 0L
                  if imn GT imx then begin
                      print, 'esi_echold_skysub: Trouble with masking!'
                      stop
                  endif
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
          print, 'esi_echold_skysub: Variance'
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

          ;; Sky Shape here!

          if qq GE BORDR then begin
              ;; BSPLINE
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Convert to 1D
              print, 'esi_echold_skysub: Grabbing sky pixels'
              skypix = where(msk EQ 1, nsky)
              srt = sort(img_arc[skypix])
              sky_wv = img_arc[skypix[srt]]
              sky_fx = img[skypix[srt]]
              sky_ivar = ivar[skypix[srt]]

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; BKPTS
              obj_wv = dblarr(sz_img[1])
              obj_mn = dblarr(sz_img[1])
              obj_mx = dblarr(sz_img[1])
              for j=0L,sz_img[1]-1 do begin
                  a = where(msk[*,j] EQ 1, na)
                  if na NE 0 then begin
                      obj_wv[j] = median(img_arc[a,j])
                      b = where(abs(img_arc[a,j]-obj_wv[j]) LT 2.,nb)
                      if nb GT 3 then begin
                          obj_mn[j] = min(img_arc[a[b],j], MAX=mx)
                          obj_mx[j] = mx
                      endif
                  endif
              endfor

              ;; Set bkpts
              subobj = obj_wv[where(obj_wv GT all_mnxwv[qq,0] AND $
                                    obj_wv LT all_mnxwv[qq,1], nobj)]
              bkpts = fltarr(nobj+1)
              bkpts[1:nobj] = (subobj + shift(subobj,1))/2.
              bkpts[0] = subobj[0] - (subobj[1]-subobj[0])/2.
              bkpts[nobj] = subobj[nobj-1] + (subobj[nobj-1]-subobj[nobj-2])/2.
              
              ;; Improve Sky lines by adding bkpts
              sz_slin = size(skylin, /dimensions)
              for k=0L,sz_slin[1]-1 do begin
                  case round(skylin[qq,k,2]) of
                      -1: 
                      0: 
                      1: begin ;; Center
                          addpt = where(subobj GT skylin[qq,k,0] AND $
                                        subobj LE skylin[qq,k,1],nadd)
                          if nadd NE 0 then bkpts = [bkpts, subobj[addpt]]
                      end
                      2: begin ;; Low/High wave
                          addpt = where(obj_mn GT skylin[qq,k,0] $
                                        AND obj_mn LE skylin[qq,k,1],nadd)
                          if nadd NE 0 then bkpts = [bkpts, obj_mn[addpt]]
                          addpt = where(obj_mx GT skylin[qq,k,0] $
                                        AND obj_mx LE skylin[qq,k,1],nadd)
                          if nadd NE 0 then bkpts = [bkpts, obj_mx[addpt]]
                      end
                      else: stop
                  endcase
              endfor
              ;; Add and Sort
              srt = sort(bkpts)
              bkpts = bkpts[srt]
              
              ;; Ends
              if bkpts[0] GT sky_wv[0] then bkpts=[sky_wv[0]-0.1,bkpts]
              if bkpts[n_elements(bkpts)-1] LT sky_wv[n_elements(sky_wv)-1] then $
                bkpts=[bkpts,sky_wv[n_elements(sky_wv)-1]+0.1]
              
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; FIT Bspline
              nord = 3L
              print, 'esi_echold_skysub: Fitting ', nsky, ' sky pixels with Bspline'
              bset = bspline_iterfit(sky_wv, sky_fx, bkpt=bkpts, nord=nord, $
                                     upper=5., lower=5., INVVAR=sky_ivar, $
                                     maxiter=5L)
              ;; Check FIT
              if keyword_set( CHK ) then begin
                  nfit = 100000L
                  x0 = all_mnxwv[qq,0]
                  xN = all_mnxwv[qq,1]
                  
                  xfit = fltarr(nfit)
                  for i=0L,nfit-1 do xfit[i] = x0 + $
                    float(i)*(xN-x0)/float(nfit)
                  yfit = bspline_valu(xfit, bset)
                  ybkpt = bspline_valu(bset.fullbkpt, bset)
                  a = where(msk EQ 2)
                  srt = sort(img_arc[a])
                  obj_w = img_arc[a[srt]]
                  obj_f = img[a[srt]]
                  x_splot, sky_wv, sky_fx, PSYM1=3, /block, XTWO=xfit, $
;                YTWO=yfit;, XTHR=obj_w, YTHR=obj_f, PSYM3=1
                  YTWO=yfit, XTHR=bset.fullbkpt, YTHR=ybkpt, PSYM3=2
              endif

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; Create sky image
              gdreg = where(msk NE 0 AND img_arc GT 0.)
;              print, 'esi_echold_skysub: Creating a sky image'
;              img_sky[gdreg] = bspline_valu(img_arc[gdreg], bset)
              
              ;; Subtract from image
              print, 'esi_echold_skysub: Subtracting sky'
              img_new[gdreg] = img[gdreg] - bspline_valu(img_arc[gdreg], bset)
              ;; Release memory
              delvarx, bkpts, gdmsk, gdreg
          endif else begin ;; POLY FIT
            ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
              ;; FIT POLY
              print, 'esi_echold_skysub: Fitting POLY... '
              if qq LT 4 then nord = 1 else nord = 2
              pfitstr = x_setfitstrct(NORD=nord, HSIG=5., LSIG=5., $
                                      NITER=2L, FLGREJ=1)
              sky_spec = fltarr(sz_img[1],2)
              for j=0L,sz_img[1]-1 do begin
                  a = where(msk[*,j] EQ 1, na)
                  b = where(msk[*,j] NE 0, nb)
                  if na GT 3 then begin
;                      fit = x_fitrej(img_arc[a,j], img[a,j], FITSTR=pfitstr, $
;                                     IVAR=ivar[a,j])
                     fit = x_fitrej(img_arc[a,j], img[a,j], FITSTR=pfitstr)
                     if fit[0] NE -1 then begin
                          img_new[b, j] = img[b, j] - $
                            x_calcfit(img_arc[b, j], FITSTR = pfitstr)
                          sky_spec[j, 0] = img_arc[b[nb/2], j]
                          sky_spec[j,1] = x_calcfit(sky_spec[j,0] $
                                                    , FITSTR = pfitstr)
                      endif else var[b,j] = -1
                  endif else begin
                      if nb NE 0 then var[b,j] = -1.
                  endelse
              endfor
          endelse
          
          
          ;; Write Sky info (bset or spectrum)
          if flg_ordr EQ 0 then begin  ;; ENTIRE SPECTRUM
              if qq GE BORDR then begin
                  if qq EQ 0 then mwrfits, bset, skyfil, /create, /silent $
                  else mwrfits, bset, skyfil, /silent 
              endif else begin
                  if qq EQ 0 then begin
                      mwrfits, fltarr(5), skyfil, /create, /silent 
                      mwrfits, sky_spec, skyfil, /silent  ; Twice for organization
                  endif else mwrfits, sky_spec, skyfil, /silent 
              endelse
          endif else begin  ;; SELECT ORDERS
              a = findfile(skyfil+'*',count=na)
              if na EQ 0 then stop
              spawn, 'cp '+skyfil+' tmp_sky.fits'
              case qq of 
                  0: begin
                      if qq GE BORDR then mwrfits, bset, skyfil, /create,/silent $
                      else begin
                          mwrfits, fltarr(5), skyfil, /create, /silent 
                          mwrfits, sky_spec, skyfil, /silent
                      endelse
                      for i=1L,9L do begin
                          tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                          mwrfits, tmp_sky, skyfil, /silent 
                      endfor
                  end
                  else: begin
                      tmp_sky = xmrdfits('tmp_sky.fits', 1, /silent)
                      ;; FIRST
                      if BORDR EQ 0 then mwrfits,tmp_sky,skyfil,/create,/silent $
                      else begin
                          mwrfits, fltarr(5), skyfil, /create, /silent 
                          mwrfits, tmp_sky, skyfil, /silent
                      endelse
                      ;; MIDDLE
                      for i=1L,qq-1 do begin
                          tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                          mwrfits, tmp_sky, skyfil, /silent 
                      endfor
                      ;; NEW
                      if qq GE BORDR then mwrfits, bset, skyfil, /silent $
                      else mwrfits, sky_spec, skyfil, /silent
                      ;; REST
                      for i=qq+1L,9 do begin
                          tmp_sky = xmrdfits('tmp_sky.fits', i+1, /silent)
                          mwrfits, tmp_sky, skyfil, /silent 
                      endfor
                  endelse
              endcase
              spawn, '\rm -f tmp_sky.fits'
          endelse
      
      endfor
      ;; Ouptut New Image
      print, 'esi_echold_skysub: Writing output to: ', imgfil
      if keyword_set( CHK ) or keyword_set( FCHK ) then $
        xatv, img_new, WVIMG=img_arc, /block, min=-50, max=200
      mwrfits, img, imgfil, head, /create, /silent
      mwrfits, var, imgfil, /silent
      mwrfits, img_new, imgfil, /silent
      ;; COMPRESS
      print, 'esi_echold_skysub: Compressing...'
      spawn, 'gzip -f '+imgfil
  endfor
  
;  DONE
  print, 'esi_echold_skysub: All done! '
  return
end
