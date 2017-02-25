;+ 
; NAME:
; x_bsplnsky   
;    Version 1.1
;
; PURPOSE:
;    Given the slitstr and slit#, subtract the sky and return 
;         the result using a bspline technique
;
; CALLING SEQUENCE:
;   x_bsplnsky, slit_fil, slit, obj_fil, [flux, wave], VAR=, $
;               SUBIMG=, PIX=, DATFIL=, ALL_RMS=, SKYIMG=, /CHK, WVMNX=
;
; INPUTS:
;   slit_fil   - File containing the slit structure
;   slit       - Slit to sky subtract
;   obj_fil    - File containing the object structure
;  [flux]      - Flux image
;  [wave]      - Wavelength image
;
; RETURNS:
;
; OUTPUTS:
;   Updates slitstr for original positions
;
; OPTIONAL KEYWORDS:
;  DATFIL=  -- FITS File containing FLUX, WAVE images
;  WVMNX=   -- Limit sky subtraction to specific wavelength range
;              [default [3900., 11000.] ]
;
; OPTIONAL OUTPUTS:
;  ALL_RMS=  -- Total RMS of the fit
;  PIX=      -- Array of pixel values 
;
; COMMENTS:
;
; EXAMPLES:
;   x_skysubslit, slitstr, map
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   10-Sep-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_bsplnsky, slit_fil, slit, obj_fil, flux, wave, VAR=var, $
                  SUBIMG=subimg, PIX=pix, DATFIL=datfil, ALL_RMS=all_rms, $
                  SKYIMG=skyimg, CHK=chk, $
                WVMNX=wvmnx


;  Error catching
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_bsplnsky, slit_fil, slit, obj_fil, [flux, wave], '
    print,  '        VAR=, SUBIMG=, PIX=, DATFIL=, ALL_RMS=  [v1.1]'
    return
  endif 


;  Optional Keywords
  if not keyword_set( FLUX ) AND not keyword_set(DATFIL) then begin
      print, 'x_bsplnsky: Need to define flux somehow!'
      stop
      return
  endif
  if not keyword_set( WAVE ) AND not keyword_set(DATFIL) then begin
      print, 'x_bsplnsky: Need to define flux somehow!'
      stop
      return
  endif

  if not keyword_set( GAIN ) then gain = 1.
  if not keyword_set( RN ) then rn = 5.
  if not keyword_set( WVMNX ) then wvmnx = [3900., 11000.]

;  SKYLIN

  skylin = fltarr(50,3)
;  skylin[0,*] = [5563, 5588, 3]
;  skylin[1,*] = [7700, 8000, 1]
;  skylin[2,*] = [8250, 9100, 1]
;  skylin[3,*] = [9350, 10100, 1]
;  Read in the data and wave info
  
  if not keyword_set( FLUX ) then flux = xmrdfits(datfil)
  if not keyword_set( WAVE ) then wave = xmrdfits(datfil,1,/silent)

;  Read in slitstruct
  if size(slit_fil,/type) EQ 7 then $
    slitstr = mrdfits(slit_fil, 1, STRUCTYP='mslitstrct', /silent) $
  else slitstr = slit_fil

;  Read in objstruct
  if size(obj_fil,/type) EQ 7 then $
    objstr = mrdfits(obj_fil, 1, STRUCTYP='specobjstrct', /silent) $
  else objstr = obj_fil

;  Find all pixels in that slit
  sz = size(flux, /dimensions)
  msk = bytarr(sz[0],sz[1])
  rnd_yedg = round(slitstr[slit].yedg_sky[*,*])
  for i=0L, sz[0]-1 do begin
      msk[i,rnd_yedg[i,0]+1:rnd_yedg[i,1]-1] = 1
  endfor

;  Mask out Objects
  obj = where(objstr.slit_id EQ slitstr[slit].id, nobj)
  for i=0L,nobj-1 do begin
      ap_low = objstr[obj[i]].aper[0] < (-1.5)
      ap_high = objstr[obj[i]].aper[1] > 1.5
      for j=0L,sz[0]-1 do begin
          ylow = round(objstr[obj[i]].trace[j]+ap_low) > rnd_yedg[j,0]
          yhigh= round(objstr[obj[i]].trace[j]+ap_high) < rnd_yedg[j,1]
          if yhigh GE ylow then msk[j,ylow:yhigh] = 2
      endfor
  endfor

; Reject bad pixels (includes slit edges)
  if keyword_set( VAR ) then begin
      badpix = where(var LE 0., COMPLEMENT=gd, nbad)
      if nbad NE 0 then msk[badpix] = 0
      ivar = fltarr(sz[0],sz[1]) -1.
      ivar[gd] = 1./var[gd]
  endif else begin
      var = (flux > 0.)*GAIN + RN^2
      ivar = 1./var
  endelse

; Set IVAR


; Subimg

  all_skypix = where(msk EQ 1)

; Set skypix
      
  sv_wv = wave[all_skypix]
  sv_sky = flux[all_skypix]
  sv_ivar= ivar[all_skypix]
  srt = sort(sv_wv)
  sky_wv = sv_wv[srt]
  sky_fx = sv_sky[srt]
  sky_ivar = sv_ivar[srt]

  ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
  ;; SETUP BKPTS
  print, 'x_bsplnsky: Setting Break points'

  wvoned = fltarr(sz[0])
  ycen = round(total(slitstr[slit].yedg_sky[*,*],2)/2.)
  for i=0L,sz[0]-1 do wvoned[i] = wave[i,ycen[i]]
  dwv = wvoned - shift(wvoned,1)
  a = where(dwv[1:sz[0]-1] GE 0.,na)
  if na NE 0 then bdpx = a[0]+1 else bdpx = sz[0]-1
      
  ;; Every in between pix
  hlf = (wvoned+shift(wvoned,1))/2.
  bkpts = hlf[1:bdpx-1]
  ;; Good ones
  gdbk = where(bkpts GT wvmnx[0] AND bkpts LT wvmnx[1])
  bkpts = bkpts[gdbk]
  ;; Improve sky lines
  sz_slin = size(skylin, /dimensions)
;  for k=0L,sz_slin[0]-1 do begin
;      case round(skylin[k,2]) of
;          -1: 
;          0: 
;          1: begin ;; Every 5 Ang
;              bkpts = [bkpts, 5.*findgen(round(skylin[k,1]-skylin[k,0])/5)+ $
;                       skylin[k,0]]
;          end
;          2: begin ;; Every ang
;              bkpts = [bkpts, findgen(round(skylin[k,1]-skylin[k,0])+1)+ $
;                       skylin[k,0]]
;          end
;          3: begin ;; Every 1/2 ang
;              bkpts = [bkpts, 0.5*findgen(round(skylin[k,1]-skylin[k,0])*2+2)+ $
;                       skylin[k,0]]
;          end
;          else: stop
;      endcase
;  endfor
  ;; Sort
  srt = sort(bkpts)
  bkpts = bkpts[srt]

;  x_splot, fin_wv, fin_fx, /block

; BSPLIN
  bset = bspline_iterfit(sky_wv, sky_fx, bkpt=bkpts, nord=3, upper=2.5, lower=2.5, $
                        INVVAR=sky_ivar, OUTMASK=outmsk)

  if keyword_set(CHK) then begin
      nfit = 50000L
      x0 = sky_wv[0] > 3300.
      xN = sky_wv[n_elements(sky_wv)-1] < 11000.
      xfit = fltarr(nfit)
      for i=0L,nfit-1 do xfit[i] = x0 + $
        float(i)*(xN-x0)/float(nfit)
      yfit = bspline_valu(xfit, bset)
      ybkpt = bspline_valu(bset.fullbkpt, bset)
      x_splot, sky_wv, sky_fx, PSYM1=2, /block, XTWO=xfit, YTWO=yfit, $
                XTHR=bset.fullbkpt, YTHR=ybkpt, PSYM_Y3=2
  endif

; SUBTRACT THE SKY

  ;; PIX
  pix = where(msk NE 0)

  ;; SKY
  sky = bspline_valu(wave[pix], bset)
  
  ;; Final image
  subimg = fltarr(sz[0], sz[1])
  subimg[pix] = flux[pix] - sky

  ;; RMS
  bfit = bspline_valu(sky_wv[where(outmsk EQ 1)], bset)
  all_rms = sqrt(total( (bfit-sky_fx[where(outmsk EQ 1)])^2 ) / $
                 (n_elements(bfit)-1.) )

  return
end
