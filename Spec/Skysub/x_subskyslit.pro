;+ 
; NAME:
; x_subskyslit   
;    Version 1.1
;
; PURPOSE:
;    Given the slitstr and slit#, subtract the sky and return 
;         the result.  Used with WFCCD (only?)
;
; CALLING SEQUENCE:
;  x_subskyslit, slit_fil, slit, obj_fil, flux, wave, VAR=, 
;                 SUBIMG=, PIX=, DATFIL=, ALL_RMS=, 
;                 SKYIMG=, REJPIX=, /DEBUG=
;
; INPUTS:
;   slit_fil  -- File for slit structure
;   slit      -- Index
;   obj_fil   -- File for object structure
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
;   20-Mar-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_subskyslit, slit_fil, slit, obj_fil, flux, wave, VAR=var, $
                  SUBIMG=subimg, PIX=pix, DATFIL=datfil, ALL_RMS=all_rms, $
                  SKYIMG=skyimg, REJPIX=rejpix, DEBUG=debug


;  Error catching
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_subskyslit, slit_fil, slit, obj_fil, [flux, wave], '
    print,  '        VAR=, SUBIMG=, PIX=, DATFIL=, /REJPIX [v1.0]'
    return
  endif 


;  Optional Keywords
  if not keyword_set( FLUX ) AND not keyword_set(DATFIL) then begin
      print, 'x_subskyslit: Need to define flux somehow!'
      stop
      return
  endif
  if not keyword_set( WAVE ) AND not keyword_set(DATFIL) then begin
      print, 'x_subskyslit: Need to define flux somehow!'
      stop
      return
  endif

  if not keyword_set( REJPER ) then rejper = 0.01

;  Read in the data and wave info
  
  if not keyword_set( FLUX ) then flux = mrdfits(datfil)
  if not keyword_set( WAVE ) then wave = mrdfits(datfil)

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
  for i=0L, sz[0]-1 do begin
      yedg = round(slitstr[slit].yedg_sky[i,*])
      msk[i,yedg[0]+1:yedg[1]-1] = 1
  endfor

;  Mask out Objects
  obj = where(objstr.slit_id EQ slitstr[slit].id, nobj)
  for i=0L,nobj-1 do begin
      ap_low = objstr[obj[i]].aper[0] < (-1.5)
      ap_high = objstr[obj[i]].aper[1] > 1.5
      for j=0L,sz[0]-1 do $
        msk[j,round(objstr[obj[i]].trace[j]+ap_low):$
            round(objstr[obj[i]].trace[j]+ap_high)] = 2
  endfor

; Reject bad pixels (includes slit edges)
  if keyword_set( VAR ) then begin
      badpix = where(var LE 0., nbad)
      if nbad NE 0 then msk[badpix] = 0
  endif


; Subimg
  subimg = fltarr(sz[0],sz[1])
  skyimg = fltarr(sz[0],sz[1])

; FIT STRUCTURES
  fitstr = { fitstrct }
  fitstr.func = 'LEGEND'
  fitstr.nord = 2
  fitstr.hsig = 1.5
  fitstr.lsig = 1.5
  fitstr.maxrej = (round(0.1*(slitstr[slit].yedg_flt[1]-$
                              slitstr[slit].yedg_flt[0]-3))) > 1
  fitstr.niter = 2
  fitstr.minpt = 3
  
  fitstr3 = fitstr
  fitstr3.nord = 3

; SET UP THE ARRAYS

  all_skypix = where(msk EQ 1)

  nsky = n_elements(all_skypix)
  sv_sky = fltarr(nsky)
  sv_wv = fltarr(nsky)
  all_rms = fltarr(sz[0])

;;;;;;;;;;;;; LOOP ON COLUMNS ;;;;;;;;
  for i=0L,sz[0]-1 do begin
      yedg = round(slitstr[slit].yedg_sky[i,*])
      ; FIND SKY PIXELS
      skypix = where(msk[i,yedg[0]:yedg[1]] EQ 1, nskypix)
      if nskypix EQ 0 then begin
          slitpix = where(msk[i,yedg[0]:yedg[1]] NE 0, nslitpix)
          if nslitpix NE 0 then begin
              subimg[i,slitpix+yedg[0]] = flux[i,slitpix+yedg[0]] 
              skyimg[i,slitpix+yedg[0]] = 0.
          endif
          continue
      endif
          
      sky_wv = reform(wave[i,skypix+yedg[0]])
      sky_fx = reform(flux[i,skypix+yedg[0]])

      ; SORT
;      srt = sort(sky_wv)
;      sky_wv = sky_wv[srt]
;      sky_fx = sky_fx[srt]

      ; FIT the Sky Pixels
      if nskypix LE 5 then begin  ; Deal with very few sky pixels
          ; MEDIAN
          sky = median(sky_fx)
          if nskypix GT 1 then $
          all_rms[i] = sqrt(total((sky_fx - replicate(sky,nskypix))^2 )$
                            /(nskypix-1.)) $
          else all_rms[i] = sqrt(sky)
      endif else begin 
          ; FIT
          sky = x_fitrej(sky_wv, sky_fx, FITSTR=fitstr, RMS=rms)
          if nskypix GT 15 then begin
              ; FIT 
              sky3 = x_fitrej(sky_wv, sky_fx, FITSTR=fitstr3, RMS=rms3) 
          ; Take nord=2 unless nord=3 is much better
              if (rms3-rms)/(rms+rms3)/2 LT -0.3 then begin
                  all_rms[i] = rms3
                  gdfit = fitstr3 
              endif else begin  ; NORD=2
                  gdfit = fitstr
                  all_rms[i] = rms
              endelse
          endif else begin
              gdfit = fitstr
              all_rms[i] = rms
          endelse 
          ; DEBUG
          if keyword_set( DEBUG ) AND i GE 780 then begin
              clr = getcolor(/load)
              plot, sky_wv, sky_fx, color=clr.black, background=clr.white,$
                psym=1
              oplot, sky_wv, sky, color=clr.red
              if keyword_set(sky3) then oplot, sky_wv, sky3, color=clr.blue
              stop
          endif
      endelse
      ; Calculate the sky at each pixel
      slitpix = where(msk[i,yedg[0]:yedg[1]] NE 0, nslitpix)
      if nskypix GT 5 then begin
          slit_wv = reform(wave[i,slitpix+yedg[0]])
          all_sky = x_calcfit(slit_wv, FITSTR=gdfit)
      endif else all_sky = replicate(sky,nslitpix)
      ; Subtract from all pix
      subimg[i,slitpix+yedg[0]] = flux[i,slitpix+yedg[0]] - all_sky
      skyimg[i,slitpix+yedg[0]] = all_sky
  endfor

  ; PIX
  pix = where(msk NE 0)
  return
end
      
