;+ 
; NAME:
; x_setslits   
;    Version 1.0
;
; PURPOSE:
;    Finds slit edges on the 'flattened' flat
;      WARNING -- This program has only been tested for the WFCCD
;
; CALLING SEQUENCE:
;   
;   x_setslits, slitstr, flatimg
;
; INPUTS:
;   slitstr     - Slit structure
;   flat        - Flat image or fits file
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   x_setslits, slitstr, flatimg
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_setslits, flat, slitstr, YSTRT=ystrt, DEBUG=debug


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'x_setslits, flat, slitstr [v1.0]'
    return
  endif 


;  Optional Keywords

  if not keyword_set( YSTRT ) then ystrt = 400L
  if not keyword_set( MXSHFT ) then mxshft = 150L

; Allow flat to be fits file

  dat = x_readimg(flat, /fscale)
  sz_img = size(dat, /dimensions)

; Reset yedg for median

  nslit = n_elements(slitstr)

; Create the sawtooth

  smsh = djs_median(dat[ystrt-10:ystrt+10,*],1)
  saw = shift(smsh,1) - shift(smsh,-1)
  npix = n_elements(saw)

; Create a 'spectrum' of the slits
  slit_spec = fltarr(sz_img[1])
  for q=0L,nslit-1 do begin
      tp_slit = round(slitstr[q].yedg[1]) < (sz_img[1]-4)
      ; Make peak
      slit_spec[tp_slit-1] = 0.3
      slit_spec[tp_slit] = 1.
      slit_spec[tp_slit+1] = 0.3
      ; Bottom
      bt_slit = 0 > round(slitstr[q].yedg[0])
      ; Make valley
      slit_spec[bt_slit] = -0.3
      slit_spec[bt_slit-1] = -1.
      slit_spec[bt_slit-2] = -0.3
  endfor
  
; Trim the edges
  saw[0:5] = 0.
  saw[npix-5:npix-1] = 0. 

; FFT
  fft_saw = fft(saw)
  fft_spec= fft(slit_spec)

; CORRELATE
  corr = fft( fft_saw * conj(fft_spec), /inverse)

; ANS
  ans = double(corr)

;;;;;;;;
;; DEBUG
  if keyword_set( DEBUG ) then begin
      x_splot, ans, /block
      stop
  endif

; Get SHIFT
  mx1 = max(ans[0:mxshft-1], shift1)
  mx2 = max(ans[sz_img[1]-mxshft-1:sz_img[1]-1], shift2)
  if mx2 GT mx1 then shift = shift2+sz_img[1]-mxshft-1 else shift = shift1
  if shift GT (sz_img[1]/2) then shift = shift - sz_img[1]

; Find the peaks
  x_fndpeaks, abs(saw), center, NSIG=50., PEAK=peak, /thin, $
    /force, pkwdth=1, /silent
;  tsaw = saw
;  tpeak = peak
;  x_splot, abs(saw), XTWO=center, YTWO=abs(saw[peak]), PSYM_Y2=2, /block
;  x_fndpeaks, abs(saw), center, NSIG=50., PEAK=peak, /thin, /norecent
;  npk = n_elements(peak)

; Sort the slit tops from top to bottom
;  srt = sort(slitstr.yedg[1])
;  srt = reverse(srt)
  
; Reverse peak as necessary (and flip saw)
;  if abs(peak[0]-slitstr[srt[0]].yedg[1]) GT $
;    abs(peak[npk-1]-slitstr[srt[0]].yedg[1]) then begin
;;      peak = reverse(peak)
;      saw = -saw
;  endif

; Find tops and bottoms
  top = peak[where(saw[peak] GT 0, ntop)]
  bot = peak[where(saw[peak] LT 0, nbot)]

  ; Padding
  bot = [-10000, bot, -10000, -10000]
  top = [-10000, top, -10000, -10000]
  
; LOOP ON SLITS

  dslit = fltarr(9)
  cslit = fltarr(9)
  for ii=0L,nslit-1 do begin
      ; Get nearest top and bot
      mn = min(abs(slitstr[ii].yedg[1]+shift - top), itop)
      mn = min(abs(slitstr[ii].yedg[0]+shift - bot), ibot)

      ; D values
      dslit[0] = top[itop]-bot[ibot]
      dslit[1] = top[itop]-bot[ibot+1]  
      dslit[2] = top[itop+1]-bot[ibot]
      dslit[3] = top[itop-1]-bot[ibot]
      dslit[4] = top[itop]-bot[ibot-1]
      dslit[5] = top[itop+1]-bot[ibot-1]
      dslit[6] = top[itop-1]-bot[ibot+1]
      dslit[7] = top[itop+1]-bot[ibot+1]
      dslit[8] = top[itop-1]-bot[ibot-1]

      ; Centers
      cslit[0] = (top[itop]+bot[ibot])/2.
      cslit[1] = (top[itop]+bot[ibot+1])/2.  
      cslit[2] = (top[itop+1]+bot[ibot])/2.
      cslit[3] = (top[itop-1]+bot[ibot])/2.
      cslit[4] = (top[itop]+bot[ibot-1])/2.
      dslit[5] = (top[itop+1]-bot[ibot-1])/2.
      dslit[6] = (top[itop-1]-bot[ibot+1])/2.
      dslit[7] = (top[itop+1]-bot[ibot+1])/2.
      dslit[8] = (top[itop-1]-bot[ibot-1])/2.

      ; Slit info
      cs = (slitstr[ii].yedg[1]+slitstr[ii].yedg[0])/2. + shift
      ds = slitstr[ii].yedg[1]-slitstr[ii].yedg[0]

      ; Good cen
      gdcen = where(abs(cslit-cs) LT 7., ngd)
      if ngd NE 0 then begin
          ; MIN D
          mn = min(abs(dslit[gdcen]-ds), best)
          ; Set
          slitstr[ii].yedg_flt[1] = cslit[gdcen[best]] + ds/2.
          slitstr[ii].yedg_flt[0] = cslit[gdcen[best]] - ds/2.
      endif else begin
          print, 'x_setslits: No good slit.  Better check!'
          slitstr[ii].yedg_flt[1] = cs + ds/2.
          slitstr[ii].yedg_flt[0] = cs - ds/2.
      endelse
  endfor

  return
end

