;+ 
; NAME:
; x_crossarc
;     Version 1.1
;
; PURPOSE:
;    Finds shift between saved template and current arc spectra. 
;
; CALLING SEQUENCE:
;  x_crossarc, guessarc, cur_aspec, guess_spec, guess_fit, ordrs, $
;         ordr_shift, row_shift, chk=chk, sigrej=sigrej, /DEBUG
;
; INPUTS:
;   guessarc  -  Filename of IDL save file with known wavelength solution
;   cur_aspec -  Current extracted Arc spectra to be fit

;   obj_id   -  Object identifier
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;   guess_spec -  saved arc spectrum rebinned to match cur_aspec
;   guess_fit  -  polynomial wavelength fits in rebinned space
;   ordrs      -  Orders matching saved spectrum
;   ordr_shift -  Order offset between template and current
;   row_shift  -  Pixel offset between template and current
;
; OPTIONAL KEYWORDS:
;   /CHK      - Manually check steps along the way
;   /DEBUG    - Debugging
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;   restore
;   fft()
;
; REVISION HISTORY:
;   2004 Written by SB
;-
;------------------------------------------------------------------------------
pro x_crossarc, guessarc, in_cur_aspec, guess_spec, guess_fit, ordrs, $
                ordr_shift, row_shift, chk=chk, sigrej=sigrej, DEBUG=debug, $
                CCDSZ=ccdsz, ONLYG=onlyg, TRIM_FFT=trim_fft

  if  N_params() LT 7  then begin 
      print,'Syntax - x_crossarc, guessarc, cur_aspec, guess_spec, ' + $
        'guess_fit, ordrs, '+$
        'ordr_shift, row_shift, chk=, sigrej=, DEBUG=, CCDSZ=, /ONLYG [v1.0]'
      return
  endif 
  
   if NOT keyword_set(SATUR) then satur = 30000.
   if NOT keyword_set(sigrej) then sigrej=2
   if NOT keyword_set(CCDSZ) then ccdsz=[2048L,4096]

   restore, guessarc
   ;; Rename
   guess_spec = temporary(sv_aspec)
   guess_fit = all_arcfit
   guess_fit.hsig = sigrej
   guess_fit.lsig = sigrej
   guess_fit.flg_rej = 1
   guess_fit.niter = 3
   ordrs = guess_ordr

   ;; KLUDGE2
   ;CHK=1
   ;;cur_aspec = cur_aspec[*,45:126]
   ;;cur_aspec = cur_aspec[*,48+lindgen(77)]
   ;cur_aspec = cur_aspec[*,22+lindgen(77)]
   ;;guess_spec = guess_spec[*, 20:*]

   ncol = (size(guess_spec))[1]
   nrow = (size(guess_spec))[2]

   if keyword_set(TRIM_FFT) then begin
      szc = size(in_cur_aspec, /dimen)
      if nrow GT szc[1] then stop ;; Try another kludge!
      cur_aspec = in_cur_aspec[*, 0:nrow-1]
   endif else cur_aspec = in_cur_aspec
   

   ccol = (size(cur_aspec))[1]
   cbin = round(float(ccdsz[1])/ccol)

   

   ;;    Need a resize?
   if ccol  NE ncol then begin
      print, 'x_crossarc: Resizing the array'
       temp_spec = rebin(guess_spec,ccdsz[1], nrow)
       guess_spec = rebin(temp_spec[0:cbin*ccol-1,*], ccol, nrow)
       change_lines = 1
   endif

   if keyword_set( ONLYG ) then begin
       row_shift =0L
       ordr_shift = 0L
       cc_peak = 1000.
       cc_err = 1.
   endif else begin

   ;; Cross correlate (FFT) against the archived 1D arc spectrum
;;   guess_fft = fft([guess_spec[*], guess_spec[*]*0])
;;  cur_fft    = fft([cur_aspec[*], cur_aspec[*]*0])
;; cc        = float(fft( cur_fft * conj(guess_fft) ))
       sg = smooth(guess_spec[*] < SATUR,5)  
       guess_fft = fft([sg, sg*0])
       cg = (smooth(cur_aspec[*] < SATUR,5)  ) > 0
       cur_fft    = fft([cg , cg*0])
       cc        = float(fft( cur_fft * conj(guess_fft) ))
       cc_smooth = shift(smooth(shift(cc, n_elements(cg)), 4*ncol), $
                         -1L*n_elements(cg))
       ;stop
       cc = cc - cc_smooth

       if keyword_set( DEBUG ) then stop
       ntot_fft = n_elements(cc)
       cc_peak = max(cc,total_shift)
       cc_err = sqrt(mean(cc^2))
       print, 'x_crossarc: Max peak in full correlation ', $
              cc_peak,' compared to error of ', cc_err
   endelse

   ;; UVES KLUDGE!
;   stop
;   total_shift = -2940L
;   cc_err = 0.000001
   

   ;; MIKE specific
   ;KLUDGE=1
   if (cc_peak LT cc_err*100) and not keyword_set(KLUDGE) then begin  ;; JXP debug
       print, 'Maybe using New blue CCD, dispersion direction is reversed'
       rev_guess = reverse(guess_spec)
       rfft      = fft([rev_guess[*], rev_guess[*]*0])
       rev_cc    = float(fft( cur_fft * conj(rfft) ))
       rev_cc_smooth = shift(smooth(shift(rev_cc, n_elements(cg)), 4*ncol), $
                             -1L*n_elements(cg))
       rev_cc = rev_cc - rev_cc_smooth
       
       ntot_fft = n_elements(rev_cc)
       cc_rev_peak = max(rev_cc,total_shift_rev)
       cc_rev_err = sqrt(mean(rev_cc^2))
       print, 'Max peak in full correlation ',cc_rev_peak,  $
         ' compared to error of ', cc_rev_err
       
       if cc_rev_peak GT 10.*cc_rev_err AND $
         cc_rev_peak*cc_err GT cc_peak*cc_rev_err then begin
           print, 'yes the dispersion orientation is reversed'
           print, 'swapping guess_spec and guess_fit'
           guess_spec = rev_guess
           total_shift = total_shift_rev
           change_lines = -1
       endif else print, 'Nope keeping normal direction'
   endif

   if keyword_set(change_lines) then begin
       for i=0, n_elements(sv_lines)-1 do begin
           ngd = sv_lines[i].nlin
           if ngd GT 0 then begin
               opix = ccdsz[1]*sv_lines[i].pix[0:ngd-1]/ncol
               
               sv_lines[i].pix[0:ngd-1] = change_lines EQ -1 ? $
                 (ccdsz[1]/cbin-1) -opix/cbin : opix/cbin 
               
           endif
       endfor 
       
       for i=0L, n_elements(guess_fit)-1 do begin
           if ptr_valid(guess_fit[i].ffit) then begin
               ffit = *(guess_fit[i].ffit)
               guess_fit[i].nrm[0] = ccdsz[1] * guess_fit[i].nrm[0]/ncol/cbin
               guess_fit[i].nrm[1] = ccdsz[1] * guess_fit[i].nrm[1]/ncol/cbin
               if change_lines EQ -1 then begin
                   ffit = ffit * (1.0 -2.0 * (lindgen(20) mod 2))
                   *(guess_fit[i].ffit) = ffit
                   guess_fit[i].nrm[0] = (ccdsz[1]-1)/cbin $
                     - guess_fit[i].nrm[0]
               endif 
           endif
       endfor
   endif

   if keyword_set(ONLYG) then return

   if total_shift GT ntot_fft/2 then total_shift = total_shift - ntot_fft

   ;; Find the shifts between the image and the archived image
   ordr_shift = round(double(total_shift) / ccol)
   row_shift = total_shift - ordr_shift*ccol
   
   if keyword_set(chk) then begin
       plot, lindgen(ntot_fft)-ntot_fft/2, shift(cc,ntot_fft/2), $
         xr=[total_shift-100, total_shift+100]
       oplot, [total_shift,total_shift], [0,10*cc_peak]
       xyouts, [total_shift,total_shift]+20, [0.9,0.8]*cc_peak, $
         strcompress(['Orders offset by '+ string(ordr_shift), $
                      ' And pixel offset by '+string(row_shift)])
    endif
   
   return
end 
