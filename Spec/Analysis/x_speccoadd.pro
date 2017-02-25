;+ 
; NAME:
; x_speccoadd
;   Version 1.1
;
; PURPOSE:
;    Routine written to coadd spectra.  I highly recommend using
;    x_combspec instead.  The routine here is old and not recently
;    tested!
;
; CALLING SEQUENCE:
;   dat = x_speccoadd(spec)
;
; INPUTS:
;   spec       - Fits file or data
;
; RETURNS:
;   dat       - Data in fits file or the input data
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  FSCALE      - Data is float
;  DSCALE      - Data is float
;
; OPTIONAL OUTPUTS:
;  HEAD        - Header
;
; COMMENTS:
;
; EXAMPLES:
;   dat = x_speccoadd('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-Feb-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_speccoadd, files, fmt, delt, outfil=outfil

;
  if  N_params() LT 3  then begin 
    print,'Syntax - ' + $
             'x_speccoadd, indata, fmt, delt, OUTFIL=[V1.1]'
    return
  endif 

  stop  ;; This code is obsolete

;  Optional Keywords

  if not keyword_set( OUTFIL ) then outfil = 'tmp.fits'

;

  if delt LT 0. then begin
      flg_vel = 1 
      delt = -1.* delt
  endif else flg_vel = 0

;  Read in the data

  print, 'x_speccoadd: Reading the files...'
  case fmt of 
      1: begin  ; Fits files (flux, sig, wave)
          nfil = n_elements(files)
          flux = fltarr(200000L, nfil)
          sig = fltarr(200000L, nfil)
          wave = dblarr(200000L, nfil)
          allwvmn = dblarr(nfil)
          allwvmx = dblarr(nfil)
          npix = lonarr(nfil)
          for q=0L,nfil-1 do begin
              sxhread,files[q],head
              npix[q] = fxpar(head, 'NAXIS1')
              flux[0:npix[q]-1,q] = mrdfits(files[q], 0, head, /silent)
              sig[0:npix[q]-1,q] = mrdfits(files[q], 1, /silent)
              wave[0:npix[q]-1,q] = mrdfits(files[q], 2, /silent, /dscale)
              ; SORT!
              srt = sort(wave[0:npix[q]-1,q])
              flux[0:npix[q]-1,q] = flux[srt,q]
              sig[0:npix[q]-1,q] = sig[srt,q]
              wave[0:npix[q]-1,q] = wave[srt,q]
              ; wvmin, wvmax
              allwvmn[q] = wave[0,q]
              allwvmx[q] = wave[npix[q]-1,q]
          endfor
      end
  endcase

;  Create variance
  if not keyword_set( var ) then begin
      var = sig*0.
      gdsig = where(sig GT 0.)
      var[gdsig] = (sig[gdsig])^2
  endif

;  Find endpoints
  wvmn = min(allwvmn)
  wvmx = max(allwvmx)

  print, 'x_speccoadd: Creating final wave array with wvmnx=', wvmn, wvmx
; Set npix, bwv, fwv

  if flg_vel EQ 0 then begin  ; Linear binning
      nwpix = round((wvmx-wvmn)/delt) + 1L
      bwv = wvmn - delt/2. + dindgen(nwpix+1)*delt
      fwv = wvmn + dindgen(nwpix)*delt
  endif

; Create tmp arrays for final array
  tmpfx = fltarr(nwpix)
  tmpfrac = fltarr(nwpix)
  tmpvar = fltarr(nwpix)
  svfx = fltarr(nwpix,nfil)
  svwt = fltarr(nwpix,nfil)
  finfx = fltarr(nwpix)
  finsig = fltarr(nwpix)
  finwt = fltarr(nwpix)

;  Loop!

  print, 'x_speccoadd: Big Loop!'
  for jj=0L, nfil-1 do begin
      print, 'x_speccoadd: Rebinning Image ', jj, ' of ', nfil-1
      ; Reset
      tmpfx[*] = 0.
      tmpfrac[*] = 0.
      tmpvar[*] = 0.
      ; Loop on pixels
      for q=0L, npix[jj]-1 do begin
          ; Dont include bad pixels!
          if var[q,jj] LE 0 then continue
          ; Set lower wavelength of pixel q
          if q GT 0L then lw = (wave[q,jj] + wave[q-1,jj])/2. $
            else lw = wave[q,jj] + (wave[q,jj]-wave[q+1,jj])/2.
          ; Set higher wavelength of pixel q
          if q LT npix[jj]-1 then hw = (wave[q,jj] + wave[q+1,jj])/2. $
            else hw = wave[q,jj] + (wave[q,jj]-wave[q-1,jj])/2.
      
          ; Fill up the Final array
          ; No overlap
          if hw LE bwv[0] OR lw GE bwv[nwpix] then continue
          ; Find pixel that bw is within
          if lw LT bwv[0] then i1 = 0 else $
            i1 = where(lw LE shift(bwv,-1) AND lw GT bwv)
          ; Same for hw
          if hw GT bwv[nwpix] then i2 = nwpix-1 else $
            i2 = where(hw LE shift(bwv,-1) AND hw GT bwv)
          j1 = i1[0]
          j2 = i2[0]

          ; Now Sum up
          for kk=j1,j2 do begin
              frac = ( (hw < bwv[kk+1]) - (lw > bwv[kk]) ) / (hw-lw)
              tmpfx[kk] = tmpfx[kk] + frac * flux[q,jj]
	              tmpvar[kk] = tmpvar[kk] + frac * var[q,jj]  ; Only weight by frac!!!!
              tmpfrac[kk] = tmpfrac[kk] + frac  
          endfor
      endfor
      ; Add it up!
      for kk=0L, nwpix-1 do begin
          if tmpvar[kk] GT 0. AND tmpfrac[kk] GT 0. then begin
              ftmp = tmpfx[kk] / tmpfrac[kk]
              vtmp = tmpvar[kk] / (tmpfrac[kk])^2
              svfx[kk,jj] = ftmp
              svwt[kk,jj] = 1./vtmp
          endif
      endfor
  endfor

; Final Values!	

  print, 'x_speccoadd: Final Weighting'
  for kk=0L, nwpix-1 do begin
      finfx[kk] = 0.
      finwt[kk] = 0.
      finsig[kk] = -1.
      ; Loop on the spectra
      for jj=0L,nfil-1 do begin
          if svwt[kk,jj] GT 0. then begin
              finfx[kk] = finfx[kk] + svfx[kk,jj] * svwt[kk,jj]
              finwt[kk] = finwt[kk] + svwt[kk,jj]
          endif 
      endfor
  endfor
  gdpx = where(finwt NE 0.)
  finfx[gdpx] = finfx[gdpx] / finwt[gdpx]
  finsig[gdpx] = 1./ sqrt(finwt[gdpx])

; Outfil

  print, 'x_speccoadd: Writing to ', outfil
  mwrfits, finfx, outfil, head, /create
  mwrfits, finsig, outfil
  mwrfits, fwv, outfil

  return
end
