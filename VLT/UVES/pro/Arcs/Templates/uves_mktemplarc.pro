;+ 
; NAME:
; uves_mktemplarc   
;     Version 1.1
;
; PURPOSE:
;    Process (bias subtract, flatten) the arc files.  In addition, the
;    code first calls uves_arcxyoff which uses a 2D FFT to determine
;    the offset between the Template Arc and the current arc due to
;    thermal expansion in the instrument.
;   
;    This file also includes the routine uves_procarc_sngl which
;    allows the processing of a single Arc image given the filename of
;    the Raw image.
;
; CALLING SEQUENCE:
;   
;  uves_procarc, uves, setup, obj_id, [chip], ATEMPL=, /CLOBBER,
;  FLATFIL=
;
; INPUTS:
;   uves     -  MIKE structure
;   setup    -  Integer defining setup
;   obj_id   -  Object identifier
;   [chip]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  One processed, combined Arc Image (e.g. Arcs/Arc_mb0539.fits)
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite exisiting arc image
;   ATEMPL -- Index of the Template Arc image (default: 0L)
;   FLATFIL -- Filename of the milky flat (pixel to pixel correction)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   uves_procarc, uves, 1, 1, 1, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_getfil
;  uves_arcxyoff
;  uves_proc
;
; REVISION HISTORY:
;   13-Aug-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_mktemplarc, root, ordrs, outfil, TRIM=trim, FLIPO=flipo, NORV=norv

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'uves_mktemplarc, root, ordrs, outfil, TRIM=  [v1.1]'
      return
  endif 
  
  ;; Data files
  fil = findfile(root+'*', count=nfil)
  if nfil EQ 0 then begin
      print, 'uves_mktemplarc: No files!'
      return
  end
  sign = 1
  if ordrs[1] LT ordrs[0] then sign = -1
  guess_ordr = ordrs[0] + sign*lindgen(nfil)

  all_arcfit = replicate({fitstrct},50)
  if keyword_set(FLIPO) then begin
      i0 = nfil-1
      i1 = 0L
      is = -1
  endif else begin
      i0 = 0L
      i1 = nfil-1
      is = 1
  endelse

  ;; Loop
  jj = 0L
  for qq=i0,i1,is do begin
      ;; Read
      readcol, fil[qq], idx, wv, fx, format='L,D,F'
      if keyword_set(TRIM) then begin
          wv = wv[0:trim-1]
          fx = fx[0:trim-1]
      endif
      if qq EQ i0 then begin
          npix = n_elements(wv)
          sv_aspec = fltarr(npix,50)
      endif

      ;; Reverse
      if not keyword_set(NORV) then begin
          wv = reverse(wv)
          fx = reverse(fx)
      endif 

      ;; Fit
      fitstr = x_setfitstrct(NORD=4, LSIG=2., HSIG=2., NITER=3, $
                             MINPT=5, MAXREJ=10, /FLGREJ)
      fit = x_fitrej(idx, alog10(wv), fitstr=fitstr)

      ;; Save
      all_arcfit[jj] = fitstr
      sv_aspec[*,jj] = fx
      jj = jj + 1 
  endfor

  save, sv_aspec, all_arcfit, guess_ordr, filename=outfil
      
  ;; ALL DONE
  print, 'uves_procarc: All Done! ', outfil
  return

end

