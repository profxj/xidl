;+ 
; NAME:
; hires_thruput
;    Version 1.1
;
; PURPOSE:
;    Estimates the throughput of HIRES using empirical measurements
;
; CALLING SEQUENCE:
;  thru = hires_thruput( wave, center, iorder, fsr, ICD=)
;
; INPUTS:
;  wave=  -- Wavelength to calculate throughputs at
;  center=  -- Blaze centroid
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;  OUTDIR=  -- Name of output directory
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  showfits
;  querydss
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP based on HIRES S2N code
;-
;------------------------------------------------------------------------------
function hires_thru_newccd, wave

  dat = [ $
          [3153.90,  11.21], $
          [3183.67,   7.03], $
          [3213.44,   9.32], $
          [3243.22,   8.04], $
          [3272.99,  7.56], $
          [3302.76,  6.77], $
          [3332.53,  6.09], $
          [3362.30,  5.75], $
          [3426.54,  3.76], $
          [3462.23,  3.97], $
          [3497.92,  3.61], $
          [3533.61,  3.43], $
          [3569.30,  3.14], $
          [3604.99,  2.83], $
          [3640.68,  2.61], $
          [3676.37,  2.49], $
          [3712.06,  2.19], $
          [3831.73,  1.89], $
          [3876.29,  1.82], $
          [3920.85,  1.75], $
          [3965.41,  1.69], $
          [4054.16,  1.86], $
          [4104.22,  1.84], $
          [4154.28,  1.79], $
          [4204.34,  1.76], $
          [4254.40,  1.71], $
          [4304.46,  1.70], $
          [4354.52,  1.66], $
          [4404.58,  1.65], $
          [4454.64,  1.63], $
          [4569.19,  1.52], $
          [4635.40,  1.53], $
          [4701.61,  1.53], $
          [4767.82,  1.54], $
          [4834.03,  1.53], $
          [4900.24,  1.52], $
          [4966.45,  1.49], $
          [5032.66,  1.47], $
          [5098.87,  1.43], $
          [5165.08,  1.40], $
          [5399.52,  1.34], $
          [5491.01,  1.32], $
          [5582.50,  1.31], $
          [5673.99,  1.27], $
          [5765.48,  1.24], $
          [5856.97,  1.23], $
          [5940.00,  1.37], $
          [6046.00,  1.37], $
          [6152.00,  1.34], $
          [6258.00,  1.32], $
          [6364.00,  1.27], $
          [6479.53,  1.11], $
          [6606.56,  1.09], $
          [6733.59,  1.07], $
          [6860.62,  1.09], $
          [6987.65,  1.08], $
          [7424.34,  1.08], $
          [7589.31,  1.13], $
          [7754.28,  1.16], $
          [7919.70,  1.37], $
          [8103.85,  1.45], $
          [8288.00,  1.51], $
          [8485.40,  1.58], $
          [8697.50,  1.67], $
          [8909.60,  1.78], $
          [9378.60,  2.22], $
          [9639.10,  2.16], $
          [9899.60,  2.18], $
          [10160.10,  2.04] ]

  ccd_boost = interpol(dat[1,*], dat[0,*], wave)
  return, ccd_boost
end

function hires_blaze, wave, center, fsr

;
;
;     I/O - communicate % of blaze function peak at lambda
;
;..............................................................................
;
;
  gamma = !PI * (center - wave)/fsr
;
  blaze = replicate(1., n_elements(wave))
  nzro = where(gamma NE 0., nz)
  if nz NE 0 then blaze[nzro] = (sin(gamma[nzro])/gamma[nzro])^2
  
  return, blaze 
end 


;;;;;;;;;;;;;;;;;;;;;
function hires_thruput, wave, center, iorder, fsr, ICD=icd, $
                        BLAZE=blaze, FLG=flg

;
;  HIRES throughput estimates revised by S.Vogt, 23 July 93
;  I/O - queries for cross disperser number (1 or 2) 
;      - communicates spectrograph throughput efficiency
;..............................................................................
;
  if not keyword_set(icd) then icd = 0 
 
;	dimension xwave(NWAVE),xthru(NCD, NORDER, NWAVE)
  xwave = [3000, 3200, 3500, 3800, 4000, $
           4500, 5000, 6000, 7000, 8000, 9500]
  allthru = [ [1.e-4, 0.003],    $ ; 3000 A
            [2.e-4, 0.008], $    ; 3200 A
            [0.008, 0.020], $    ; 3500 A
            [0.020, 0.019], $    ; 3800OA A
            [0.035, 0.018], $    ; 4000 A
            [0.061, 0.009], $    ; 4500 A
            [0.077, 0.003], $    ; 5000 A
            [0.080, 3.e-4], $    ; 6000 A
            [0.068, 3.e-4], $    ; 7000 A
            [0.051, 3.e-4], $    ; 8000 A
            [0.017, 3.e-4]]       ; 9500 A
  sz = size(allthru, /dimens)

  nwv = n_elements(wave)
  thru = fltarr(nwv)
  for qq=0L,nwv-1 do begin
      xthru = reform(allthru[iorder[qq]-1,*],sz[1])

      thru1 = interpol(xthru, xwave, center[qq])
      if keyword_set(BLAZE) then $
        thru[qq] = hires_blaze(wave[qq], center[qq], fsr[qq])*thru1 $
      else thru[qq] = thru1
  endfor

  ;; New CCD?
  if keyword_set(flg) then begin
      ccd_boost = hires_thru_newccd(wave)
      thru = thru * ccd_boost
  endif
 
  return, thru
end

