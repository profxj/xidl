;+ 
; NAME: 
; x_lyman_ew_spline   
;    Version 1.1
;
; PURPOSE:
;    Generate a spline solution for the EW of Lyman series lines for a
;    given Doppler parameter
;
; CALLING SEQUENCE:
;   
;
; INPUTS:
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
; x_lyman_ew_spline, 24.
;
; PROCEDURES/FUNCTIONS CALLED:
;  
;
; REVISION HISTORY:
;   May-2011 Written by JXP
;-
;------------------------------------------------------------------------------
pro x_lyman_ew_spline, bval, EW_FIL=ew_fil

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'x_lyman_ew_spline, bval, EW_FIL=, [v1.0]'
    	return
  endif 

  if not keyword_set(EW_FIL) then begin
     ew_fil = 'EW_SPLINE_b'+strtrim(round(bval),2)+'.fits'
  endif
  
  c = x_constants()

  lseries= [1215.6701, 1025.7223, 972.5368, $
            923.1504d, 926.2257, 930.7483, 937.8035, 949.7431, $
            920.9631, 919.3514, 918.1294, 917.1806, 916.429, 915.824, $
            915.329, 914.919, 914.576, 914.286, 914.039, $
            913.826, 913.641, 913.480, 913.339, 913.215, 913.104, $
            913.006, 912.918, 912.839, 912.768, 912.703, 912.645  ]
  srt = reverse(sort(lseries))
  lseries = lseries[srt]
  nls = n_elements(lseries)
  print, 'N lines = ', nls
  for jj=0L,nls-1 do begin
     if jj EQ 0 then lines=x_setline(lseries[0]) else $
        lines = [lines, x_setline(lseries[jj])]
  endfor

  ;; Write
;  printcol, lindgen(nls), lines.ion
  ;mwrfits, lines, 'Lyman_lines.fits', /create
  ;save, lines, filen='Lyman_lines.idl'

  ;; EW SPLINE
  nspl = 100L
  NHI = 11.0 + 11*findgen(nspl)/(nspl-1)
  tmp = { $
        wrest: 0.d, $
        NHI: NHI, $
        EW: fltarr(nspl), $
        splint: fltarr(nspl) $
        }


  EW_SPLINE = replicate(tmp, nls)
  EW_SPLINE.wrest = lines.wrest

  nvel = 60001L
  velo = -30000.d + findgen(nvel) ; km/s
  dvel = 1. ; km/s
  uval = velo / bval

  ;; Loop
  for qq=0L,nls-1 do begin

     ;; Wave array
     dwv = dvel * lines[qq].wrest / (c.c/1e5)  ;; Ang

     ;; Voigt
     vd = bval/ (lines[qq].wrest * 1.0d-13)  ;; Frequency
     a = lines[qq].gamma / (12.56637 * vd)
     vgt = voigt(a, uval)

     ;; tau
     tau = 0.014971475*lines[qq].f*vgt/vd  ;; Normalized to N_HI = 1 cm^-2

     ;; Flux
     tau_array = tau # (10.d^NHI)
     fx = exp(-1.*tau_array)

     ;; EW
     EW = total(1.-fx, 1) * dwv
     EW_SPLINE[qq].EW = EW
     
     ;; Spline
     EW_SPLINE[qq].splint = spl_init(NHI, EW, /double)
  endfor

  ;; Write
  mwrfits, EW_SPLINE, EW_FIL, /create

  return
end

     
