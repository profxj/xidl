;+ 
; NAME:
; cldy_calctcoll
;  V1.1
;
; PURPOSE:
;    Given a pair of ions and the ratio between those ions, calculate
;   the temperature for gas with temperature T that gives that ratio
;
; CALLING SEQUENCE:
;   
;   cldy_calctcoll, ion1, ion2, rtio, tans, INFIL=, /PLOT
;
; INPUTS:
;   ion1 -- [Z,i]
;   ion2 -- [Z,i]
;   rtio -- Ratio between ion1 and ion2
;   tans -- Temperature which gives that ratio
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  OFFS - 'Error' in rtio to determine temperature range [default: 1.]
;  PLOT - Create a plot
;  COLL - CIE grid
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
; cldy_calctcoll, [7,5], [8,6], -0.9, tans, /plot
;
; PROCEDURES CALLED:
;  getabnd
;  prs_cldycoll
;
; REVISION HISTORY:
;   04-Feb-2004 Written by JXP
;-
;------------------------------------------------------------------------------
pro cldy_calctcoll, ion1, ion2, rtio, tans, INFIL=infil, PLOT=plot, O $
                    FFS=offs, COLL=coll

  if (N_params() LT 4) then begin 
    print,'Syntax - ' + $
             'cldy_calctcoll, ion1, ion2, rtio, tans, INFIL=, /PLOT, OFFS= (v1.1)'
    return
  endif 

  if not keyword_set( OFFS ) then offs = 1.
  if not keyword_set( INFIL ) then $
    infil = getenv('XIDL_DIR')+'/Cloudy/cloudy_collisions.fits'

  ;; Read in structure
  prs_cldycoll, coll, infil

  ;; Grab abundances  
  getabnd, elm, ion1[0], abd1, flag=1
  getabnd, elm, ion2[0], abd2, flag=1

  ;; Calculate model ratios

  mod_rt = coll.X[ion1[0],ion1[1]] - coll.X[ion2[0], ion2[1]] + abd1 - abd2

  tans = interpol(coll.T, mod_rt, rtio)

  if keyword_set( PLOT ) then begin
      plot, coll.T, mod_rt, xstyle=1, ystyle=1, psym=10, yrange=[rtio-offs, rtio+offs], /xlog
      oplot, [tans, tans], [-9e9, 9e9], linestyle=1
  endif


  return
end
