;+ 
; NAME:
; dla_abund_vpfit
;  V1.2
;
; PURPOSE:
;    Parse a VPFIT output file and fill up relevant things.
;
; CALLING SEQUENCE:
;   dla_abund_vpfit, dla, nn
;
; INPUTS:
;   dla -- DLA structure
;   nn  -- Index of the structure
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   12-Jul-2012 Written by JXP
;- 
;------------------------------------------------------------------------------
pro dla_abund_vpfit, dla, nn

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'dla_calcmtl, dla, nn [v1.1]'
    return
  endif 

  ;; Parse
  x_parse_vpfit, strtrim(dla[nn].vpfit_fil,2), vp_strct

  ;; Sum up CI
  iCI = where( strmatch(strmid(vp_strct.ion,0,2), 'CI') AND $
               (strmatch(strmid(vp_strct.ion,3,1), '') OR $
                strmatch(strmid(vp_strct.ion,3,1), '*')), nCI)
  if nCI GT 0 then begin
     x_logclm, linear_N, linear_sig, vp_strct[iCI].N, vp_strct[iCI].Nsig, /REVER
     tot_N = total(linear_N)
     tot_sig = sqrt( total(linear_sig^2))
     ;;
     x_logclm, tot_N, tot_sig, log_CI, siglog_CI
     dla[nn].flg_CI = 1
     dla[nn].CI = log_CI
     dla[nn].sig_CI = siglog_CI
  endif

  ;; Sum up H2 
  iH2 = where( strmatch(strmid(vp_strct.ion,0,2), 'H2'), nH2)
  if nH2 GT 0 then begin
     x_logclm, linear_N, linear_sig, vp_strct[iH2].N, vp_strct[iH2].Nsig, /REVER
     tot_N = total(linear_N)
     tot_sig = sqrt( total(linear_sig^2))
     ;;
     x_logclm, tot_N, tot_sig, log_H2, siglog_H2
     dla[nn].flg_H2 = 1
     dla[nn].H2 = log_H2
     dla[nn].sig_H2 = siglog_H2
  endif
     
  return
end
