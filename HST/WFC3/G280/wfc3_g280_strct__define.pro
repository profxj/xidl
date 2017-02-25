;+ 
; NAME:
;  wfc3_g280_strct__define
;
; PURPOSE:
;   Defines the structure for the reduction of WFC3_G280 data
;
; CALLING SEQUENCE:
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
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Jun-2016 Written by MN
;------------------------------------------------------------------------------

pro wfc3_g280_strct__define, wfc3_struct

  tmp={wfc3_g280_strct, $
       name:'',$
       ra:0D,$
       dec:0D,$
       exptime:0D,$
       img_fil:'',$
       spec_fil:'',$
       xguess:0D,$
       yguess:0D,$
       x0:0D,$
       y0:0D,$
       chip:0L,$
       trace_xa:lonarr(1000),$
       trace_ya_orig:dblarr(1000),$
       trace_ya_fit:dblarr(1000),$
       trace_offseta:dblarr(1000),$
       trace_sigmaa:dblarr(1000),$
       trace_sigma_fita:dblarr(1000),$
       median_trace_offseta:0D,$
       wavea:dblarr(1000),$
       cnta:0L,$
       trace_xc:lonarr(1000),$
       trace_yc_orig:dblarr(1000),$
       trace_yc_fit:dblarr(1000),$
       trace_offsetc:dblarr(1000),$
       trace_sigmac:dblarr(1000),$
       trace_sigma_fitc:dblarr(1000),$
       median_trace_offsetc:0D,$
       wavec:dblarr(1000),$
       cntc:0L,$ 
       countsa:dblarr(1000),$
       vara:dblarr(1000),$
       skya:dblarr(1000),$
       countsc:dblarr(1000),$
       varc:dblarr(1000),$
       skyc:dblarr(1000),$
       fluxa:dblarr(1000),$
       flux_siga:dblarr(1000),$
       fluxc:dblarr(1000),$
       flux_sigc:dblarr(1000),$
       badpixa:dblarr(1000),$
       badpixc:dblarr(1000),$
       cnt:0L,$ 
       wave:dblarr(1000),$
       flux:dblarr(1000),$
       sig:dblarr(1000),$
       sky:dblarr(1000),$
       extra:''$
      }
  
end
