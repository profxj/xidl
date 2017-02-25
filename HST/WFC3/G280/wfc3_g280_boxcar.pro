;+ 
; NAME:
;  wfc3_g280_boxcar
;
; PURPOSE:
;   Simple algorithm to perform boxcar extraction on the WFC3/G280 data
;
; CALLING SEQUENCE:
;   
;  wfc3_g280_boxcar, wfc3_g280_strct, ii, specim, var
;
; INPUTS:
;   wfc3_g280_strct -- the wfc3_g280 structure
;   ii -- the index of the object in the structure
;   specim -- Spectral image (should be sky subtracted)
;   var    -- Variance image
;
; RETURNS:
;
; OUTPUTS:
;  Updates the wfc3_g280 structure to include the boxcar flux.
;
; OPTIONAL KEYWORDS:
;  BOX_SIZE= -- Total width in pixels of the boxcar [default: 9L]
;  BEAM= -- Beam to reduce
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
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Major Update (structure + multibeam) by MN
;------------------------------------------------------------------------------

pro wfc3_g280_boxcar, wfc3_g280_strct, ii, specim, var, $
                      BOX_SIZE=box_size, BEAM=beam

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
          'wfc3_g280_boxcar, specim, var, trace, BOX_SIZE= [v1.0]'
    return
  endif 
  
  if not keyword_set(BOX_SIZE) then box_size = 9L

  case beam of
     0: begin
        cnt=wfc3_g280_strct(ii).cnta
        trace_x=wfc3_g280_strct(ii).trace_xa
        trace_y_fit=wfc3_g280_strct(ii).trace_ya_fit
     end
     1: begin
        cnt=wfc3_g280_strct(ii).cntc
        trace_x=wfc3_g280_strct(ii).trace_xc
        trace_y_fit=wfc3_g280_strct(ii).trace_yc_fit
     end
     else: stop
  endcase
  
  ;; Simple addition
  dim=size(specim,/dim)
  yoff = (box_size-1)/2
  count = dblarr(cnt)
  tvar = dblarr(cnt)

  for kk=0L,cnt-1 do begin
     xval = trace_x[kk]
     yval = round(trace_y_fit[kk])
     if xval ge dim(0) or xval lt 0 or $
        yval-yoff+BOX_SIZE-1 ge dim(1) or $
        yval-yoff lt 0 then continue
     count[kk] = total(specim[xval, yval-yoff+lindgen(BOX_SIZE)])
     tvar[kk] = total(var[xval, yval-yoff+lindgen(BOX_SIZE)])
  endfor

  ;; Save
  case beam of
     0: begin
        wfc3_g280_strct(ii).countsa=count
        wfc3_g280_strct(ii).vara=tvar
     end
     1: begin
        wfc3_g280_strct(ii).countsc=count
        wfc3_g280_strct(ii).varc=tvar
     end
     else: stop
  endcase

end
