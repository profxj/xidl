;+ 
; NAME:
; slls_stat
;
; PURPOSE:
;    Given a DLA struct and gz list, determine the indices of those
;    DLA satisfying the statistical sample.
;
; CALLING SEQUENCE:
;   slls_stat
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
;   Nov-2004 Written by JXP
;-
function slls_stat, llsstr, gzstr, BAL=bal, ALL=all, IDXA=idxa

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'indx = slls_stat(llsstr, gzstr, ALL=, IDXA=), [v1.1]'
    return, -1
  endif 
  if not keyword_set(ALL) then lls = where(llsstr.NHI LT 20.3, nlls) $
  else begin
      nlls = n_elements(llsstr)
      lls = lindgen(nlls)
  endelse
  msk_smpl = bytarr(nlls)  ; 0=bad, 1=good
  if arg_present(BAL) then bal = intarr(nlls)
  if arg_present(IDXA) then idxa = lonarr(nlls)

  for ii=0L,nlls-1 do begin
      qq = lls[ii]
      ;; RA, DEC
      x_radec, llsstr[qq].qso_ra, llsstr[qq].qso_dec, rad, decd
      idx = where(abs(gzstr.ra-rad) LT 0.001 AND $
                  abs(gzstr.dec-decd) LT 0.001 AND $
                  abs(gzstr.zem-llsstr[qq].qso_zem) LT 0.03, nidx)
      case nidx of 
          0: 
          1: begin
              ;; z windows
              gd = where((llsstr[qq].zabs-gzstr.z1[*,idx])*$
                (llsstr[qq].zabs-gzstr.z2[*,idx]) LT 0., ngd)
              if ngd GT 1 then stop
              ;; NHI cut dependent on ESI vs ECHELLE
              flg_s = 0
              tmps=gzstr.flg_smpl[idx]
              case tmps[0] of
                  1: if llsstr[qq].NHI GE 19.3 then flg_s = 1
                  2: if llsstr[qq].NHI GE 19. then flg_s = 1
                  else: stop
              endcase
              if NGD NE 0 AND flg_s then msk_smpl[qq] = 1B
              ;; Other
              if arg_present(BAL) then bal[qq] = gzstr.flg_bal[idx]
              if arg_present(IDXA) then idxa[ii] = idx
          end
          else: stop
      endcase
  endfor

  ;; Return
  gd = where(msk_smpl,ngd)
  if ngd EQ 0 then return, -1 else begin
      if arg_present(IDXA) then idxa = idxa[gd]
      return, lls[gd]
  endelse

end
