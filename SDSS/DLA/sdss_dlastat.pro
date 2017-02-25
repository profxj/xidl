;+ 
; NAME:
; sdss_dlastat
;
; PURPOSE:
;    Given a DLA struct and gz list, determine the indices of those
;    DLA satisfying the statistical sample.
;
; CALLING SEQUENCE:
;   sdss_dlastat
;
; INPUTS:
;  dlastr -- Structure of SDSS DLAs
;  gzstr  -- Structure containing g(z) info
;  [vprox] -- Proximity velocity (to avoid in calculation)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /ALL   -- Return all DLA and candidates
;  BUFF=  -- Require DLA occur BUFF km/s to the red of z1
;  VMAX=  -- Require the DLA occur at less than VMAX km/s of z_em
;            (for PDLA searches)
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
function sdss_dlastat, dlastr, gzstr, vprox, BAL=bal, ALL=all, IDXA=idxa, $
                       VMAX=vmax, BUFF=buff, VMIN=vmin

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'indx = sdss_dlastat(dlastr, gzstr, [vprox], CLR=, ALL=, IDXA=), [v1.1]'
    return, -1
  endif 
  if not keyword_set(BUFF) then buff = 3000.
  if not keyword_set(VMIN) then vmin = 0.
  if not keyword_set(ALL) then dla = where(dlastr.NHI GE 20.2999, ndla) $
  else begin
      ndla = n_elements(dlastr)
      dla = lindgen(ndla)
  endelse
  msk_smpl = bytarr(ndla)  ; 0=bad, 1=good
  if arg_present(BAL) then bal = intarr(ndla)
  if arg_present(IDXA) then idxa = lonarr(ndla)

  ;; Proximates
  if keyword_set(vprox) then begin
      c = x_constants()
      vb = (vprox + buff) / (c.c / 1e5)
      ovb = 1. - vb
      vc = (vprox) / (c.c / 1e5)
      ovc = 1. - vc
      zmax = 9e99
  endif

  for ii=0L,ndla-1 do begin
      qq = dla[ii]
      ;; RA, DEC
      if strlen(strtrim(dlastr[qq].qso_ra)) EQ 0 then stop
      x_radec, dlastr[qq].qso_ra, dlastr[qq].qso_dec, rad, decd
      idx = where(abs(gzstr.ra-rad) LT 0.001 AND $
                  abs(gzstr.dec-decd) LT 0.001 AND $
                  abs(gzstr.zem-dlastr[qq].qso_zem) LT 0.03, nidx)
      case nidx of 
          0: 
          1: begin
              if keyword_set(VPROX) then begin   ;; PDLAs
                  ;; Check on buffer
                  if gzstr.z1[idx] LE x_relvel(gzstr.zem[idx],(vprox+buff)) then begin
                      ;; Velocity
                      zmin = x_relvel(gzstr.zem[idx],vprox) 
                      if keyword_set(VMAX) then zmax = x_relvel(gzstr.zem[idx],vmax) 
                      
                      if dlastr[qq].zabs GE zmin AND dlastr[qq].zabs LT zmax AND $
                        gzstr.flg_bal[idx] NE 2 then msk_smpl[ii] = 1B
                      if arg_present(BAL) then bal[qq] = gzstr.flg_bal[idx]
                      if arg_present(IDXA) then idxa[ii] = idx
                  endif
              endif else begin  ;; Intervening DLAs
                  zmin = x_relvel(gzstr.z1[idx],vmin)
                  ;;
                  if dlastr[qq].zabs GE zmin AND $
                    dlastr[qq].zabs LE gzstr.z2[idx] AND $
                    gzstr.flg_bal[idx] NE 2 then msk_smpl[ii] = 1B
              endelse
              if arg_present(BAL) then bal[qq] = gzstr.flg_bal[idx]
              if arg_present(IDXA) then idxa[ii] = idx
          end
          else: stop
      endcase
  endfor

  ;; CLR
  if keyword_set( CLR ) then begin
      noclr = where(dlastr[dla].sdss_plate LT 716)
      msk_smpl[noclr] = 0B
  endif

  ;; Return
  gd = where(msk_smpl,ngd)
  if ngd EQ 0 then return, -1 else begin
      if arg_present(IDXA) then idxa = idxa[gd]
      return, dla[gd]
  endelse

end
