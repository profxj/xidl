;+ 
; NAME:
; x_tweakarc   
;     Version 1.1
;
; PURPOSE:
;  Allow the user to interactively (with a GUI) fiddle with the 1D
;  wavelength solution of a arc-line spectrum.  This is designed for
;  Echelle observations only.
;
; CALLING SEQUENCE:
;  x_tweakarc, twkfil, ordrs, templfil, _EXTRA=extra, $
;               LINLIST=linlist, QAFIL=qafil, OSTR_FIL=ostr_fil
;
; INPUTS:
;   twkfil      -- File containing the 1D arc-line fits
;   ordrs       -- Orders to tweak [physical order numbers]
;   [templfil]  -- name of archived wavelength solution to use as at
;                  template
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;  LINLIST=  -- Name of spectral line list
;  OSTR_FIL= -- Name of file containing the order structure.  Required
;              for QA
;  FLIP_TEMPL= -- Flip the template file
;
; OPTIONAL OUTPUTS:
;  QAFIL=  -- Name of file to write QA info
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES/FUNCTIONS CALLED:
;   x_fitarc
;
; REVISION HISTORY:
;   Summer-2005 Created by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro x_tweakarc, twkfil, ordrs, templfil, _EXTRA=extra, NORD=nord, $
                LINLIST=linlist, QAFIL=qafil, OSTR_FIL=ostr_fil, $
                FLIP_TEMPL=flip_templ


;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'x_tweakarc, twkfil, ordrs, [templfil], QAFIL=, LINLIST=, /FLIP_TEMPL [v1.0]'
      return
  endif 
  
;;  Optional Keywords
  if not keyword_set( LINLIST ) then $
    linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 

  ;; Open line list
  x_arclist, linlist, lines


  ;; Open templfil
  if keyword_set( TEMPLFIL ) then begin
      restore, templfil
      
      if not keyword_set(GUESS_ORDR) then begin
          guess_ordr = extra.guess_ordr
      endif
      
      tguess = guess_ordr
      tfit = all_arcfit
      tspec = sv_aspec
      if keyword_set(FLIP_TEMPL) then begin 
         ;; MIKE blue
         if NOT keyword_set(CCDSZ) then ccdsz=[2048L,4096]
         ccol = (size(sv_aspec))[1]
         cbin = round(float(ccdsz[1])/ccol)
         tspec = reverse(sv_aspec)
         for i=0L, n_elements(all_arcfit)-1 do begin
            if ptr_valid(all_arcfit[i].ffit) then begin
               ffit = *(all_arcfit[i].ffit)
;               all_arcfit[i].nrm[0] = ccdsz[1] * all_arcfit[i].nrm[0]/ncol/cbin
;               all_arcfit[i].nrm[1] = ccdsz[1] * all_arcfit[i].nrm[1]/ncol/cbin
               ;; Flip
               ffit = ffit * (1.0 -2.0 * (lindgen(20) mod 2))
               *(all_arcfit[i].ffit) = ffit
               all_arcfit[i].nrm[0] = (ccdsz[1]-1)/cbin $
                                      - all_arcfit[i].nrm[0]
            endif 
         endfor
         tfit = all_arcfit
      endif
  endif

  ;; Open tweak file
  restore, twkfil

  nordr = n_elements(ordrs)
  for qq=0L,nordr-1 do begin
      if keyword_set(TEMPLFIL) then at = where(tguess EQ ordrs[qq],nat) $
        else nat = 0
      ak = where(guess_ordr EQ ordrs[qq], nak)
      if nak + nat EQ 2 then begin
          print, 'x_tweakarc: Tweaking order ', ordrs[qq]
          ;; FIT by hand
          x_identify, sv_aspec[*,ak], calib, MSPEC=tspec[*,at], $
          MFITSTR=tfit[at], LINELIST=linlist, OUTLIN=outlin, _EXTRA=extra
          ;; Save
          ngd = n_elements(outlin)
          sv_lines[ak].nlin = ngd
          sv_lines[ak].pix[0:ngd-1] = outlin.pix
          sv_lines[ak].wv[0:ngd-1] = outlin.wave
          ;; Cant deal with additional rejected lines for now  JXP
          if all_arcfit[ak].nord NE 0 then tmpfit = all_arcfit[ak] else $
            tmpfit = tfit[at]
          if keyword_set(NORD) then tmpfit.nord = NORD $
          else tmpfit.nord = (calib.nord > 2)
          print, 'x_tweakarc: Using NORD = ', tmpfit.nord
          fit = x_fitrej(outlin.pix, alog10(outlin.wave), FITSTR=tmpfit)
          ;; Save
          all_arcfit[ak] = tmpfit
      endif else begin
          if nak EQ 1 then begin
              print, 'x_tweakarc: Checking order ', ordrs[qq]
              lines.flg_plt = 0
              for kk=0L,sv_lines[ak].nlin-1 do begin
                  a = where(abs(sv_lines[ak].wv[kk] - lines.wave) LT 0.01, na)
                  if na NE 0 then begin
                      lines[a].flg_plt = 1
                      lines[a].pix = sv_lines[ak].pix[kk]
                  endif
              endfor
              svlines = lines
              x_identify, sv_aspec[*,ak], calib, LINELIST=linlist, $
                OUTLIN=outlin, INLIN=lines, _EXTRA=extra
              rejpt = where(lines.flg_plt EQ 0 AND svlines.flg_plt EQ 1)
              ;; Save gdfit (for QA)
              gdfit = where(lines.flg_plt EQ 1, ngd)
              rejstr[ak].ngdf = ngd
              if ngd NE 0 then begin
                  rejstr[ak].gdfpt[0:ngd-1] = gdfit
                  rejstr[ak].gdfpx[0:ngd-1] = lines[gdfit].pix
              endif
              ;; Save
              ngd = n_elements(outlin)
              sv_lines[ak].nlin = ngd
              sv_lines[ak].pix[0:ngd-1] = outlin.pix
              sv_lines[ak].wv[0:ngd-1] = outlin.wave

              ;; Rejected points (for QA)
              if rejpt[0] NE -1 then begin
                  rejstr[ak].nrej = 0
                  rejstr[ak].nrej = n_elements(rejpt)
                  rejstr[ak].rejpt[0:rejstr[ak].nrej-1] = rejpt
                  rejstr[ak].rejwv[0:rejstr[ak].nrej-1] = $
                    lines[gdfit[rejpt]].wave
              endif

              ;; Cant deal with additional rejected lines for now  JXP
              calib.flg_rej = 1
              calib.lsig = 2.
              calib.hsig = 2.
              calib.maxrej = 10L
              fit = x_fitrej(outlin.pix, alog10(outlin.wave), FITSTR=calib)
              all_arcfit[ak] = temporary(calib)
          endif
      endelse
  endfor

  save, guess_ordr, sv_lines, sv_aspec, all_arcfit, rejstr, filename=twkfil
  ;; 
  resolve_routine, 'x_fitarc', /no_recompile, /either, /COMPILE_FULL_FILE
  if keyword_set( QAFIL ) then begin
      if not keyword_set(OSTR_FIL) then return
      ordr_str = xmrdfits(ostr_fil, 1, /silent)
      nordr = n_elements(ordr_str)
      x_psopen, qafil, /maxs
      clr = getcolor(/load)
      !p.multi=[0,3,2]
      
      ;; Plot to ps file  (this does not include edited lines)
      sz_arc = size(sv_aspec,/dimensions)
      for ii=0L,nordr-1 do begin
          if rejstr[ii].ngdf EQ 0L then continue
          gdfit = rejstr[ii].gdfpt[0:rejstr[ii].ngdf-1]
          ;; Rejpt
          if rejstr[ii].nrej NE 0 then rejpt = $
            rejstr[ii].rejpt[0:rejstr[ii].nrej-1] $
          else rejpt = -1L
          ;; Dwv
          wv = 10^x_calcfit(dindgen(sz_arc[1]),fitstr=all_arcfit[ii])
          if not keyword_set(BCKWD) then dwv = median(wv - shift(wv,1)) else $
            dwv = median(wv - shift(wv,-1))
          ;; Subroutine 
          fit=10^x_calcfit(double(rejstr[ii].gdfpx[0:rejstr[ii].ngdf-1]), $
                           FITSTR=all_arcfit[ii])
          x_fitarc_ps, WV=lines[gdfit].wave, $
            fit=fit, $
            REJ=rejpt, $
            DWV=dwv, $
            ORDR=ordr_str[ii].order,$
            RMS=all_arcfit[ii].rms*sv_lines[ii].wv[0]*alog(10.) / dwv, $ ; pix
            FORDR=all_arcfit[ii].nord
      endfor
      
      x_psclose
      !p.multi=[0,1,1]
      spawn, 'gzip -f '+qafil
  endif

; All done
  print, 'x_fitarc: All done!'

  return
end
