;+ 
; NAME:
; esi_echfitarc   
;     Version 1.1
;
; PURPOSE:
;    1) Take the arc image from esi_echmkarc 
;    2) Extract 1D spectra (1 per order)
;    3) Determine pixel offset from archived fit
;    4) Automatically identify a set of lines
;    5) Create a fit to these lines, identify all lines and fit again
;    6) Write arc solutions (one per order) to a fits file 
;
; CALLING SEQUENCE:
;   
;  esi_echfitarc, esi, slit, /INTER, LINLIST=, /CHK, /CLOBBER,
;  SIGREJ=, /DEBUG, ORDRS=, PIXSHFT=, /PINTER
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  IDL fit file (one per order)  (e.g. Arcs/ArcECH_##fit.idl)
;
; OPTIONAL KEYWORDS:
;   /PINTER   - Perform fit only for the pre-identified lines
;              (recommended)
;   /INTER    - Identify lines interactively and then fit
;   LINLIST=  -  Arc line list (default:
;                                $XIDL_DIR/Spec/Arcs/Lists/ESI_ech.lst)
;   GUESSARC= -  Initial guess for wavelength solution (default:
;                                $XIDL_DIR/ESI/CALIBS/ESI_arcfit.idl)
;   /CHK      - Manually check steps along the way
;   /DEBUG    - More intensive checking
;   ORDRS=    - Only fit select orders (default: [0L,9L] )
;   PIXSHFT=  - Manually set pixel shift from calib file (deafult:
;                Let the program determine this using FFT formalism)
;   SIGREJ=   - Rejection sigma for outliers in arc line fitting
;              (default: 2.)
;   /CLOBBER  - Overwrite previous fits
;   /CUAR     - CuAr lamps only!
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   There is some fine tuning in a couple of orders to help with the
;    fits.  Only setup for 1x1 binning right now
;
; EXAMPLES:
;   esi_echfitarc, esi, 0.5, /CHK, /PINTER, /CLOBBER
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   21-Aug-2002 Streamlined + Added ps output
;   01-Feb-2003 Polished (JXP)
;   16-Sep-2004 Modified to fit alog10 wavelengths by JXP
;-
;------------------------------------------------------------------------------

pro esi_echfitarc_ps, flg, svdecomp, filnm, WV=wv, FIT=fit, REJPT=rejpt,$
                      GSWV=gswv, ORDR=ordr, RMS=rms

  COMMON ps_common, old_dname, old_pfont, file_name, opened, color, ledger

  ;; FLG
  case flg of
      0: begin  ; Open ps file
          device, get_decomposed=svdecomp
          device, decompose=0
          ;; Open
          ps_open, file=filnm, font=1, /color
      end
      1: begin
          clr = getcolor(/load)
          ;; All points
          plot, wv, fit-wv, psym=1, $
            charsize=1.8, $
            background=clr.white, color=clr.black, $
            xtitle='Wave', ytitle='Residual (Ang)', $
            xmargin=[12,2], ymargin=[6,2], yrange=[-0.2, 0.2]
          ;; MASK 
          if rejpt[0] NE -1 then $
            oplot, wv[rejpt], fit[rejpt]-wv[rejpt], psym=2, color=clr.red
          ;; Order
          xyouts, 0.2, 0.9, 'Order = '+string(15-ordr, FORMAT='(i2)'), /normal, $
            charsize=2.0
          xyouts, 0.2, 0.83, 'RMS = '+string(rms, FORMAT='(f5.3)'), /normal, $
            charsize=2.0
      end
      9: begin ; CLOSE
          if keyword_set(opened) then begin
              if opened EQ 1 then ps_close, /noprint, /noid
          endif
          ;; Device
          if keyword_set( SVDECOMP ) then device, decomposed=svdecomp
      end
      else:
  endcase
          
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echfitarc, esi, slit, INTER=inter, LINLIST=linlist, $
                   CHK=chk, CLOBBER=clobber, SIGREJ=sigrej, DEBUG=debug,$
                   ORDRS=ordrs, PIXSHFT=pixshft, PINTER=pinter, CUAR=cuar,$
                   GUESSARC=guessarc, NORD=nord, CBIN=cbin, RBIN=rbin, $
                   PKSIG=pksig, MAXORD=maxord, GLOG=glog, $
                   SEDG_FIL=sedg_fil


;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echfitarc, esi, slit, /INTER, ORDRS=, /DEBUG, /CHK '
      print, '     PIXSHIFT=, /CUAR, CBIN=, RBIN=, /GLOG [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( ORDRS ) then begin
      ordrs = [0L, 9L] 
      flg_ordr = 0
  endif else flg_ordr = 1

  if not keyword_set(MXSHFT) then mxshft = 25
  if not keyword_set( MAXORD ) then maxord = 4
  if not keyword_set( PKSIG ) then pksig = 7.
  if not keyword_set( SIGREJ ) then sigrej = 2.
  if not keyword_set( TRCNSIG ) then trcnsig = 3.
  if not keyword_set(MEDWID) then medwid = 5L
  if not keyword_set( RADIUS ) then radius = 2.0
  if not keyword_set(LINLIST) then begin
      if keyword_set( CUAR ) then $
        linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ESI_CuArech.lst' $
      else linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ESI_ech.lst'
  endif
  if not keyword_set( GUESSARC ) then begin
      if keyword_set( CUAR ) then $
        guessarc = getenv('ESI_CALIBS')+'/ECH_CuArarcfit.idl' $
      else guessarc = getenv('ESI_CALIBS')+'/ECH_arcfit.idl' 
  endif
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1

; Slit name
;  c_s = esi_slitnm(slit)

; Set Peak Width
  case slit of
      0.3: pkwdth = 4L
      0.5: pkwdth = 5L
      0.75: pkwdth = 6L
      1.0: pkwdth = 7L
      else: stop
  endcase

; Check for outfil
  outfil = esi_getfil('arc_fit', SLIT=slit, cbin=cbin, rbin=rbin, /name)
  if flg_ordr EQ 0 then begin
      a = findfile(outfil, count=na)
      if na NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'esi_echfitarc: Arc fits file exists!  Use /CLOBBER to overwrite'
          return
      endif
  endif 

; Open line list
  x_arclist, linlist, lines

; Open ps file
  ;; CLOSE first
  esi_echfitarc_ps, 9, svdevice, psfil
  if not keyword_set( CHK ) AND not keyword_set( PINTER ) then begin
      psfil = esi_getfil('arc_psfil', SLIT=slit, cbin=cbin, rbin=rbin, /name)
      esi_echfitarc_ps, 0, svdevice, psfil
  endif
  
; Grab Arc IMG

  arc_fil = esi_getfil('arc_fil', SLIT=slit, cbin=cbin, rbin=rbin, /name)
  a = findfile(arc_fil+'*', count=na)
  if na EQ 0 then begin
      print, 'esi_echfitarc: Arc ', arc_fil, ' does not exist. Run esi_echmkarc!'
      return
  endif
  print, 'esi_echfitarc: Reading arc: ', arc_fil
  arc_img = xmrdfits(arc_fil, 0, head, /silent)
  sz_arc = size(arc_img, /dimensions)

; Open Slit file
  if not keyword_set(SEDG_FIL) then $
    slit_edg = esi_getfil('sedg_fil', SLIT=slit, cbin=cbin, rbin=rbin) $
  else slit_edg = xmrdfits(sedg_fil, 0, /silent)
  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)

  rnd_edg = round(slit_edg)

; Open Guess (If not interactive)
  if not keyword_set( INTER ) then begin
      if x_chkfil(guessarc+'*') EQ 0 then begin
          print, 'esi_echfitarc:  Put ', guessarc, 'in the CALIBS file!!'
          stop
      endif
      restore, guessarc
      guess_spec = temporary(sv_aspec)
      guess_fit = temporary(all_arcfit)
      x_arclist, linlist, lines, /GDONLY
      ;; SORT
      srt = sort(lines.wave)
      lines = lines[srt]
      guess_fit.hsig = 3.
      guess_fit.lsig = 3.
      guess_fit.flg_rej = 1
      guess_fit.niter = 2
  endif

; Open Old file as necessary
  if flg_ordr EQ 1 then restore, outfil

; Loop on Orders

  msk = bytarr(sz_arc[1])
  tmp = fltarr(2*MEDWID+1L, sz_arc[1])
  ;; Create output as necessary
  if flg_ordr EQ 0 then begin
      sv_aspec = fltarr(sz_arc[1], 10)
      fittmp = { fitstrct }
      all_arcfit = replicate(fittmp, 10)
      lintmp = { $
                 pix: dblarr(90), $
                 wv: dblarr(90), $
                 nlin: 0 $
               }
      sv_lines = replicate(lintmp, 50)
  endif
      
; LOOP
  for qq=ordrs[0], ordrs[1] do begin 
      print, 'esi_echfitarc: Fitting order ', string(15L-qq, FORMAT='(i3)')

      ;; Zero out
      lines[*].flg_plt = 0
      tmp[*]=0.
      ;; Extract a 1D spectrum
      gd = where(slit_cen[*,qq] - MEDWID GE 0L, ngd)
      for q=0L,ngd-1 do begin
          j = gd[q]
          tmp[*,j] = arc_img[slit_cen[j,qq]-MEDWID:slit_cen[j,qq]+MEDWID,j]
      endfor
      sv_aspec[*, qq] = djs_median(tmp > (-10.0), 1) 

      ;; CuAr
      if keyword_set( CUAR ) AND qq EQ 9L then begin
          all_arcfit[qq] = guess_fit[qq]
          ;; Output
          print, 'esi_echfitarc: Writing fit to ', outfil, ' and returning'
          save, sv_aspec, all_arcfit, sv_lines, filename=outfil
          break
      endif
          
          
      if keyword_set(INTER) then begin ;; INTERACTIVE
          x_identify, sv_aspec[*,qq], finfit, LINELIST=linlist
          stop
      endif else begin  ;; FIT
          ;; MASK
          msk = bytarr(sz_arc[1])
          case qq of 
              0: msk[1400:3790] = 1B
              1: begin
                  msk[0:200] = 1B
                  msk[800:sz_arc[1]-1] = 1B
              end
              2: begin
                  msk[0:200] = 1B
                  msk[700:3900] = 1B
              end
              3: begin
                  msk[0:200] = 1B
                  msk[500:sz_arc[1]-1] = 1B
              end
              4: msk[400:sz_arc[1]-1] = 1B
              else: msk[0:sz_arc[1]-1] = 1B
          endcase

          ;; ADD SOME LINES BY HAND!
          case qq of 
              3: begin  ;; ORDER 12
                  ;; 4726
                  a = where(long(lines.wave) EQ 4726)
                  lines[a].pix = 210.4
                  lines[a].flg_plt = 1
              end
              else:
          endcase

          ;; Get template
          if keyword_set( DEBUG ) then stop
          if not keyword_set(PIXSHFT) then begin
              x_templarc, sv_aspec[*,qq], lines, guess_fit[qq], /FFT, MSK=msk, $
                MOCK_FFT=fft(guess_spec[*,qq]), SHFT=shft, ALL_PK=all_pk, $
                PKWDTH=pkwdth, /THIN, FORDR=9, LOG=GLOG, MXSHFT=MXSHFT
          endif else begin
              shft = pixshft[0]
              x_templarc, sv_aspec[*,qq], lines, guess_fit[qq], MSK=msk, $
                SHFT=shft, ALL_PK=all_pk, PKWDTH=pkwdth, /THIN, FORDR=9, LOG=GLOG
          endelse
          print, 'esi_echfitarc: Shifted ', shft, ' pixels'

          ;; Good lines
          gdfit = where(lines.flg_plt EQ 1, ngd)  
        if ngd LE 5 then begin
;          if ngd LT 5 then begin
              print, 'esi_echfitarc: Insufficient lines for AUTO!!'
              stop
          endif

          ;; FIRST FIT
          tmp_fit = {fitstrct}
          copy_struct, guess_fit[qq], tmp_fit, EXCEPT_TAGS=['FFIT','NRM','RMS']
          tmp_fit.flg_rej = 1 
          tmp_fit.niter = 2 
          tmp_fit.maxrej = 5 
          tmp_fit.minpt = 5 
          tmp_fit.hsig = sigrej
          tmp_fit.lsig = sigrej
          if keyword_set(MAXORD) then tmp_fit.nord = tmp_fit.nord < MAXORD

          ;;
          fin_fit = tmp_fit
          if not keyword_set( NORD ) then begin
             fin_fit.nord = 5L
             ;; added by JFH Sep 07, 2015
             ;; there was higher order structure in the fits for 0.75 slit, and
             ;; increasing the order reduced the RMS
             ;; Not sure why we were only doing this for the
             ;; 1.0" slit, as I suspect it helps the fits for all slits
             ;;case slit of
             ;;   1.0: fin_fit.nord = 5L
             ;;   else:
             ;;endcase
          endif else begin
              fin_fit.nord = nord[qq]
          endelse

          ;; TMP fit 
          fit = x_fitrej(lines[gdfit].pix, alog10(lines[gdfit].wave), $
                         FITSTR=tmp_fit, REJPT=rejpt)
          if fit[0] EQ -1 then begin
              print, 'esi_echfitarc: AUTO Failed!!'
              stop
          endif
;          print, 'esi_echfitarc: RMS = ', tmp_fit.rms
          
          ;; Grab new lines
          lines.flg_plt = 0
          x_templarc, sv_aspec[*,qq], lines, tmp_fit, MSK=msk, PKWDTH=pkwdth, $
            /THIN, FORDR=9, /LOG, PKSIG=pksig
          
          ;; ADD SOME LINES BY HAND!
          case qq of 
              3: begin  ;; ORDER 12
                  ;; 4726
                  a = where(long(lines.wave) EQ 4726)
                  lines[a].pix = 210.4
                  lines[a].flg_plt = 1
              end
              else:
          endcase

          gdfit = where(lines.flg_plt EQ 1, ngd)

          ;; CHK 
          if keyword_set( CHK ) then begin
              if keyword_set( DEBUG ) then $
                printcol, lines[gdfit].pix, lines[gdfit].wave
              x_prspeaks, sv_aspec[*,qq], lines[gdfit].pix, /block
              if keyword_set( DEBUG ) then stop
          endif

          ;; FIT
          if not keyword_set( PINTER ) then begin ;; AUTO FIT
              fit = x_fitrej(lines[gdfit].pix, alog10(lines[gdfit].wave), $
                             FITSTR=fin_fit, REJPT=rejpt, GDPT=gdpt)
              if fit[0] EQ -1 then begin
                  print, 'esi_echfitarc: Lowering order by 1'
                  fin_fit.nord = fin_fit.nord - 1
                  fit = x_fitrej(lines[gdfit].pix, alog10(lines[gdfit].wave), $
                                 FITSTR=fin_fit, REJPT=rejpt, GDPT=gdpt)
                  if fit[0] EQ -1 then stop
              endif
              ;; Save lines
              ngd = n_elements(gdpt)
              sv_lines[qq].nlin = ngd
              sv_lines[qq].pix[0:ngd-1] = lines[gdfit[gdpt]].pix
              sv_lines[qq].wv[0:ngd-1] = lines[gdfit[gdpt]].wave
          endif else begin ;; INTER FIT in x_identify
              x_identify, sv_aspec[*,qq], id_fit, LINELIST=linlist, $
                INLIN=lines, OUTLIN=outlin
              ngd = n_elements(outlin)
              sv_lines[qq].nlin = ngd
              sv_lines[qq].pix[0:ngd-1] = outlin.pix
              sv_lines[qq].wv[0:ngd-1] = outlin.wave
              ;; Final auto fitting in log10 space!
              fin_fit.nord = id_fit.nord
              fit = x_fitrej(outlin.pix, alog10(outlin.wave), FITSTR=fin_fit)
          endelse

          wv = 10^(x_calcfit(findgen(sz_arc[1]),FITSTR=fin_fit))
          dwv = median(wv-shift(wv,1))
          print, 'esi_echfitarc: RMS(pix) = ', $
            fin_fit.rms*sv_lines[qq].wv[0]*alog(10.)/dwv, ' for ', $
            strtrim(n_elements(gdfit),2), ' good lines'
          all_arcfit[qq] = temporary(fin_fit)

          ;; Plot to ps file  (this does not include edited lines)
          if not keyword_set( CHK ) AND not keyword_set( PINTER ) then $
            esi_echfitarc_ps, 1, WV=lines[gdfit].wave, FIT=10^fit, REJ=rejpt, $
            gswv=x_calcfit(lines[gdfit].pix, FITSTR=guess_fit[qq]), ORDR=qq, $
            RMS=all_arcfit[qq].rms*sv_lines[qq].wv[0]*alog(10.)
          
      endelse
      ;; Output
      print, 'esi_echfitarc: Writing fit to ', outfil
      save, sv_aspec, all_arcfit, sv_lines, filename=outfil
  endfor


; Close ps file
  if not keyword_set( CHK ) then begin
      esi_echfitarc_ps, 9, svdevice
  endif

  return
end
