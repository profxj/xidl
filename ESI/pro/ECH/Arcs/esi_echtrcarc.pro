;+ 
; NAME:
; esi_echtrcarc   
;     Version 1.1
;
; PURPOSE:
;    Takes the arc, creates a fit and traces the lines.
;
; CALLING SEQUENCE:
;   
;  esi_echtrcarc, esi, slit, LINLIST=, /CHK, /AUTO, ORDR=
;
; INPUTS:
;   esi     -  ESI structure
;   slit    -  Slit width  (0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   /AUTO     - Perform steps automatically (recommended)
;   ORDR=     - Only fit select order [one at a time]
;   /CHK      - Manually check steps along the way
;   /DEBUG    - More intensive checking
;   LINLIST=  - Arc line list (default:
;                                $XIDL_DIR/Spec/Arcs/Lists/ESI_ech.lst)
;   GUESSARC= -  Initial guess for wavelength solution (default:
;                                $XIDL_DIR/ESI/CALIBS/ESI_arcfit.idl)
;   /NO11KLDG - Turn off the kludge in order 11.  Only do this if the
;       arc image is very good (long CuAr exposures)
;   PIXSHFT=  - Manually set pixel shift from calib file (deafult:
;                Let the program determine this using FFT formalism)
;   PKSIG=    - Number of sigma for peaks to use in tracing 
;               (default: 5.)
;   /CUAR     - CuAr lamps only!
;   /KLUDGE   - Force kludging of bottom end
;   /ONED     - Identify wavelengths using 1D solutions only (not
;               recommended)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echtrcarc, esi, 0.5, /CHK, /AUTO
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   12-Aug-2002 Written by JXP
;   21-Aug-2002 Streamlined + ps output
;   01-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------

pro esi_echtrcarc_ps, filnm, YVAL=yval, XVAL=XVAL, $
                      ORDR=ordr

  COMMON ps_common, old_dname, old_pfont, file_name, opened, color, ledger

  ;; Check open
  if keyword_set(opened) then begin
      if opened EQ 1 then ps_close, /noprint, /noid
  endif
  ;; Device
  device, get_decomposed=svdecomp
  device, decompose=0
  ;; Open
  ps_open, file=filnm, font=1, /color
  clr = getcolor(/load)
  maxy = max(yval,min=miny)
  ;; All points
  plot, xval, yval, psym=1, $
    background=clr.white, color=clr.black, $
    xtitle='x', ytitle='y', yrang=[miny,maxy], ystyle=1, xstyle=1
  ;; Order
  xyouts, 0.2, 1.03, 'Order = '+string(15-ordr, FORMAT='(i2)'), /normal, $
    charsize=2.0
  ps_close, /noprint, /noid
  device, decomposed=svdecomp
          
  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echtrcarc, esi, slit, LINLIST=linlist, CHK=chk, AUTO=auto, ORDR=ordr $
                   , DEBUG = debug, PKSIG = pksig, PIXSHFT = pixshft $
                   , NO11KLDG = no11kldg, CUAR = cuar, GUESSARC = guessarc $
                   , KLUDGE = kludge, CBIN = cbin, RBIN = rbin $
                   , LEDG = LEDG, REDG = REDG, PKWDTH=pkwdth, $
                   SEDG_FIL=sedg_fil

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echtrcarc, esi, slit, LINLIST=, /CHK, /AUTO, ORDR=, /DEBUG '
      print, '          PKSIG=, PIXSHFT=, /NO11KLDG, /CUAR, GUESSARC=, '
      print, '          /KLUDGE, CBIN=, RBIN=, PKWDTH=  [v1.1]'
      return
  endif 
  
;  Optional Keywords
  IF NOT KEYWORD_SET(LEDG) THEN LEDG = 0;9L
  IF NOT KEYWORD_SET(REDG) THEN REDG = 0 ;3L
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1
  if not keyword_set(MXSHFT) then mxshft = 25
  if not keyword_set( MAPCEN ) then mapcen = 2048L/rbin
  if not keyword_set( TRCNSIG ) then trcnsig = 3.
  if not keyword_set(MEDWID) then medwid = 5L
  if not keyword_set(PKSIG) then pksig = 5. 
  if not keyword_set( RADIUS ) then radius = 2.0
  if not keyword_set( MAP_FIL ) then $
    map_fil = getenv('ESI_CALIBS')+'/ECH_map.fits'
  ;; LINLIST
  if not keyword_set(LINLIST) then begin
      if keyword_set( CUAR ) then $
        linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ESI_CuArech.lst' $
      else linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/ESI_ech.lst'
  endif
  ;; GUESSARC
  if not keyword_set( GUESSARC ) then begin
      if keyword_set( CUAR ) then $
        guessarc = getenv('ESI_CALIBS')+'/ECH_CuArarcfit.idl' $
      else guessarc = getenv('ESI_CALIBS')+'/ECH_arcfit.idl' 
  endif
  ;; ORDR
  if not keyword_set( ORDR ) then begin
      qstrt = 0L
      if not keyword_set( CUAR ) then qend = 9L else qend = 8L
  endif else begin
      qstrt = ordr
      qend = ordr
  endelse

; Open line list
  x_arclist, linlist, lines, /GDONLY
  
; Set Peak Width
  if not keyword_set(PKWDTH) then begin
      case slit of
          0.3: pkwdth = 3L
          0.5: pkwdth = 4L
          0.75: pkwdth = 5L
          1.0: pkwdth = 6L
      endcase
  endif

; Grab Arc IMG

  c_s = esi_slitnm(slit)
  arc_fil = esi_getfil('arc_fil', SLIT=slit, cbin=cbin, rbin=rbin, /name)
  a = findfile(arc_fil+'*', count=na)
  if na EQ 0 then begin
      print, 'esi_echtrcarc: Arc ', arc_fil, ' doesnt exist. Run esi_echmkarc!'
      return
  endif
  print, 'esi_echtrcarc: Reading arc: ', arc_fil
  arc_img = xmrdfits(arc_fil, 0, head, /silent)
  arc_var = xmrdfits(arc_fil, 1, /silent)
  sz_arc = size(arc_img, /dimensions)

; AUTO
  if keyword_set(AUTO) then begin
      restore, guessarc
      guess_spec = temporary(sv_aspec)
      guess_fit = temporary(all_arcfit)
  endif

; Grab Arc Fit
  fit_fil = esi_getfil('arc_fit', SLIT=slit, cbin=cbin, rbin=rbin, /name)
  a = findfile(fit_fil, count=na)
  if na EQ 0 then begin
      print, 'esi_echtrcarc: Arc fit file ', fit_fil, ' doesnt exist.'
      print, 'esi_echtrcarc: Run esi_echfitarc first!'
      return
  endif
  print, 'esi_echtrcarc: Reading arc fit file: ', fit_fil
  restore, fit_fil              ; Variables: all_arcfit, sv_aspec

  ;; 2D
  if not keyword_set( ONED ) then begin
      fit2d_fil = esi_getfil('arc_2Dfit', SLIT=slit, $
                             cbin=cbin, rbin=rbin, /name)
      if x_chkfil(fit2d_fil+'*') EQ 0 then begin
          print, 'esi_echtrcarc: 2D Arc fit file ', fit_fil, ' doesnt exist.'
          print, 'esi_echtrcarc: Run esi_echfitarc first!'
          return
      endif
      print, 'esi_echtrcarc: Reading arc fit file: ', fit2d_fil
      arc_str = xmrdfits(fit2d_fil,1,/silent) ; Structure with poly coeff
  endif

; Open Map
  a = findfile(map_fil+'*', count=na)
  if na EQ 0 then begin
      print, 'esi_echtrcarc: Map file ', map_fil, ' doesnt exist.'
      print, 'esi_echtrcarc: Run esi_echmkmap first!'
      return
  endif
  print, 'esi_echtrcarc: Reading map file: ', map_fil
  map = xmrdfits(map_fil, /silent)

  ;; Binning
  if cbin NE 1 OR RBIN NE 1 then $
    map = rebin(map, 2048L/cbin, 4096L/rbin) / float(cbin)

  ;; Transpose the map
  map = transpose(map)

; AUTO
  if keyword_set( AUTO ) then begin
      gdx_fil = getenv('ESI_CALIBS')+ $
        '/ArcECH_'+c_s+'gdx.fits.gz'
      if x_chkfil(gdx_fil) EQ 0 then begin
          print, 'esi_echtrcarc: File ', gdx_fil, 'doesnt exist! Copy it!!'
          return
      endif
      all_gdx = xmrdfits(gdx_fil, /silent)
  endif else stop

; Open Slit file
  if not keyword_set(SEDG_FIL) then $
    slit_edg = esi_getfil('sedg_fil', SLIT=slit, cbin=cbin, rbin=rbin) $
  else slit_edg = xmrdfits(sedg_fil, 0, /silent)

  slit_cen = round((slit_edg[*,*,0] + slit_edg[*,*,1])/2.)
  rnd_edg = lonarr(size(slit_edg, /dim))
  rnd_edg[*, *, 0] = floor(slit_edg[*, *, 0])
  rnd_edg[*, *, 1] = ceil(slit_edg[*, *, 1])
  
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;  
; TRC
  for qq=qstrt,qend do begin
      print, 'esi_echtrcarc: Tracing order ', string(15L-qq, FORMAT='(i3)')
      ;; TMP IMG
      mnedg = min(rnd_edg[*,qq,0]) > 0L
      mxedg = max(rnd_edg[*,qq,1])
      tmp_img = arc_img[mnedg:mxedg, *]
      tmp_var = arc_var[mnedg:mxedg, *]
      ;; Zero out
      for j=0L,sz_arc[1]-1 do begin
          lhs = (rnd_edg[j, qq, 0]-mnedg+LEDG) > 0L ;; Funny LHS (0.5" only?)
          tmp_img[0:lhs,j] = 0. 
          tmp_var[0:lhs,j] = 0.
          rhs = (rnd_edg[j, qq, 1]-REDG-mnedg) > 0L
          tmp_img[rhs:mxedg-mnedg,j] = 0. 
          tmp_var[rhs:mxedg-mnedg,j] = 0.
      endfor
      ;; Rotate, Setup Image
      tmp_img = transpose(tmp_img)
      tmp_var = transpose(tmp_var)
      ;; Straighten
      rec_img = x_rectify(tmp_img, map[*,mnedg:mxedg], /silent)
      rec_var = x_rectify(tmp_var, map[*,mnedg:mxedg], /silent)

      ;; Inverse Variance
      rec_var[where(rec_var LE 0. OR rec_var GT 1e5)] = -1.
      ivar = 1./rec_var

      ;; Keep the good ones
      if keyword_set( AUTO ) then begin  ; AUTO
          tempmsk = bytarr(sz_arc[1])
          psig = pksig
          case qq of 
              0: tempmsk[1400:3790] = 1B
              1: begin
                  tempmsk[0:200] = 1B
                  tempmsk[800:sz_arc[1]-1] = 1B
              end
              2: begin
                  tempmsk[0:200] = 1B
                  tempmsk[700:3900] = 1B
              end
              3: begin
                  tempmsk[0:200] = 1B
                  tempmsk[500:sz_arc[1]-1] = 1B
              end
              4: begin
                  tempmsk[300:sz_arc[1]-1] = 1B
                  psig = 3.  ;; Trying to get first few lines!
              end
              else: tempmsk[0:sz_arc[1]-1] = 1B
          endcase

          ;; Find Peaks using arcs to find the shift
          if not keyword_set(PIXSHFT) then begin
              x_templarc, sv_aspec[*,qq], lines, guess_fit[qq], /FFT, $
                MSK=tempmsk, MXSHFT=mxshft, FORDR=5L, $
                MOCK_FFT=fft(guess_spec[*,qq]), SHFT=shft, ALL_PK=center, $
                PKWDTH=pkwdth, /THIN, PKSIG=psig, /SKIPLIN
          endif else begin
              if n_elements(pixshft) GT 1 then shft = pixshft[qq] $
              else shft = pixshft[0]
              x_templarc, sv_aspec[*, qq], lines, guess_fit[qq] $
                          , MSK = tempmsk, MOCK_FFT = fft(guess_spec[*, qq]) $
                          , SHFT = shft, ALL_PK = center, PKWDTH = pkwdth $
                          , /THIN, PKSIG = psig, /SKIPLIN, MXSHFT = MXSHFT
           endelse
          print, 'esi_echtrcarc: Shifted ', shft, ' pixels'
          ;; 
          msk = bytarr(n_elements(center)) 
          for i=0L,n_elements(center)-1 do begin
              if min(abs(center[i]-shft-all_gdx[*,qq])) LT 3. then msk[i] = 1B
          endfor
      endif else begin  ; INTERACTIVE
          stop
          ;; Find Peaks 
          x_fndpeaks, sv_aspec[*,qq], center, NSIG=pksig, /silent, $
            PKWDTH=pkwdth, /THIN
          print, 'esi_echtrcarc: Choose Arc lines for tracing!'
          x_prspeaks, sv_aspec[*,qq], center, msk, /block
      endelse

      gdtrc = where(msk EQ 1, npk)
      xstrt = fltarr(npk)
      ystrt = replicate(slit_cen[MAPCEN,qq]-mnedg,npk)

      ;; Center up peaks
      xstrt = trace_fweight(rec_img, center[gdtrc], ystrt, $
                            radius=1.5, invvar=ivar)
      for j=0L,9 do $
        xstrt = trace_fweight(rec_img, xstrt, ystrt, $
                              radius=radius, invvar=ivar)

      ;; Set wavelength
      if keyword_set( ONED ) then begin
          wav = 10^x_calcfit(double(xstrt), FITSTR=all_arcfit[qq])
      endif else begin
          wav = 10^(esi_echget2dwv(arc_str, double(xstrt), 15-qq))
;          wav2 = 10^x_calcfit(double(xstrt), FITSTR=all_arcfit[qq])
;          x_splot, xstrt, 10^wav-wav2, /block,psym1=1
      endelse
          
      if keyword_set( DEBUG ) and keyword_set( AUTO ) then begin
          x_prspeaks, sv_aspec[*,qq], xstrt, /block
      endif
      ;; TRACE CRUDE
      print, 'esi_echtrcarc: Tracing ', n_elements(xstrt), ' lines...'
      xcen = trace_crude(rec_img, ivar, yset=ycen, xstart=xstrt, $
                         radius=radius, ystart=ystrt, xerr=xerr, nave=5, $
                         maxshifte=0.2, maxshift0=1.0)
      ;; Tie down the bottom end!
      if (qq LE 5 AND qq GT 1) OR keyword_set( KLUDGE ) then begin
          print, 'esi_echtrcarc: Kludging the bottom end'
          sz_xcen = size(xcen, /dimensions)
          if keyword_set(ONED) then $
            wav0 = 10^x_calcfit(0.d, fitstr=all_arcfit[qq]) $
          else $
            wav0 = 10^esi_echget2dwv(arc_str, 0.d, 15-qq)
          xstrt = [xstrt, 0.]
          ystrt = [ystrt, slit_cen[MAPCEN,qq]-mnedg]
          wav = [wav, wav0]
          ;; xcen and xerr
          tmp = xcen
          tmp2 = xerr
          tmp3 = ycen
          ycen = fltarr(sz_xcen[0], sz_xcen[1]+1)
          xcen = fltarr(sz_xcen[0], sz_xcen[1]+1)
          xerr = fltarr(sz_xcen[0], sz_xcen[1]+1)
          xcen[*,0:sz_xcen[1]-1] = tmp
          xerr[*,0:sz_xcen[1]-1] = tmp2
          ycen[*,0:sz_xcen[1]-1] = tmp3
          xcen[*,sz_xcen[1]] = replicate(0., sz_xcen[0])
          ;; Only put it over the slit width at ths pointrnd_edg[j,qq,0]
          xerr[*,sz_xcen[1]] = 999.
          xerr[rnd_edg[MAPCEN,qq,0]-mnedg+9:rnd_edg[MAPCEN,qq,1]-mnedg-3, $
               sz_xcen[1]] = 0.2
          ycen[*,sz_xcen[1]] = ycen[*,0]
          sz_xcen = size(xcen, /dimensions)
      endif
      if qq EQ 4 AND not keyword_set(NO11KLDG) then begin
          print, 'esi_echtrcarc: Kludging 11 special'
          sz_xcen = size(xcen, /dimensions)
          if keyword_set( ONED ) then $
            wav5 = 10^x_calcfit(500.d, fitstr=all_arcfit[qq]) $
          else $
            wav5 = 10^esi_echget2dwv(arc_str, 500.d, 15-qq)
          xstrt = [xstrt, 0.]
          ystrt = [ystrt, slit_cen[MAPCEN,qq]-mnedg]
          wav = [wav, wav5]
          ;; xcen and xerr
          tmp = xcen
          tmp2 = xerr
          tmp3 = ycen
          ycen = fltarr(sz_xcen[0], sz_xcen[1]+1)
          xcen = fltarr(sz_xcen[0], sz_xcen[1]+1)
          xerr = fltarr(sz_xcen[0], sz_xcen[1]+1)
          xcen[*,0:sz_xcen[1]-1] = tmp
          xerr[*,0:sz_xcen[1]-1] = tmp2
          ycen[*,0:sz_xcen[1]-1] = tmp3
          xcen[*,sz_xcen[1]] = replicate(500., sz_xcen[0])
          ;; In good region only!
          xerr[*,sz_xcen[1]] = 999.
          xerr[rnd_edg[MAPCEN,qq,0]-mnedg+9:rnd_edg[MAPCEN,qq,1]-mnedg-3, $
               sz_xcen[1]] = 0.2
          ycen[*,sz_xcen[1]] = ycen[*,0]
          sz_xcen = size(xcen, /dimensions)
      endif
      trcstr = { $
                 xstrt: xstrt, $
                 ystrt: ystrt, $
                 sz_img: size(rec_img, /dimensions), $
                 xoff: mnedg, $
                 wav: wav, $
                 xerr: xerr, $
                 xcen: xcen $
               }
      sz_xcen = size(xcen,/dimensions)
      ;; CHK
      if keyword_set( CHK ) then begin
          x_splot, xcen[rnd_edg[2048L/cbin,qq,0]-mnedg:sz_xcen[0]-1,*], $
            ycen[rnd_edg[2048L,qq,0]-mnedg:sz_xcen[0]-1,*], psym1=1, /block
      endif
      ;; DEBUG
      if keyword_set( DEBUG ) then stop
      ;; OUTPUT
      ordr = 15L - qq
      outfil = esi_getfil('arc_trc', SLIT=slit, cbin=cbin, rbin=rbin, ORDR=ordr, $
                          /name)
      print, 'esi_echtrcarc: Outputing trace structure: ', outfil
      mwrfits, trcstr, outfil, /create, /silent

      ;; PSFIL
      psfil = esi_getfil('arc_trcps', SLIT=slit, cbin=cbin, rbin=rbin, ORDR=ordr, $
                          /name)
      esi_echtrcarc_ps, psfil, $
          XVAL=xcen[rnd_edg[2048L/cbin,qq,0]-mnedg:sz_xcen[0]-1,*], $
            YVAL=ycen[rnd_edg[2048L/cbin,qq,0]-mnedg:sz_xcen[0]-1,*], ORDR=qq
  endfor
  ;;
  print, 'esi_echtrcarc: All done!'
      
  return
end
