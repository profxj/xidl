;+ 
; NAME:
; wfccd_arcsol
;    Version 1.0
;
; PURPOSE:
;    Solves arc solutions for a given mask
;      Designed to do only 1 at a time
;
; CALLING SEQUENCE:
;   
;   wfccd_arcsol, wfccd, WFARC=
;
; INPUTS:
;   wfstrct     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   wfarc      -  WFCCD arc structure (fits file)
;
; OPTIONAL INPUTS:
;   w0off      - half-range of lambda zeropoint to check on initial guess
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_arcsol, wfstrct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   14-Feb-2002 Written by JXP
;   12-Jul-2002 Modified by JXP: Changed the fit structure
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Initialize ths common blocks
pro wfccd_arcsol_init, xslit, guess, nguess

  common wfccd_arcsol_common, w1fit, w2fit, w3fit

  ; 1st
  w1fit = { fitstrct }

  ; 2nd
  w2fit = { fitstrct }
  w2fit.func = 'POLY'
  w2fit.nord = 3
  ffit = [-0.00010227686d, -1.7738333d-05, 2.6352564d-05, -5.8507298d-05]
  w2fit.ffit = ptr_new(ffit)
  w2fit.nrm = [998.55432129d, 1246.80883789d]

  ; 3rd
  w3fit = { fitstrct }
  w3fit.func = 'POLY'
  w3fit.nord = 2
  ffit = [1.1878507e-07, -5.1938572e-08, 2.2261541d-08]
  w3fit.ffit = ptr_new(ffit)
  w3fit.nrm = [998.55432129d, 1246.80883789d]

  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Uninitialize
pro wfccd_arcsol_uninit

  common wfccd_arcsol_common
  delvarx, w1fit, w2fit, w3fit
  return
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Setups the Guess for a Slit given the x position (in pixels)
pro wfccd_arcsol_setup, xslit, guess, nguess, W0OFF=w0off, W2OFF=w2off

  common wfccd_arcsol_common

  ; Optional keywords
  if not keyword_set( W0OFF ) then w0off = 100.
  if not keyword_set( W2OFF ) then w2off = 2d-5
  ; Create guess
  guess = dblarr(2,4)
  nguess = lonarr(3)

  ; Zero term
  nguess[0] = round(2.*w0off/3.)  ; Loop over as many pixels as spread
  guess[0,0] = (4874.3883 + 3.1626093*xslit - w0off) 
  guess[1,0] = guess[0,0] + 2*w0off

  ; First term --  No fit for now
  nguess[1] = 20
  guess[0,1] = -3.20
  guess[1,1] = -3.13

  ; 2nd term
  val = x_calcfit(xslit, FITSTR=w2fit)

  nguess[2] = 80
  guess[0,2] = val - w2off
  guess[1,2] = guess[0,2]+2.*w2off
  
  ; 3rd term
  guess[0,3] = x_calcfit(xslit, FITSTR=w3fit)

  return
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  MAIN DRIVER
pro wfccd_arcsol, wfccd, obj_id, ARC=arc, WFARC=wfarc, SLITSTR=slitstr, $
                  OUTFIL=outfil, NOFITS=nofits, CLOBBER=clobber, W0OFF=w0off, $
                  W2OFF=w2off

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_arcsol, wfccd, mask_id, AIMG=, WFARC=, SLITSTR=, OUTFIL= '
    print, '          /NOFITS, /CLOBBER [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set(OUTFIL) then begin
      ; Key on R_
      ipos = strpos(wfccd[obj_id].arc_fil, 'R_')
      ; Append
      outfil = 'Arcs/ArcS_'+strmid(wfccd[obj_id].arc_fil,ipos+2)
  endif

; Check to see if OUTFIL exists
  a = findfile(outfil, count=count)
  if count GT 0 AND not keyword_set( CLOBBER ) then begin
      print, 'wfccd_arcsol: ', outfil, ' exists, returning'
      return
  endif

;  Initialize the fits

  wfccd_arcsol_init

;  Read in the line list
  linelist = $
    getenv('XIDL_DIR')+'/Spec/Arcs/Lists/wfccdB_HeNe.lst'
  x_arclist, linelist, lines

;  Set slit structure
  if not keyword_set( SLITSTR ) then $
    slitstr = xmrdfits(wfccd[obj_id].slit_fil,1, STRUCTYP='mslitstrct', /silent)

;  Read the Arc
  if not keyword_set( ARC ) then arc = xmrdfits(wfccd[obj_id].arc_fil, /silent)

; Find slits
  gd = where(slitstr.flg_anly NE 0, ngd)
  nslit = n_elements(slitstr)
  
;  Create structure
  tmp = { wfccdarcstr }
  wfarc = replicate(tmp, nslit)

  ;; FITSTR
      
;  Loop

  if not keyword_set( SILENT ) then $
    print, 'wfccd_arcsol: Looping on the slits...'
  for i=0L,ngd-1 do begin
      if not keyword_set( SILENT ) then $
        print, 'wfccd_arcsol: Slit ', strtrim(i,2), ' of ', strtrim(ngd-1,2), $
        ' ( '+strtrim(string(gd[i]),2)+' )'
      ;; Center of slit
      cent = total(slitstr[gd[i]].yedg_flt)/2.
;; 0 Slit kludge -- slits should be right
; if i EQ 0 AND cent GT 1850L then cent = slitstr[gd[i]].yedg_flt[0]+5
; if i EQ (ngd-1) AND cent LT 100L then cent = slitstr[gd[i]].yedg_flt[1]-5
      wfarc[gd[i]].cent = cent

      ;; Take center 5 rows
      spec = djs_median(arc[*,cent-2:cent+2], 2)
      npix = n_elements(spec)

      ;; Setup autoid
      wfccd_arcsol_setup, slitstr[gd[i]].xpos, guess, nguess, w0off=w0off, $
        w2off=w2off
      
      ; Kill RHS if GUESS < 7000A (KLUDGE)
      if guess[0,0] LT 7000. then spec[1800:npix-1] = spec[1799]

      ;; Reset Fit 
      if keyword_set( fitstr ) then delvarx, fitstr
      fitstr = {fitstrct}
      fitstr.func = 'POLY'
      fitstr.lsig = 3.
      fitstr.hsig = 3.
      fitstr.niter = 3
      fitstr.minpt = 5
      fitstr.maxrej = 5

      ;; Order (hacked to put limit on how far the right it can go)
      if (guess[0,0] GT 5800 AND guess[0,0] LT 10000.) then $
        fitstr.nord = 6 else fitstr.nord = 4

      ;; Run x_autoid
      x_autoid, spec, GUESS=guess, NGUESS=nguess, LINES=lines, /cprog, $
        /silent, FINFIT=fitstr, FLG_AID=flg_aid, STATUS=status
      if(status eq 5) then begin
          splog, 'x_autoid failed, trying again with LBUFFER=100'
          x_autoid, spec, GUESS=guess, NGUESS=nguess, LINES=lines, /cprog, $
            /silent, FINFIT=fitstr, FLG_AID=flg_aid, LBUFFER=100., $
            STATUS=status
          if(status eq 5) then begin
              splog, 'x_autoid failed, trying again with LBUFFER=150'
              x_autoid, spec, GUESS=guess, NGUESS=nguess, LINES=lines, $
                /cprog, /silent, FINFIT=fitstr, FLG_AID=flg_aid, $
                LBUFFER=150., STATUS=status
          endif
      endif 
      if(status ne 0) then begin
          message, 'x_autoid returned failure!'
      endif

      ;; Check Flag
      if flg_aid EQ 0 then begin
          print, 'wfccd_arcsol: Bad solution for slit ', i 
          print, 'wfccd_arcsol: Run wfccd_fixarc!!'
          wfarc[gd[i]].flg_anly = 0
      endif else begin ;; Save fitstr
          wfarc[gd[i]].wave = x_calcfit(findgen(npix), FITSTR=fitstr)
          wfarc[gd[i]].spec = temporary(spec)
          wfarc[gd[i]].fit = temporary(fitstr)
          wfarc[gd[i]].flg_anly = 1
      endelse
  endfor


  ;; Output fits image
  if not keyword_set( NOFITS ) then wfccd_writeastrct, wfarc, outfil
  
  ;; Free memory
  wfccd_arcsol_uninit

  if not keyword_set( SILENT ) then print, 'wfccd_arcsol: All done'
  
  return
end
