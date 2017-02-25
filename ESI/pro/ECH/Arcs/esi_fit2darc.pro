;+ 
; NAME:
; esi_fit2darc   
;     Version 1.1
;
; PURPOSE:
;    Fit a 2D solution to the arc lines
;
; CALLING SEQUENCE:
;   
;  esi_fit2darc, esi
;
; INPUTS:
;   esi     -  ESI structure
;   [slit]  -  Slit size (e.g. 0.5, 0.75, 1.0)
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   /IFLAT   - Use internal flats
;   /REDOOV  - Overwrite OV files if they exist for the flats
;   /SVOV    - Save the OV files created during this step
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_fit2darc, esi, 1.0
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_fit2darc, esi, slit, CLOBBER=clobber, CBIN=cbin, RBIN=rbin, $
                   nycoeff=nycoeff, nocoeff=nocoeff, DEBUG=debug, CHKRES=chkres
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_fit2darc, esi, slit, nycoeff=, nocoeff= [v1.0]'
      return
  endif

;  Optional Keywords
  if NOT keyword_set(nycoeff) then nycoeff = 7
  if NOT keyword_set(nocoeff) then nocoeff = 8
  if not keyword_set( CBIN ) then cbin = 1
  if not keyword_set( RBIN ) then rbin = 1

   ;;  Check for outfil
  outfil = esi_getfil('arc_2Dfit', slit=slit, /name, CHKFIL=chkf, $
                      cbin=cbin, rbin=rbin)
  if chkf NE 0 and not keyword_set( CLOBBER ) then begin
      print, 'esi_fit2darc: Arc fit file exists. ' + $
        'Continuing..'
      return
  endif

; Open arc 1D fit
   arc_info = esi_getfil('arc_fit', slit=slit, cbin=cbin, rbin=rbin, /name)
   restore, arc_info

   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; SETUP THE DATA

   gd_ordr = where(sv_lines.nlin NE 0, ngd_ordr)
   print, 'esi_fit2darc: Setting up the values (normalizing)'
   npix = round(total(sv_lines.nlin))
   t = dblarr(npix)
   ;; PIX AND WV
   cnt = 0L
   for j=0L,ngd_ordr-1 do begin
       ;; Dummy indx
       i = gd_ordr[j]
       ;; PIX
       if j EQ 0 then all_pix = [sv_lines[i].pix[0:sv_lines[i].nlin-1]] $
       else all_pix = [all_pix,sv_lines[i].pix[0:sv_lines[i].nlin-1]]
       ;; WV
       if j EQ 0 then all_wv = [sv_lines[i].wv[0:sv_lines[i].nlin-1]] $
       else all_wv = [all_wv,sv_lines[i].wv[0:sv_lines[i].nlin-1]]
       ;; Order #
       t[cnt:cnt+sv_lines[i].nlin-1] = 15-j
       cnt = cnt + sv_lines[i].nlin
   endfor

   nrm = dblarr(2)
   ;; NORMALIZE PIX
   mnx = min(all_pix, MAX=mxx)
   nrm[0] = 0.5 * (mnx + mxx)
   nrm[1] = mxx - mnx
   pix_nrm = 2. * (all_pix - nrm[0])/nrm[1]

   ;; NORMALIZE ORDER
   nrmt = dblarr(2)
   mnx = min(t, MAX=mxx)
   nrmt[0] = 0.5 * (mnx + mxx)
   nrmt[1] = mxx - mnx
   t_nrm = 2. * (t - nrmt[0])/nrmt[1]
   
   invvar = replicate(1., npix)

;  Setup the Functions
   work2d = dblarr(npix,nycoeff*nocoeff)
   worky = fchebyshev(pix_nrm[*], nycoeff)
   workt = fchebyshev(t_nrm[*], nocoeff)
   
   for i=0,nocoeff-1 do begin
       for j=0,nycoeff-1 do begin
           work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
       endfor
   endfor

   ;; Do the matrix algebra
   work2di = transpose(work2d * (invvar[*] # replicate(1,nocoeff*nycoeff)))
   alpha = work2di # work2d
   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p
   res = cholsol(alpha,p,beta, /double)
   wv_mod = dblarr(npix)
   wv_mod[*] = work2d # res

   ;; Get RMS
   gd_wv = where(invvar GT 0.0, ngd)
   msk = bytarr(npix)
   msk[gd_wv] = 1B
   ;; REJECT
   djs_iterstat, (wv_mod-alog10(all_wv)), sigrej=3.0, mask=msk
   gd = where(msk EQ 1B, complement=bad, ncomplement=nbad)

   ;; RESET invvar
   if nbad NE 0 then begin
       print, 'esi_fit2darc_work: Rejecting ', nbad, ' of ', n_elements(all_wv), $
         ' lines'
       invvar[bad] = 0.
   endif


   ;; Do the matrix algebra
   work2di = transpose(work2d * (invvar[*] # replicate(1,nocoeff*nycoeff)))
   alpha = work2di # work2d
   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p
   res = cholsol(alpha,p,beta, /double)
   wv_mod = all_wv * 0.0
   wv_mod[*] = work2d # res


   ;; Finish
   gd_wv = where(invvar GT 0.0, ngd)
   resid = 10^wv_mod[gd_wv] - all_wv[gd_wv]
   rms = sqrt( total( resid^2 ) / float(ngd))
   print, 'esi_fit2darc: RMS = ', rms, ' Ang'
;   if keyword_set(CHKRES) then x_splot, all_wv[gd_wv], resid, psym1=1, /block
   if keyword_set(CHKRES) then x_splot, resid, /block

   ;; DEBUG
   if keyword_set( DEBUG ) then begin
       x_splot, all_wv[gd], 10^wv_mod[gd]-all_wv[gd], psym1=1, $
         XTWO=all_wv[bad], YTWO=10^wv_mod[bad]-all_wv[bad], psym2=2, /block
       stop
   endif


;  OUTPUT
   out_str = { $
               nrm: nrm, $
               nrmt: nrmt, $
               ny: nycoeff, $
               no: nocoeff, $
               res: res }
   
   ;;  OUTPUT
   print, 'esi_fit2darc: Writing '+outfil
   mwrfits, out_str, outfil, /create

   print, 'esi_fit2darc: All done'
   return

end

