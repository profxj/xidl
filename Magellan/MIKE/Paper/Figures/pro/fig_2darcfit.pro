;+ 
; NAME:
; fig_2darcfit
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in x_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  x_fit2darc.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
;   
;  x_fit2darc
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
; A fits file containing the 2D solution.  Named something like
; 'Arcs/Fits/Arc_mb0539_fit2D.fits' 
;
; OPTIONAL KEYWORDS:
;   NOCOEFF   - Number of coefficients to use in the x-direction
;               (default: 6 for blue, 7 for red)
;   NYCOEFF   - Number of coefficients to use in the y-direction
;               (default: 4)
;   /CLOBBER  - Overwrite any previous solution
;   /DEBUG    - debug
;   /CHKRES   - Plot the residuals
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   print, x_fit2darc(arcfil, ordr_str, arc_info)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;   Feb-2005 Ported to XIDL by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fig_2darcfit, arcfil, ordr_str, arc_info, nycoeff=nycoeff, $
                  nocoeff=nocoeff, CLOBBER=clobber, SZ=sz, $
                  DEBUG=debug, CHKRES=chkres, out_str=out_str, $
                  YVAL=yval, ORIG=orig, SIGREJ=sigrej
                   
  if not keyword_set(PSFIL) then psfil='fig_2darcfit'
  if not keyword_set(ARCFIL) then arcfil = $
    getenv('MIKE_PAP')+'Arcs/Arc_mb0005.fits'
  if not keyword_set(OSTRFIL) then ostrfil = $
    getenv('MIKE_PAP')+'Flats/OStr_B_01.fits'
  ordr_str = xmrdfits(ostrfil, 1, /silent)
  if not keyword_set(ARC_INFO) then arc_info = $
    getenv('MIKE_PAP')+'Arcs/Fits/mb0005_fit.idl'
  lsz = 1.3

;
;  if  N_params() LT 3  then begin 
;      print,'Syntax - ' + $
;        'rslt = x_fit2darc( arcfil, ordr_str, arc_info, OUT_STR=) [v1.1]'
;      return, -1
;  endif

   if NOT keyword_set(nycoeff) then nycoeff = 4
   if NOT keyword_set(nocoeff) then nocoeff=5
   if NOT keyword_set(SIGREJ) then sigrej = 10.
;   if NOT keyword_set(nycoeff) then nycoeff = 4
;   if NOT keyword_set(nocoeff) then nocoeff=6

   ;; Arc info
   restore, arc_info

   
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; 
; SETUP THE DATA

   gd_ordr = where(sv_lines.nlin NE 0 AND ordr_str.flg_anly NE (-1), ngd_ordr)
   print, 'x_fit2darc: Setting up the values (normalizing)'
   npix = round(total(sv_lines.nlin))
   omsk = lonarr(npix)
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
       if keyword_set( ORIG ) then begin
           if j EQ 0 then all_wv = [sv_lines[i].wv[0:sv_lines[i].nlin-1]] $
           else all_wv = [all_wv,sv_lines[i].wv[0:sv_lines[i].nlin-1]]
           stop
       endif else begin
           if j EQ 0 then all_wv = $
             [sv_lines[i].wv[0:sv_lines[i].nlin-1]]*ordr_str[i].order $
           else all_wv = [all_wv,sv_lines[i].wv[0:sv_lines[i].nlin-1] $
                          *ordr_str[i].order]
       endelse
       ;; t
       if not keyword_set( YVAL ) then begin
           ;; Order #
           t[cnt:cnt+sv_lines[i].nlin-1] = ordr_str[i].order
       endif else begin
           ;; x position
           xcen = (ordr_str[i].lhedg + ordr_str[i].rhedg)/2.
           np = n_elements(ordr_str[i].lhedg)
           spl = spl_init(findgen(np), xcen)
           ;; Evaluate
           t[cnt:cnt+sv_lines[i].nlin-1] = $
             spl_interp( findgen(np), xcen, $
                         spl, sv_lines[i].pix[0:sv_lines[i].nlin-1]) 
           
;           t[cnt:cnt+sv_lines[i].nlin-1] = $
;             sv_lines[i].pix[0:sv_lines[i].nlin-1] * ordr_str[i].order
       endelse
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
   worky = flegendre(pix_nrm[*], nycoeff)
   workt = flegendre(t_nrm[*], nocoeff)
   
   for i=0,nocoeff-1 do begin
       for j=0,nycoeff-1 do begin
           work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
       endfor
   endfor

   ;; Do the matrix algebra
   work2di = transpose(work2d * (invvar[*] # replicate(1,nocoeff*nycoeff)))
   alpha = work2di # work2d
   beta = work2di # all_wv[*]
;   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p
   res = cholsol(alpha,p,beta, /double)
   wv_mod = dblarr(npix)
   wv_mod[*] = work2d # res

   ;; Get RMS
   gd_wv = where(invvar GT 0.0, ngd)
   msk = bytarr(npix)
   msk[gd_wv] = 1B
   ;; RESID
;   resid = (wv_mod[gd_wv] - all_wv[gd_wv]);/t[gd_wv]  ; Ang
;   fin_rms = sqrt( total( resid^2 ) / float(ngd))
;   print, 'x_fit2darc: RMS = ', fin_rms, ' Ang*Order#'
;   if keyword_set(CHKRES) then x_splot, all_wv[gd_wv], resid, psym1=1, /block
;   if keyword_set(CHKRES) then x_splot, resid, /block

   ;; REJECT
;   djs_iterstat, (wv_mod-alog10(all_wv)), sigrej=2.5, mask=msk
   djs_iterstat, (wv_mod-all_wv), sigrej=sigrej, mask=msk
   gd = where(msk EQ 1B, complement=bad, ncomplement=nbad)

   ;; RESET invvar
   if nbad NE 0 then begin
       print, 'x_fit2darc_work: Rejecting ', $
         nbad, ' of ', n_elements(all_wv), ' lines'
       invvar[bad] = 0.
   endif

   ;; Do the matrix algebra
   work2di = transpose(work2d * (invvar[*] # replicate(1,nocoeff*nycoeff)))
   alpha = work2di # work2d
   beta = work2di # (all_wv[*])
;   beta = work2di # (alog10(all_wv[*]))
   choldc, alpha, p
   res = cholsol(alpha,p,beta, /double)
   wv_mod = all_wv * 0.0
   wv_mod[*] = work2d # res


   ;; Finish
   gd_wv = where(invvar GT 0.0, ngd)
   resid = (wv_mod[gd_wv] - all_wv[gd_wv]);/t[gd_wv]  ; Ang
;   resid = 10^wv_mod[gd_wv] - all_wv[gd_wv]
   fin_rms = sqrt( total( resid^2 ) / float(ngd))
   print, 'x_fit2darc: RMS = ', fin_rms, ' Ang*Order#'
;   if keyword_set(CHKRES) then x_splot, all_wv[gd_wv], resid, psym1=1, /block
   if keyword_set(CHKRES) then x_splot, resid, /block


   ;; OUTPUT Structure
   out_str = { $
               nrm: nrm, $
               nrmt: nrmt, $
               ny: nycoeff, $
               no: nocoeff, $
               res: res }
   

   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
   ;; QA

   if not keyword_set(SZ) then begin
       head = xheadfits(arcfil, /silent)
       sz = lonarr(2)
       sz[0] = sxpar(head,'NAXIS1')
       sz[1] = sxpar(head,'NAXIS2')
   endif
   
   ;; 2D
   x_psopen, psfil+'a.ps', /maxs
   clr = getcolor(/load)
       
   ;; NORMALIZE PIX
   all_pix = dindgen(round(sz[1]/8.))*8.
   npix = n_elements(all_pix)
   pix_nrm = 2. * (all_pix - nrm[0])/nrm[1]
   worky = flegendre(pix_nrm[*], nycoeff)
       
   ;; Main plot
   !p.multi=[0,1,1] 
   mn = min(wv_mod/t, max=mx)
;       mn = min(10^wv_mod, max=mx)
   plot, [0.], [0.], color=clr.black, $
         background=clr.white, charsize=1.5, xrange=[mn-30., mx+30.], $
         yrange=[0.,sz[1]+140], xstyle=1, ystyle=1, ytitle='CCD Row', $
         xtitle='Wavelength (Ang)', xmargin=[8,2], ymargin=[5,2], /nodata
   
   nordr = n_elements(ordr_str)
   for jj=0L,nordr-1 do begin
       ;; NORMALIZE ORDER
       ii = ordr_str[jj].order
       tsub = replicate(float(ii), npix)
       t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
       
       ;; work2d and wv
       work2d = dblarr(npix,nycoeff*nocoeff)
       workt = flegendre(t_nrm[*], nocoeff)
       
       for i=0,nocoeff-1 do begin
           for j=0,nycoeff-1 do begin
               work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
           endfor
       endfor
       
       wv = dblarr(npix)
;           wv[*] = work2d # out_str.res 
       wv[*] = work2d # out_str.res / ordr_str[jj].order
       
       mn = min(wv, max=mx)
       
       oplot, wv, all_pix, color=clr.blue
           
       ;; Resid
       pts = where(t EQ ii, npts)
       if npts NE 0 then begin
           nlin = sv_lines[jj].nlin
           sres = (wv_mod[pts] - all_wv[pts])/ordr_str[jj].order
           oplot, wv_mod[pts]/t[pts] + sres*500., $
                  sv_lines[jj].pix[0:nlin-1],  $
                  psym=1, color=clr.black
           ;; RMS
;           rms = sqrt( total( sres^2 ) / float(nlin))
       endif
           
   endfor
   xyouts, 0.5, 0.91, '(a) 2D Fit (Res x500)   nx='+strtrim(nocoeff,2)+ $
           ' ny='+strtrim(nycoeff,2)+' ' + $
           '  RMS='+string(fin_rms,format='(f6.4)')+'Ang*Order#', $
           color=clr.black, charsize=1.5, /normal, alignment=0.5
           
   ;;;;;;;;;;
   ;; Individual plots
   x_psclose

   x_psopen, psfil+'b.ps', /maxs
   clr = getcolor(/load)
   nordr = n_elements(ordr_str)
   !p.multi=[0,3,2]
   nordr = 6
   for jj=0L,nordr-1 do begin
       ;; NORMALIZE ORDER
       ii = ordr_str[jj].order
       tsub = replicate(float(ii), npix)
       t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
       
       ;; work2d and wv
       work2d = dblarr(npix,nycoeff*nocoeff)
       workt = flegendre(t_nrm[*], nocoeff)
       
       for i=0,nocoeff-1 do begin
           for j=0,nycoeff-1 do begin
               work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
           endfor
       endfor
       
       wv = dblarr(npix)
;           wv[*] = work2d # out_str.res
       wv[*] = work2d # out_str.res / tsub
       
;           mn = min(10^wv, max=mx)
;           plot, all_pix, 10^wv, color=clr.black, $
       mn = min(wv, max=mx)
       plot, all_pix, wv, color=clr.black, $
             background=clr.white, charsize=2.3, yrange=[mn, mx], $
             xrange=[0.,sz[1]], xstyle=1, ystyle=1, xtitle='CCD Row', $
             ytitle='Wavelength (Ang)', xmargin=[8,1], ymargin=[4,0.5], /nodata
       
       ;; Fit
       oplot, all_pix, wv, color=clr.blue
;           oplot, all_pix, 10^wv, color=clr.blue
       
       ;; Resid
       pts = where(t EQ ii, npts)
       rms = 9.99
       if npts NE 0 then begin
           nlin = sv_lines[jj].nlin
;               sres = 10^wv_mod[pts] - all_wv[pts] 
;               oplot, sv_lines[jj].pix[0:nlin-1],  10^wv_mod[pts] + sres*100., $
           sres = (wv_mod[pts] - all_wv[pts])/t[pts]
           oplot, sv_lines[jj].pix[0:nlin-1],  wv_mod[pts]/t[pts] $
                  + sres*100., $
                  psym=1, color=clr.black
           ;; Rej
           rej = where(invvar[pts] LE 0., nrej, ncomplement=nnorej, $
                       complement=norej)
           if nrej NE 0 then $
             oplot, [(sv_lines[jj].pix[0:nlin-1])[rej]],  $
                    [(all_wv[pts]/t[pts] + sres*100)[rej]], $
                    psym=2, color=clr.red
;                       [(10^wv_mod[pts] + sres*100.)[rej]], $
           ;; RMS
           if nnorej NE 0 then $
             rms = sqrt( total( sres[norej]^2 ) / float(nnorej-1)) $
           else rms = 9.99
       endif
           
       ;; Label
       ylbl = mn + (mx-mn)*0.08*(findgen(5)+1) 
       
       xyouts, sz[1]*0.05, ylbl[2], 'Order = '+strtrim(ii,2), $
               color=clr.black, charsize=lsz
       dwv = abs(wv[0]-wv[npix-1])/float(all_pix[npix-1])
;           dwv = abs(10^wv[0]-10^wv[npix-1])/float(all_pix[npix-1])
       xyouts, sz[1]*0.05, ylbl[0], '!9Dl!X = '+string(dwv,format='(f6.4)'), $
               color=clr.black, charsize=lsz
       xyouts, sz[1]*0.05, ylbl[1], 'RMS(pix) = '+$
               string(rms/dwv,format='(f4.2)'), color=clr.black, charsize=lsz
   endfor
   
   x_psclose
   !p.multi=[0,1,1]

   return

end

