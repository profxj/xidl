;+ 
; NAME:
; x_fit2darc   
;     Version 1.1
;
; PURPOSE:
;  To fit the arc lines identified in x_fitarc as a fucntion of
;  their y-centroid and order number.  The main routine is in
;  x_fit2darc.  The fit is a simple least-squares with one
;  round of rejection.
;
; CALLING SEQUENCE:
; x_fit2darc, arcfil, ordr_str, arc_info, nycoeff=nycoeff, $
;                    nocoeff=nocoeff, CLOBBER=clobber, SZ=sz, $
;                    DEBUG=debug, CHKRES=chkres, out_str=out_str, $
;                    QAFIL=qafil, YVAL=yval, ORIG=orig, SIGREJ=sigrej $
;                    , pixrms = pixrms
;   
;
; INPUTS:
;  arc_fil  -- Name of arc file
;  ordr_str -- Order strucure describing the echelle footprint
;  arc_info -- Name of file created by x_fitarc
;
; RETURNS:
;
; OUTPUTS:
; A fits file containing the 2D solution.  Named something like
; 'Arcs/Fits/Arc_B0010_fit2D.fits' 
;
; OPTIONAL KEYWORDS:
;   NOCOEFF   - Number of coefficients to use in the x-direction
;               (default: 5)
;   NYCOEFF   - Number of coefficients to use in the y-direction
;               (default: 4)
;   /CLOBBER  - Overwrite any previous solution
;   /DEBUG    - debug
;   /CHKRES   - Plot the residuals
;   /YVAL     - Use the x position on the CCD [not recommended]
;   /ORIG     - Use wavelengths for the basis not wavelength*order#
;               [not recommended]
;   SIGREJ=   - Value for sigma rejection in the fit [default: 10]
;
; OPTIONAL OUTPUTS:
;  PIXRMS=  -- RMS of the fit in pixels (for each echelle order)
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
function m_fit2darc, arcfil, ordr_str, arc_info, nycoeff=nycoeff, $
                     nocoeff=nocoeff, CLOBBER=clobber, SZ=sz, $
                     DEBUG=debug, CHKRES=chkres, out_str=out_str, $
                     QAFIL=qafil, YVAL=yval, ORIG=orig, SIGREJ=sigrej $
                     , pixrms = pixrms
                    
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = x_fit2darc( arcfil, ordr_str, arc_info, OUT_STR=) [v1.1]'
      return, -1
  endif

   if NOT keyword_set(nycoeff) then nycoeff = 5
   if NOT keyword_set(nocoeff) then nocoeff=6
   if NOT keyword_set(SIGREJ) then sigrej = 2.

   ;; Arc info
   restore, arc_info

   x_psclose
   
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
       sv_lines[i].nlin = min([sv_lines[i].nlin,n_elements(sv_lines[i].pix)])
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
   t = t[0:cnt-1]
   nrmt = dblarr(2)
   mnx = min(t, MAX=mxx)
   nrmt[0] = 0.5 * (mnx + mxx)
   nrmt[1] = mxx - mnx
   t_nrm = 2. * (t - nrmt[0])/nrmt[1]
   
   invvar = replicate(1., npix)

;  Setup the Functions
;   work2d = dblarr(npix,nycoeff*nocoeff)
   work2d = dblarr(cnt,nycoeff*nocoeff)
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
   choldc, alpha, p, /double
   res = cholsol(alpha,p,beta, /double)
;   wv_mod = dblarr(npix)
   wv_mod = dblarr(cnt)
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
   print, 'x_fit2darc: RMS = ', fin_rms, ' Ang * Order#'
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
   if keyword_set( QAFIL ) then begin

       if not keyword_set(SZ) then begin
           head = xheadfits(arcfil, /silent)
           sz = lonarr(2)
           sz[0] = sxpar(head,'NAXIS1')
           sz[1] = sxpar(head,'NAXIS2')
       endif

       x_psopen, qafil, /maxs
       clr = getcolor(/load)

       
;       head = xheadfits(arcfil, /silent)
;       sz = lonarr(2)
;       sz[0] = sxpar(head, 'NAXIS1')        
;       sz[1] = sxpar(head, 'NAXIS2')
       ;; NORMALIZE PIX
       all_pix = dindgen(round(sz[1]))
;       all_pix = dindgen(round(sz[1]/8.))*8.
       npix = n_elements(all_pix)
       pix_nrm = 2. * (all_pix - nrm[0])/nrm[1]
       worky = flegendre(pix_nrm[*], nycoeff)
       
       ;; Main plot
       !p.multi=[0,1,1] 
       mn = min(wv_mod/t, max=mx)
;       mn = min(10^wv_mod, max=mx)
       
       plot, [0.], [0.], color=clr.black, $
         background=clr.white, charsize=1.5, xrange=[mn-30., mx+30.], $
         yrange=[0.,sz[1]], xstyle=1, ystyle=1, ytitle='Row', $
         xtitle='Wavelength ('+ STRING(197B)+')' , xmargin=[11,2], ymargin=[5,2], /nodata
       
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
;           if jj EQ 19 then stop
           
           mn = min(wv, max=mx)
;           mn = min(10^wv, max=mx)
;       all_x = (ordr_str[jj].lhedg + $
;         ordr_str[jj].rhedg) / 2.
           
           oplot, wv, all_pix, color=clr.blue
;           oplot, 10^wv, all_pix, color=clr.blue
           
           ;; Resid
           pts = where(t EQ ii, npts)
           if npts NE 0 then begin
               nlin = sv_lines[jj].nlin
;               sres = 10^wv_mod[pts] - all_wv[pts] 
;               oplot, 10^wv_mod[pts] + sres*500., $
               sres = (wv_mod[pts] - all_wv[pts])/ordr_str[jj].order
               oplot, wv_mod[pts]/t[pts] + sres*1., $
                 sv_lines[jj].pix[0:nlin-1],  $
                 psym=1, color=clr.black
               ;; RMS
;           rms = sqrt( total( sres^2 ) / float(nlin))
           endif
           
           ;; Label
           xyouts, 0.5, 0.96, 'Arc 2D Fit (Res x500) nx='+strtrim(nocoeff,2)+ $
                   ' ny='+strtrim(nycoeff,2)+' ' + $
                   'RMS='+string(fin_rms,format='(f7.4)')+ STRING(197B)+'*Order#', $
             color=clr.black, charsize=1.5, /normal, alignment=0.5
       endfor
           
   ;;;;;;;;;;
       ;; Individual plots
       nordr = n_elements(ordr_str)
       !p.multi = [0, 3, 2]
       pixrms = fltarr(nordr)
       for jj=0L,nordr-1 do begin
           ;; NORMALIZE ORDER
           ii = ordr_str[jj].order
           tsub = replicate(float(ii), npix)
;           tsub = replicate(float(ii), cnt)
           t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
           
           ;; work2d and wv
           work2d = dblarr(npix,nycoeff*nocoeff)
;           work2d = dblarr(cnt,nycoeff*nocoeff)
           workt = flegendre(t_nrm[*], nocoeff)
           
           for i=0,nocoeff-1 do begin
               for j=0,nycoeff-1 do begin
                   work2d[*,j*nocoeff+i] = worky[*, j] * workt[*,i]
               endfor
           endfor
           
           wv = dblarr(npix)
;           wv[*] = work2d # out_str.res
           wv[*] = work2d # out_str.res / tsub

;           stop
           
;           mn = min(10^wv, max=mx)
;           plot, all_pix, 10^wv, color=clr.black, $
           mn = min(wv, max=mx)
           plot, all_pix, wv, color=clr.black, $
             background=clr.white, charsize=1.5, yrange=[mn, mx], $
             xrange=[0.,sz[1]], xstyle=1, ystyle=1, xtitle='Row', $
             ytitle='Wavelength ('+STRING(197B)+')', xmargin=[11,2], ymargin=[5,1], /nodata
           
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
                 + sres*1., $
                 psym=1, color=clr.black
               ;; Rej
               rej = where(invvar[pts] LE 0., nrej, ncomplement=nnorej, $
                           complement=norej)
               if nrej NE 0 then $
                 oplot, [(sv_lines[jj].pix[0:nlin-1])[rej]],  $
                        [(all_wv[pts]/t[pts] + sres*1.)[rej]], $
                        psym=2, color=clr.red
;                       [(10^wv_mod[pts] + sres*100.)[rej]], $
               ;; RMS
               if nnorej NE 0 then $
                 rms = sqrt( total( sres[norej]^2 ) / float(nnorej-1)) $
               else rms = 9.99
           endif
           
           ;; Label
;           if side EQ 2 then ylbl = mn + (mx-mn)*0.08*(findgen(5)+1) $
;           else ylbl = mx - (mx-mn)*0.08*(findgen(5)+1) 
           ylbl = mx - (mx-mn)*0.08*(findgen(5)+1) 
           
           xyouts, sz[1]*0.05, ylbl[0], 'Order = '+strtrim(ii,2), $
             color=clr.black, charsize=1.5
           dwv = abs(wv[0]-wv[npix-1])/float(all_pix[npix-1])
;           dwv = abs(10^wv[0]-10^wv[npix-1])/float(all_pix[npix-1])
           xyouts, sz[1]*0.05, ylbl[2], '!9Dl!X = '+string(dwv,format='(f6.4)'), $
             color=clr.black, charsize=1.5
           xyouts, sz[1]*0.05, ylbl[1], 'RMS(pix) = '+$
             string(rms/dwv,format='(f5.2)'), color=clr.black, charsize=1.5
           pixrms[jj] = rms/dwv
       endfor
           
       x_psclose
       !p.multi=[0,1,1]

       replace_title = '"' + '%%Title: '+qafil + ' ' +systime() + '"'
       ps_replacetitle, replace_title, qafil

       spawn, 'gzip -f '+qafil
   endif

   return, 'Success'

end

