;+ 
; NAME:
; mike_fittflat   
;     Version 1.1
;
; PURPOSE:
;  To create a 2D solution which describes the order curvature.  This
;  solution is derived from the individual traces created by
;  mike_trcflat and saved within the Trace structure.  The 2D fitting
;  algorithm is a simple least-squares algorithm.  The code then
;  attempts to extrapolate the solution for orders which are partially
;  on the CCD.  The code then makes a guess for the physical order
;  number which is not particularly accurate right now.  Finally, the
;  order structure (a key input for the MIKE pipeline) is written to
;  disk.
;
; CALLING SEQUENCE:
;   
;  mike_fittflat, mike, setup, [side], /DEBUG, INNY=, INNT=, LHEDGE=,
;                  /CLOBBER
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   INNY   -- Number of coefficients for fitting in vertical
;             direction.  (Default: 7)
;   INNT   -- Number of coefficients for fitting in horizontal
;             direction.  (Default: 6)
;   /DEBUG -- Turn debugging on
;   LHEDG  -- Used primarily for debugging
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_fittflat, mike, 1, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro mike_fittflat, mike, setup, side, DEBUG=debug, INNY=inny, INNT=innt, $
                   LHEDGE=lhedg, CLOBBER=clobber

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_fittflat, mike, setup, [side], /CLOBBER, /DEBUG, INNY=, ' + $
        'INNT=  [v1.1]'
      return
  endif

; Optional keywords
   if not keyword_set( SIDE ) then side = [1L,2L]

;  Loop on side
   for ii=0L,n_elements(side)-1 do begin
       qq = side[ii]
       ;; SIDE
       if qq EQ 1 then begin
           print, 'mike_fittflat: Tracing BLUE trace flat'
           if NOT keyword_set(INNY) then nycoeff = 7 else nycoeff=inny[ii]
           if NOT keyword_set(INNT) then ntcoeff = 6 else ntcoeff=innt[ii]
       endif else begin
           print, 'mike_fittflat: Tracing RED trace flat'
           if NOT keyword_set(INNY) then nycoeff = 7 else nycoeff=inny[ii]
           if NOT keyword_set(INNT) then ntcoeff = 6 else ntcoeff=innt[ii]
       endelse

       ;; Order structure
       ordr_fil = mike_getfil('ordr_str', setup, SIDE=qq, CHKFIL=chkf, /name)
       if CHKF NE 0 AND $
         not keyword_set( CLOBBER ) then begin
           print, 'mike_trcflat: Order structure exists. Moving on!', ordr_fil
           continue
       endif

       ;; Trace
       trc_str = mike_getfil('tflat_str', setup, SIDE=qq, FIL_NM=trc_fil, $
                             CHKFIL=chkf)
       if CHKF EQ 0 then begin
           print, 'mike_trcflat: Trace doesnt exist. ' + $
             'Run mike_trcflat first!', trc_fil
           continue
       endif

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; BEGIN with LH edges
       print, 'mike_fittflat: Fitting LH edges'
       neg = where(trc_str.flg EQ 1,nneg)
       xcen = trc_str.xfit[*,neg]

       ntrace = (size(xcen))[2]
       npix   = (size(xcen))[1]
       xerr = trc_str.xerr[*,neg] 

       nleft = ntrace
       
       ;; Assume xstrt are already sorted
       lnrm0 = (ntrace-1.)
       smallt = (2. * findgen(ntrace)/lnrm0) - 1
       ;; JXP BUG
       lbasis = flegendre(smallt, ntcoeff)
;       lbasis = fchebyshev(smallt, ntcoeff)
;       rbasis = fchebyshev(smallt, ntcoeff)

       trace2d, 1.0d*xcen, xerr, lhedg, mask=mask_lhs, $
            nycoeff=nycoeff, ntcoeff=ntcoeff, res=lh_res, sigrej=5

       ;; Get RMS
       gd_xcen = where(xerr LT 0.2, ngd)
       rms = sqrt( total( (lhedg[gd_xcen]-xcen[gd_xcen])^2 ) / float(ngd))
       ;; DEBUG
       if keyword_set( DEBUG ) then begin
           x_splot, lhedg[gd_xcen]-xcen[gd_xcen], /block
           stop
       endif
       print, 'mike_fittflat: LH edge RMS = ', rms
       lhrms = rms
       
   ;; Should probably add an iteration of rejection here

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; RH edges

       print, 'mike_fittflat: Fitting RH edges'
       pos = where(trc_str.flg EQ 2, npos)
       xcen = trc_str.xfit[*,pos]
       xerr = trc_str.xerr[*,pos]
       
       ntrace = (size(xcen))[2]
       npix   = (size(xcen))[1]
       
       ;; Setup basis vectors
       invvar = 1.0/xerr^2 * (xerr LT 90.0)
       y = ((2*findgen(npix) - npix) / npix) # replicate(1,ntrace)
       nrm0 = (ntrace-1.)/2.
       nrm1 = ntrace-1.
       smallt = 2. * (findgen(ntrace) - nrm0)/nrm1 
       rbasis = flegendre(smallt, ntcoeff)
;       rbasis = fchebyshev(smallt, ntcoeff)

       trace2d, 1.0d*xcen, xerr, rhedg, mask=mask_rhs, $
          nycoeff=nycoeff, ntcoeff=ntcoeff, res=rh_res, sigrej=5

       ;; Get RMS
       gd_xcen = where(xerr LT 0.2, ngd)
       rms = sqrt( total( (rhedg[gd_xcen]-xcen[gd_xcen])^2 ) / float(ngd))
       print, 'mike_fittflat: RH edge RMS = ', rms
       if keyword_set( DEBUG ) then begin
           x_splot, rhedg[gd_xcen]-xcen[gd_xcen], /block
           stop
       endif
       rhrms = rms

   ;; Should probably add an iteration of rejection here

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  ACCOUNTING

       ;; Find RH values for last order if necessary 
       ntrc = (size(trc_strxcen))[2]
       if max(trc_str.xfit[trc_str.smshrow,pos]) LT $
         trc_str.xfit[trc_str.smshrow,neg[nneg-1]] then begin
           if qq EQ 2 then begin
               print, 'mike_fittflat: Extrapolating RH edge of last order'
               ntrace = npos + 1
               ;; Set basis vectors
               y = ((2*dindgen(npix) - npix) / npix) # replicate(1,ntrace)
               smallt = 2. * (dindgen(ntrace) - nrm0)/nrm1 
               t = smallt ## replicate(1,npix)
               ;; now the functions
               work2d = fltarr(npix*ntrace,nycoeff*ntcoeff)
               ;; JXP Bug fix  22Feb05
               worky = flegendre(y[*], nycoeff)
               workt = flegendre(t[*], ntcoeff)
;               worky = fchebyshev(y[*], nycoeff)
;               workt = fchebyshev(t[*], ntcoeff)
               for i=0,ntcoeff-1 do begin
                   for j=0,nycoeff-1 do begin
                       work2d[*,j*ntcoeff+i] = worky[*, j] * workt[*,i]
                   endfor
               endfor
               ;; Finish
               rhedg = fltarr(npix,ntrace)
               rhedg[*] = work2d # rh_res
           endif else begin  ; Dont bother with blue side
               print, 'mike_fittflat: Partial blue order ignored'
           endelse
       endif

       ;; Find LH of first order if necessary
       if min(trc_str.xfit[trc_str.smshrow,pos]) LT $
         trc_str.xfit[trc_str.smshrow,neg[0]] then stop
;           print, 'mike_fittflat: Extrapolating LH edge of first order'
;           ntrace = nneg+1
;           ;; Set basis vectors
;           y = ((2*dindgen(npix) - npix) / npix) # replicate(1,ntrace)
;           smallt = 2. * (dindgen(ntrace) - nrm0)/nrm1 
;           t = smallt ## replicate(1,npix)
;           ;; now the functions
;           work2d = fltarr(npix*ntrace,nycoeff*ntcoeff)
;;           worky = fchebyshev(y[*], nycoeff)
;;           workt = fchebyshev(t[*], ntcoeff)
;           worky = flegendre(y[*], nycoeff)
;           workt = flegendre(t[*], ntcoeff)
;           for i=0,ntcoeff-1 do begin
;               for j=0,nycoeff-1 do begin
;                   work2d[*,j*ntcoeff+i] = worky[*, j] * workt[*,i]
;               endfor
;           endfor
;           ;; Finish
;           lhedg = fltarr(npix,nneg+1)
;           lhedg[*] = work2d # lh_res
;       endif

       ;;  NOTE:  IGNORING Upper RH corner of chip for now!!

       tmp = { $
               order: 0L, $     ;  Physical order #
               flg_anly: 0, $   ;  1=ok, 0=ng
               xcen: 0., $      ;  Center of order at ycen
               ycen: 0., $      ;  Smash row
               ymin: 0L, $      ;  Minimum used in y direction (-> -1)
               ymax: npix, $    ;  Maximum used in y direction (-> +1)
               lcen: 0.d, $     ;  Wavelength correspondoning to xcen,ycen
               lhedg: fltarr(npix), $ ; LH edge
               rhedg: fltarr(npix),  $ ; RH edge
               arc_m: fltarr(npix), $
               profile0 : fltarr(251),$ ; Slit profile
               profile1 : fltarr(251),$ ; Slit profile
               lhc : fltarr(nycoeff),$ ; LH coefficients
               rhc : fltarr(nycoeff) $ ; RH coefficients
             }
       
       smrow = trc_str.smshrow
       
       ord_str = replicate(tmp, nleft)
       ord_str.ycen = float(smrow)

       buffer = trc_str.swidth/3.      
       ;; Loop on LH edges
       for q=0L,nleft-1 do begin
           ;; Just use dummy number for now
           if qq EQ 1 then ord_str[q].order = trc_str.orderstrt+q $
           else ord_str[q].order = trc_str.orderstrt-q
           
           ;; Identify RH edge
           a = where( abs(rhedg[smrow,*]-lhedg[smrow,q]-$
                         trc_str.swidth) LT buffer $
                      AND rhedg[smrow,*] GT lhedg[smrow,q],na)
           if na NE 1 then begin
             if na GT 1 then begin
               print, 'Mutliple matches to RH edge'
               stop
             endif

             if na EQ 0 then begin
               if q EQ nleft-1 then begin
                 print, 'mike_fittflat: No match to last LH edge, dropping orders by 1'
                 nleft = nleft-1
                 ord_str = ord_str[0:nleft-1]
                 break
               endif else begin
                  print, 'mike_fittflat: No match and not even to the edge of CCD'
                  print, 'mike_fittflat: If your trace flat has very low counts to the blue, you may continue..'
                  print, 'mike_fittflat: Note, the code will not extract down there!'
                  stop
                  nleft = q
                  ord_str = ord_str[0:nleft-1]
                  break
               endelse
             endif 
           endif
            
           ;; xcen
           ord_str[q].xcen = (lhedg[smrow,q]+rhedg[smrow,a[0]])/2.
           
           ;; Fill in edges
           ord_str[q].lhedg = lhedg[*,q]
           ord_str[q].rhedg = rhedg[*,a[0]]
           
           ;; Add in smooth y legendre coefficients
           
           ord_str[q].lhc = reform(lh_res,ntcoeff,nycoeff) ## lbasis[q,*]
           ord_str[q].rhc = reform(rh_res,ntcoeff,nycoeff) ## rbasis[a[0],*]
       endfor

       ;; Assign order
       ord_str.order = lindgen(n_elements(ord_str)) + 1000
 
       ;; Write to fits
       print, 'mike_fittflat: Writing order structure ', ordr_fil
       mwrfits, ord_str, ordr_fil, /create

       ;; QA
;       restore, 'tmp.idl'
       idx = where(mike.setup EQ setup AND mike.side EQ qq AND mike.flg_anly NE 0)
       sz = lonarr(2)
       sz[0] = round(4096. / mike[idx[0]].rowbin)
       sz[1] = round(2048. / mike[idx[0]].colbin)

       qafil = mike_getfil('qa_fittflat', setup, SIDE=qq)
       x_psopen, qafil, /maxs
       clr = getcolor(/load)
       ;; LHS 
       !p.multi=[0,1,1]
       xcen = trc_str.xcen[*,neg]
       ycen = findgen(npix) # replicate(1.,nleft)
       gd = where(mask_lhs EQ 1B, complement=bad, ncomplement=nbad)

       plot, ycen[gd], 10*(xcen[gd]-lhedg[gd])+lhedg[gd], $
         psym=3, color=clr.black, thick=5, xstyle=1, ystyle=1, $
         background=clr.white, charsize=1.5, xrange=[0., sz[0]], yrange=[0.,sz[1]]
       !p.thick=1
       if nbad NE 0 then $
         oplot, ycen[bad], 10*(xcen[bad]-lhedg[bad])+lhedg[bad], $
         psym=3, color=clr.red
       for jj=0L,nleft-1 do begin
           oplot, findgen(sz[0]), lhedg[*,jj], color=clr.green
       endfor
       xyouts, 0.5, 0.96, 'LHS CEN (10x resid)'+'  RMS = ' $
         +string(lhrms,format='(f7.5)'), color=clr.blue, $
         charsize=1.5, /normal, alignment=0.5
       ;; FIT
       xcen = trc_str.xfit[*,neg]
       ycen = findgen(npix) # replicate(1.,nleft)
       gd = where(mask_lhs EQ 1B, complement=bad, ncomplement=nbad)

       plot, ycen[gd], 50*(xcen[gd]-lhedg[gd])+lhedg[gd], $
         psym=3, color=clr.black, xstyle=1, ystyle=1, $
         background=clr.white, charsize=1.5, xrange=[0., sz[0]], yrange=[0.,sz[1]]
       !p.thick=1
       if nbad NE 0 then $
         oplot, ycen[bad], 50*(xcen[bad]-lhedg[bad])+lhedg[bad], $
         psym=3, color=clr.red
       for jj=0L,nleft-1 do begin
           oplot, findgen(sz[0]), lhedg[*,jj], color=clr.green
       endfor
       xyouts, 0.5, 0.96, 'LHS FIT (50x resid)'+'  RMS = '+string(lhrms,format='(f7.5)'), color=clr.blue, charsize=1.5, /normal, alignment=0.5

       ;; RHS  (positive peaks)
       xcen = trc_str.xcen[*,pos]
       ycen = findgen(npix) # replicate(1.,npos)
       gd = where(mask_rhs EQ 1B, complement=bad, ncomplement=nbad)

       plot, ycen[gd], 10*(xcen[gd]-rhedg[gd])+rhedg[gd], $
         psym=3, color=clr.black, thick = 5, xstyle=1, ystyle=1, $
         background=clr.white, charsize=1.5, xrange=[0., sz[0]], yrange=[0.,sz[1]]
       if nbad NE 0 then $
         oplot, ycen[bad], 10*(xcen[bad]-rhedg[bad])+rhedg[bad], $
         psym=3, color=clr.red
       for jj=0L,npos-1 do begin
           oplot, findgen(sz[0]), rhedg[*,jj], color=clr.green
       endfor
       xyouts, 0.5, 0.96, 'RHS CEN (10x resid)'+'  RMS = '+string(rhrms,format='(f7.5)'), color=clr.blue, charsize=1.5, /normal, alignment=0.5
       ;; FIT
       xcen = trc_str.xfit[*,pos]
       ycen = findgen(npix) # replicate(1.,npos)
       gd = where(mask_lhs EQ 1B, complement=bad, ncomplement=nbad)

       plot, ycen[gd], 50*(xcen[gd]-rhedg[gd])+rhedg[gd], $
         psym=3, color=clr.black, thick=5, xstyle=1, ystyle=1, $
         background=clr.white, charsize=1.5, xrange=[0., sz[0]], yrange=[0.,sz[1]]
       !p.thick=1
       if nbad NE 0 then $
         oplot, ycen[bad], 50*(xcen[bad]-rhedg[bad])+rhedg[bad], $
         psym=3, color=clr.red
       for jj=0L,npos-1 do begin
           oplot, findgen(sz[0]), rhedg[*,jj], color=clr.green
       endfor
       xyouts, 0.5, 0.96, 'RHS FIT (50x resid)'+'  RMS = '+string(rhrms,format='(f7.5)'), color=clr.blue, $
         charsize=1.5, /normal, alignment=0.5


       ;; Close
       x_psclose
       !p.multi=[0,1,1]
       spawn, 'gzip -f '+qafil
       
   endfor

   ;; All done
   print, 'mike_fittflat: All done'
   return

end

