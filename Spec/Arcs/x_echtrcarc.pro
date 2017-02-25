;+ 
; NAME:
; x_echtrcarc   
;     Version 1.1
;
; PURPOSE:
;    To trace the arc lines in each order (individually) and fit a
;    straight line to each one.  The following steps are taken:
;    1.  Scattered light is removed from the image
;    2.  All significant arc lines are identified (5 sigma)
;    3.  trace_crude is used to trace the lines 
;    4.  trace_crude is reapplied to only those lines which are
;    entirely in the order
;    5.  xy2traceset is used to fit a straight line to each arc line
;    6.  Only the good lines are saved for 2D fitting in a structure
;    which is written to disk
;
; CALLING SEQUENCE:
; x_echtrcarc, arc_fil, ordr_str, out_fil, SZCCD=, $
;                            INIO=, ALL_XSET=, /CLOBBER, $
;                            QAFIL=
;   
; INPUTS:
;  arc_fil  -- Name of arc file for tracing
;  ordr_str -- Order strucure describing the echelle footprint
;  out_fil  -- Name of FITS file describing the traces
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   INIO      - Initial order to trace (for debugging)
;   /CLOBBER  - Overwrite previous fits
;  SZCCD -- Dimensions of the ccd
;
; OPTIONAL OUTPUTS:
; ALL_XSET --  A structure containing the info
; QAFIL    --  Filename for QA output
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Apr-2003 Written by SB
;
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function x_echtrcarc, arc_fil, ordr_str, out_fil, SZCCD=szccd, $
                             INIO=inio, ALL_XSET=all_xset, CLOBBER=clobber, $
                             QAFIL=qafil
;                             XOFF=xoff
;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = x_echtrcarc(arc_fil, ordr_str, out_fil, ' + $
        'SLOPECOEFF=, INIO=, /CLOBBER, ALL_XSET=, XIMAGE=, QAFIL= ) [v1.0]'
      return, -1
  endif 

;  Optional Keywords
  if not keyword_set( INIO ) then inio = 0L
  if not keyword_set( SAT ) then sat = 50000.
  if not keyword_set( SZCCD ) then szccd = [2048L,4096L]
;  if not keyword_set( XOFF ) then xoff =  0.

  ;; Ordstr
  nordr = n_elements(ordr_str)
          
  ;; Arc file
  print, 'x_echtrcarc: Reading arc: ', arc_fil
  arc = xmrdfits(arc_fil, 0, head, /silent)
  arcivar = xmrdfits(arc_fil, 1, /silent)
  sz_arc = size(arc, /dimensions)
          
  ;; Process
;  if ARG_PRESENT(ximage) then ximage=transpose(arc*0.0)
          
  norders = n_elements(ordr_str)
  norig_col = (size(arc))[1]
          
  medrow = 9L
  frac = (findgen(medrow)+0.5)/medrow
  nxcoeff=2 
;  slopecoeff = fltarr(nxcoeff, norders)
          

  ;; Make all_xset
  tmp = { $
          func: 'legendre', $
          xmin: 0., $
          xmax: 0., $
          ngood: 0L, $
          xguess: fltarr(200), $
          coeff: dblarr(2,200) $
        }
  all_xset = replicate(tmp,norders)


  ;; QA
  tmp = { $
          nslit: 0L, $
          ntrc: 0L, $
          msk: bytarr(300L,300L), $
          x1: fltarr(300L,300L), $
          xfit: fltarr(300L, 300L) $
        }
  sv_qa = replicate(tmp, norders)
          
  for iord=inio,norders-1 do begin
      
      ;; Set residual limit
      if keyword_set(BLUE) then begin
          if ordr_str[iord].order GT 105 then res_limit = 1. $
          else res_limit = 0.25
      endif else res_limit = 0.25

      if ordr_str[iord].flg_anly EQ 0 then begin
         print, 'x_echtrcarc: Skipping order ', ordr_str[iord].order
         continue
      endif
                  
      ;; Trace
      print, 'Tracing order...', iord, ordr_str[iord].order, $
              format=('(a,i6, i6, $)')

      ordr_width = (ordr_str[iord].rhedg - ordr_str[iord].lhedg)
      min_width = min(ordr_width) ;/2 ; JXP kludge -- Remove!!
;      if iord EQ (norders-1) then stop
      
      ;; Offset
      ;; Calculate xoff
      if keyword_set(SHFTPRM) then begin
	stop
;          xoff = mike_shifti(shftprm, ORDRS=ordr_str[iord].order)
      endif else xoff = 0.
      print, xoff, format=('(f9.4, $)')

      if abs(xoff) GT 10 then begin
        print, 'Offset is larger than 10 pixels, this could be a problem'
        stop
      endif

      mncen = 0.5*(ordr_str[iord].lhedg + ordr_str[iord].rhedg) + $
        xoff
      tempcen =  (ordr_str[iord].lhedg + xoff) # (1-frac) + $
        (ordr_str[iord].rhedg + xoff) # frac
      
      radius = min_width/medrow/2.
;      print, 'Working on extract_boxcar...', format='(a,$)'
      rec_tmp = extract_boxcar(arc, tempcen, radius=radius)/(2.0*radius)
      bad_extract = where(finite(rec_tmp) EQ 0)
      if bad_extract[0] NE -1 then rec_tmp[bad_extract] = 0
;     
;   Let's take off the scattered background...
;
              
      ncol = (size(rec_tmp))[1]
      x = findgen(ncol)
      width = 3L*ncol/256L      ; breakpoint separation
      xbin = (float(szccd[1])/ncol)
;      print, 'Working on djs_iterstat...',format='(a,$)'
      in = where(rec_tmp GT radius AND rec_tmp LT (sz_arc[1]) - radius, nin)
      if nin GT 10 then begin 
        djs_iterstat, 1.0d*rec_tmp[in], sigma=sigma, sigrej=2
        if (sigma LE 0 OR finite(sigma) EQ 0) then $
                 sigma = median(abs(rec_tmp[in])+1)
       endif else begin
          print, 'No good points here?', ordr_str[iord].order, xoff
          stop
       endelse
;      print, 'Finished'
              
      ;; Deal with saturation 
      saturated = (convol(1.0*(rec_tmp GT sat), fltarr(xbin*4)+1) GT 0)
      invvar = (x*0 + 1/sigma^2) * (saturated EQ 0)
     
      bset = bspline_iterfit(x, djs_median(rec_tmp,2), invvar=invvar, $
                             everyn=width, upper=2, lower=2, /groupbadpix, $
                             maxrej=3, niter=20, $
                             yfit=background, /silent)
      
      work = rec_tmp - background # replicate(1,medrow)
      invtemp = invvar # replicate(1,medrow)
      
      check = where(invtemp GT 0,ncheck)
      if ncheck LT 10 then continue
      
      djs_iterstat, work[check], median=medwork, sigma=sigwork 
      
;      print, 'Working on trace_crude...',format='(a,$)'
   
      xpeak = find_npeaks(total(work,2), nfind=100, ypeak=ypeak, minsep=3) 
      significant = where(ypeak GT 5.0*sqrt(medrow)*sigwork,ns) 
     
      if ns EQ 0 then begin
         print, ' Could not find any significant peaks'
         continue
      endif
 
      xstart = trace_crude(work, invtemp, nave=1, xstart=xpeak[significant], $
                           xerr=xerr, maxshifte=xbin) 

      print, ',  finished with', ns, ' peaks', format=('(a, i6, a, $)')

;
;   Now pick out the region which contains the whole order
;
          
      mincol = long(min(ordr_str[iord].lhedg + xoff)) > 0
      maxcol = long(max(ordr_str[iord].rhedg + xoff)+1) < (norig_col-1)
      nn = maxcol - mincol +1
      xcol = (lindgen(nn) + mincol) ## replicate(1,ncol)
      ll = (ordr_str[iord].lhedg + xoff)  # replicate(1,nn)
      uu = (ordr_str[iord].rhedg + xoff) # replicate(1,nn)
      
      submask = (xcol GE ll AND xcol LE uu) 
      sub = (transpose(arc[mincol:maxcol,*]) - background $
             # replicate(1,nn)) * submask
      
      meancoeff = 0.5*(ordr_str[iord].lhc + ordr_str[iord].rhc)
      meancoeff[0] = meancoeff[0] - mincol
      
;      print, 'Working on trace_crude again...',format='(a,$)'
      
      xguess = (xstart[medrow/2,*])[*]
      xnorm = (2.0*xguess-ordr_str[iord].ymin)/ $
        (ordr_str[iord].ymax-ordr_str[iord].ymin) -1
      ystart = flegendre(xnorm[*], n_elements(meancoeff)) # meancoeff
      
      ;; Require ystart are positive and smaller than the image
      ;; JXP -- bug fix 10feb04
      gd = where(ystart GT 0. AND ystart LT (size(sub))[2])
      xguess= xguess[gd]
      ystart= ystart[gd]
      
      xarc = trace_crude(sub, submask, yset=yarc, xerr=xerr, radius=2, $
                    ystart=ystart, xstart=xguess, maxshifte=xbin) 
      
      x1 = xarc
      ;; Extra 'tweaking' with trace_fweight
      for i=0,9 do x1 = trace_fweight(sub, x1, yarc, $
                                      radius=2, invvar=submask, xerr=xerr)
      
;;
;;   First try to fit each arcline trace individually w/ xy2traceset
;;
      xivar = 1.0/xerr^2 * (xerr LT 90)
      slit_frac = 2*(yarc - mncen[x1] + mincol)/ $ 
                    ordr_width[x1]
      xy2traceset, yarc, x1, xset, ncoeff=2, yfit=xfit, invvar=xivar, $
        outmask=outmask, maxdev=0.4, /silent
      traceset2xy, xset, transpose(ystart), xpos
      djs_iterstat, xset.coeff[1,*], mask=o1, med=medslope

      if ns EQ 1 then o1 = 1
;;
;;   find good lines...
;;
      tx = total(xivar GT 0,1)
      averesidual = sqrt(total((x1-xfit)^2 *(xivar GT 0),1)/(tx + (tx EQ 0)))
      goodcen = total(xerr LT 90,1)
 
      good_lines = where((averesidual LT res_limit) AND o1  AND tx NE 0 $
                          AND goodcen GT 0.75*min_width AND $
                           sub[xguess,ystart] LT sat, ngood)
      print, ', ', ngood, ' good', format='(a, i3, a)'
     
      if ngood EQ 0 then continue

      x_frac = (2.0*xfit[*,good_lines]-ordr_str[iord].ymin)/ $
        (ordr_str[iord].ymax-ordr_str[iord].ymin) -1

      x_basis = flegendre(x_frac[*], n_elements(meancoeff))
      l = x_basis # ordr_str[iord].lhc
      r = x_basis # ordr_str[iord].rhc


      ycen = 1.0*yarc[*,good_lines]
      ycen[*] = ycen -  (l+r)/2. + mincol
      limit = min_width/2. + 1

      xy2traceset, ycen, x1[*,good_lines], xset, ncoeff=2, $
            yfit=xfit_final, invvar=xivar[*,good_lines], $
            outmask=outmask, maxdev=0.2, $
            xmin=-1.0*limit, xmax=1.0*limit, /silent

      ;; Save
      all_xset[iord].xmin = xset.xmin
      all_xset[iord].xmax = xset.xmax
      all_xset[iord].xguess[0:ngood-1] = xguess[good_lines]
      all_xset[iord].coeff[0,0:ngood-1] = xset.coeff[0,*]
      all_xset[iord].coeff[1,0:ngood-1] = xset.coeff[1,*]
      all_xset[iord].ngood = ngood

      ;; Save QA
      szqa = size(x1, /dimensions)
      sv_qa[iord].nslit = szqa[0]
      sv_qa[iord].ntrc = ngood
      sv_qa[iord].msk[0:szqa[0]-1,0:ngood-1] = outmask
      sv_qa[iord].x1[0:szqa[0]-1,0:ngood-1] = x1[*,good_lines]
      sv_qa[iord].xfit[0:szqa[0]-1,0:ngood-1] = xfit_final

  endfor

  ;; Output
  print, 'x_echtrcarc: Writing '+out_fil
  mwrfits, all_xset, out_fil, /create
  mwrfits, sv_qa, out_fil
  spawn, 'gzip -f '+out_fil

  ;; QA
  if keyword_set(QAFIL) then begin
      x_psopen, qafil, /portrait
      clr = getcolor(/load)
      
      ;; Big loop
      for kk=-1L,999 do begin
          if kk EQ -1 then begin
              i1 = 0L
              i2 = norders-1
          endif else begin
              i1 = kk*5L
              if i1 GT norders-1 then break
              i2 = ((kk+1)*5L-1) < (norders-1)
          endelse
          if not keyword_set(BCKWD) then off = 0.5 else off = -0.5
          plot, [0.], [0.], /nodata, background=clr.white, charsize=1.5, $
            xrange=[ordr_str[i1].order-off, ordr_str[i2].order+off], thick=5, $
            yrange=[0.,sz_arc[1]], xstyle=1, ystyle=1, xtitle='Order', $
            ytitle='ypix', xmargin=[11,2], ymargin=[5,2]
          
          xyouts, 0.5, 0.97, 'QA Tracearc'+qafil+' (20x RES)', color=clr.black, $
            /normal, alignment=0.5
          
          for jj=i1,i2 do begin
              ;; x
              mxgd = 0L
              for ii=0L,sv_qa[jj].ntrc-1 do begin
                  good = where(sv_qa[jj].msk[0:sv_qa[jj].nslit-1,ii] EQ 1B, ngood)
                  if ngood NE 0 then begin
                      npt = good[ngood-1] - good[0] + 1
                      MXGD = MXGD > npt
                  endif 
              endfor
              
              nrm0 = 0.5 * mxgd
              nrm1 = mxgd
              !p.thick = 1
              for ii=0L,sv_qa[jj].ntrc-1 do begin
                  good = where(sv_qa[jj].msk[0:sv_qa[jj].nslit-1,ii] EQ 1B, ngood, $
                               complement=bad, ncomplement=nbad)
                  if ngood NE 0 then begin
                      xval = ordr_str[jj].order + (float(good)-nrm0-good[0])/nrm1
                      ;; Fit
                      oplot, [xval], [sv_qa[jj].xfit[good,ii]], color=clr.blue
                      ;;  Points
                      res = sv_qa[jj].x1[*,ii] - sv_qa[jj].xfit[*,ii]
                      oplot, [xval], [sv_qa[jj].xfit[good,ii] + 20*res[good]], $
                        psym=1, color=clr.black, symsize=0.5
                      ;; Rejected
                      if nbad NE 0 then begin
                          xbad = ordr_str[jj].order + (float(bad)-nrm0-good[0])/nrm1
                          plt = where(bad GT good[0] AND bad LT good[ngood-1], nplt)
                          if nplt NE 0 then $
                            oplot, [xbad[plt]], $
                            [sv_qa[jj].xfit[bad[plt],ii]+20*res[bad[plt]]], $
                            psym=1, color=clr.red, symsize=0.5
                      endif
                  endif
              endfor
          endfor
      endfor
      x_psclose
      !p.multi=[0,1,1]
      spawn, 'gzip -f '+qafil
  endif

  return, out_fil
end

