;+ 
; NAME:
; mike_arcalign
;     Version 1.1
;
; PURPOSE:
;    Finds the offset between the center of arc lines and the trace to
;    the trace flats.  The program is largely based on the routines
;    called in mike_slitflat.
;
; CALLING SEQUENCE:
;   
;  mike_arcalign, mike, setup, side, [/chk, TEMPL=]
;
; INPUTS:
;   mike   -  MIKE structure
;   setup  -  Setup identifier 
;   side   -  Blue (1), or Red (2) 
;
; RETURNS:
;
; OUTPUTS:
;   Updates the arc_xyoff tag in the mike structure
;
; OPTIONAL KEYWORDS:
;  /CHK  - Show the profiles and fits order by order
;  TEMPL= -- Template arc file (assumed to be the first good one of
;            the night)
;  /CLOBBER -- Overwrite any previous solution
;  NSMPL=  -- Number of orders to sample (x2)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_arcalign, mike, 1, 1, /CHK
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Jun-2004 Written by JXP (based on mike_slitflat by SB)
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
function mike_arcalign_work, arc, arci, ordri, chk=chk, NSMPL=nsmpl, $
                             ORGOFF=orgoff, FITPRM=fitprm, FLG=flg

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'rslt = mike_arcalign_work( arc, arci, ordri, /CHK, FITPRM=, /FLG, NSMPL= ) [v1.1]'
      return, -1
  endif 

  if NOT keyword_set(slit_samp) then slit_samp = 0.06
  if NOT keyword_set(NSMPL) then nsmpl = 8L
  if NOT keyword_set(ORGOFF) then orgoff = 0.

;    Colors
  clr = getcolor(/load)

  ;; Save
  ordr_str = ordri

  ;; Run arc align for a quick offset
  if not keyword_set( FLG ) and not keyword_set( ORGOFF ) then begin
      dxdy = mike_arcalign_work(arc, arci, ordr_str, CHK=chk, NSMPL=1L, /FLG)
      orgoff = dxdy[0]
  endif else dxdy = [0., 0.]

  ;; Apply the offset and then perform a more extensive calculation
  ordr_str.lhedg = ordr_str.lhedg + dxdy[0]
  ordr_str.rhedg = ordr_str.rhedg + dxdy[0]

  ;; Check for arc solution
  if ordr_str[0].arc_m[0] EQ 0. then begin
      print, 'mike_arcalign_work: Run mike_setarcm first!'
      return, -1
  endif
     
  sz        = size(arc, /dimen)
  n_orders = n_elements(ordr_str)

;
;    Mapping out orders and gaps, this is kind of slow
;
;     print, 'Calling mike_ordermask...', format='(a,$)'
     maskimage = mike_ordermask(sz[0], sz[1], ordr_str, trim=0.5)

;     print, 'Done.'
     print, ' i order             indx    frac     col      row      slope'

     ;; Offset the order structure

     ;; Zero out top and bottom
     lowb = round(1*sz[1]/4.)
     upb =  round(3*sz[1]/4.)

     arci[*,upb:*] = 0.
     arci[*,0:lowb-1] = 0.
     maskimage[*,upb:*] = -999.
     maskimage[*,0:lowb-1] = -999.


     ncol= (size(maskimage))[1]
     nrow= n_elements(ordr_str[0].lhedg)
     nord = n_elements(ordr_str)


   ;; Loop on the Orders
   slit_cen = fltarr(nord)
   sv_dx = fltarr(nord)
   sv_dy = fltarr(nord)
   q0 = nord/2 - nsmpl
   q1 = nord/2 + nsmpl

   arc_spec = extract_arc(arc, ordr_str)
   arc_med  = median_row(arc_spec, 25)
   med0 = djs_median(arc_spec, 1)
   ycol = dindgen(sz[1])



;   for q=q1-1,q1 do begin

   for q=q0,q1 do begin

       od = ordr_str[q].order
       inorder = where((maskimage EQ od  OR $
                       maskimage EQ -1*od OR $
                       maskimage EQ -1*(od-1)) AND arci GT 0, nin)

      if nin LE nrow then begin
          print, '...not enough pixels in order'
          continue
      endif 

       scat_set = bspline_iterfit(ycol, arc_spec[*,q], nbkpt=1, $
           invvar=(arc_spec[*,q] - arc_med[*,q]) LT 2*med0[q], $
           maxrej=5, /groupbadpix, upper=med0[q]/2., lower=med0[q])


      print, q, ordr_str[q].order, format='(i2, i4, $)'
      ystart = 1.0d*(inorder / ncol)
      ordrcen =  (ordr_str[q].lhedg[ystart] + ordr_str[q].rhedg[ystart])/2.0
      xstart = 1.0d*(inorder mod ncol)  

      ;; Calculate wavelengths (and the slope at each x,y pair)
      ywave = x_qckwav(xstart - ordrcen, ystart, ordr_str[q].arc_m, $
                        arc_slope=arc_slope)
      slit_frac = frac_order(ordr_str[q], xstart, ywave)

      y = ((2*ystart - nrow) / nrow)


      N = n_elements(inorder)
      xsort = sort(slit_frac)
      arc_norm = arc[inorder] - bspline_valu(ywave, scat_set)
      arci_norm = arci[inorder]
     
;
;    good will NOT work on dual slit... we need to find a way to get
;    pixel
;
      good = where(abs(slit_frac) LT 0.7 AND $
                   xstart GT 4 AND xstart LT ncol-4 AND $
                   arci_norm GT 0,ngood)
;
      everyn = long(0.6*ngood/(max(ywave[good]) - min(ywave[good]) +1))
      ys = sort(ywave[good])
      nfit = long(ngood/10) - 1
      pick_these = lindgen(10, nfit) + 5
      ywave_compress = total(ywave[good[ys[pick_these]]],1)/10.
      arc_norm_compress = djs_median(arc_norm[good[ys[pick_these]]],1)
      linterp, ywave_compress, arc_norm_compress, ywave, arc_fit
 
;      bset = bspline_iterfit(ywave[good[ys]], arc_norm[good[ys]],$
;                      everyn=everyn, /groupbadpix, maxrej=1, $
;                      invvar=arci_norm[good[ys]], upper=10, lower=10, /silent)

;      if total(bset.coeff) EQ 0 then begin
;         print, "Bspline fit failed, skipping "
;         continue
;      endif

;      arc_fit  = bspline_valu(ywave, bset)

;
;  Now it's time for the slit profile
;

      profile = arc_norm / (arc_fit + (arc_fit EQ 0)) * (arc_fit GT 10)
      profile_ivar = arci_norm * arc_fit^2 * $
        (profile GT -0.2 AND profile LT 1.5 AND arc_fit GT 10)
      w = where(profile_ivar[xsort] GT 0,nw)
 
;      if nw GT 200 then begin
;;        profile_smooth = median(profile[xsort[w]],200)
;        profile_ivar[xsort[w]] = profile_ivar[xsort[w]] * $
;              (abs(profile[xsort[w]] - profile_smooth) LT 0.2)
;      endif

 
      outlier = where( profile_ivar[xsort] GT 0., noutlier)
      if noutlier LT 10 then begin
         print, "Not enough useful pixels in this order"
         continue
      endif

      totalivar = total(profile_ivar[xsort]) / noutlier
      slit_samp_use = slit_samp 
      sfxo = slit_frac[xsort[outlier]]

;      nfit = long(noutlier / slit_length[0]) - 1
;      profile_compress = djs_avsigclip( profile[xsort[outlier[lindgen(slit_length[0],nfit) + slit_length[0]/2]]],1)
;      slit_frac_compress = total( slit_frac[xsort[outlier[lindgen(slit_length[0],nfit) + slit_length[0]/2]]],1) / slit_length[0]
;      profile_ivar_compress = total( profile_ivar[xsort[outlier[lindgen(slit_length[0],nfit) + slit_length[0]/2]]],1) 

 
       profile_set = bspline_iterfit(sfxo, profile[xsort[outlier]], $
                                    invvar=profile_ivar[xsort[outlier]], $
                                    upper=10, lower=10, /relative, $
                                    /groupbadpix, maxrej=10, bkpt=bkpts0, $
                                    bkspace=slit_samp_use, yfit=profile_fit, $
                                    outmask=profile_mask, /silent) 

      ;; Renormalize center
      cen = where(profile_fit GT 0.5 AND abs(sfxo) LT 0.7, ncen)

      if ncen LT 10 then begin
         print, "Not enough useful pixels in this order"
         continue
      endif

      fitcen = x_setfitstrct(NORD=1, /flgrej, niter=4, hsig=2.5, lsig=2.5)
      fit = x_fitrej(sfxo[cen], profile_fit[cen], fitstr=fitcen)
      all_scl = x_calcfit( sfxo, FITSTR=fitcen)
      profile_fit = profile_fit / all_scl
      ;; Center
      lhs = where(profile_fit LT 0.4 AND sfxo LT -0.5,nlhs)
      if nlhs EQ 0 then svl = min(sfxo) else svl = sfxo[lhs[nlhs-1]]
      rhs = where(profile_fit LT 0.4 AND sfxo GT 0.5, nrhs)
      if nrhs EQ 0 then svr=max(sfxo) else svr = sfxo[rhs[0]]
      slit_cen[q] = (svl+svr)/2.
      
      if keyword_set( CHK ) then begin
          plot, sfxo, profile_fit, yr=[0.,1.2], xr=[-1.4,1.4], $
             /xs, xtitle='Fractional slit position (from +/- 1 at edges)', $
             ytitle='Slit response'
          oplot, [svl,svl], [-9e9,9e9], color=clr.red
          oplot, [svr,svr], [-9e9,9e9], color=clr.red
          oplot, [slit_cen[q],slit_cen[q]], [-9e9,9e9], color=clr.blue
      endif

      ;; Save dx, dy
      ycen = sz[1]/2
      slen = ordr_str[q].rhedg[ycen] - ordr_str[q].lhedg[ycen]
      theta = atan(ordr_str[q].arc_m[ycen])
      
      sv_dx[q] = slen * (slit_cen[q]/2.) * cos(theta)
      sv_dy[q] = slen * (slit_cen[q]/2.) * sin(theta)

      print, ': Slit center', q, slit_cen[q], sv_dx[q], sv_dy[q], $
            (*fitcen.ffit)[1], format='(a, i6, 4(f9.4))'
  endfor

  ;; Order separation
  if NSMPL GT 2L then begin
      qval = q0 + lindgen(q1-q0+1)
;      ocen = (ordr_str.lhedg[ycen] + ordr_str.rhedg[ycen]) / 2.
;      osep = (shift(ocen,-1) - shift(ocen,1))/2.
;  slen = ordr_str[qval].rhedg[ycen] - ordr_str[qval].lhedg[ycen]
;  stop
      fitstr = x_setfitstrct(NITER=3L, /flgrej, NORD=1, hsig=2., lsig=2.)
      fit = x_fitrej(ordr_str[qval].order, sv_dx[qval]+orgoff, fitstr=fitstr)
      fitprm = fit_convxptox(fitstr)
      
      if keyword_set( CHK ) then begin
          plot, ordr_str[qval].order, sv_dx[qval]+orgoff, ps=1, $
            xtitle='Echelle Diffraction Order', $
            ytitle='Trace shift in binned pixels', /yno

          oplot, ordr_str.order, poly(ordr_str.order, fitprm), color=clr.red
      endif
      ;; ladfit
;      fit2 = ladfit(ordr_str[qval].order, sv_dx[qval]+orgoff)
;      stop
  endif

  return, [ median(sv_dx[q0:q1]) + dxdy[0], median(sv_dy[q0:q1]) + dxdy[1]]
end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_arcalign, mike, setup, side, indx, chk=chk, CLOBBER=clobber, $
                   NSMPL=nsmpl

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mike_arcalign, mike, setup, side, [indx], /chk, /clobber, NSMPL= [v1.1]' 
      return
  endif

   if not keyword_set( chk  ) then chk=0 $
   else  window, 0, title='mike_arclign'

   if not keyword_set( TEMPL  ) then templ = 0L
;   if not keyword_set( INDX  ) then indx = 0L

   ;; Grab all of the arcs
   arcs = where(mike.type EQ 'ARC' AND mike.flg_anly NE 0 $
                AND mike.setup EQ setup $
                AND mike.side EQ side, narc)

   ;; Open the order structure
   ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
   
   ;; Check to see if there is an arc solution, if not get one
   if ordr_str[0].arc_m[0] EQ 0. then begin
       rawfil = mike[arcs[templ]].rootpth+mike[arcs[templ]].img_root
       rslt = mike_allarc_sngl(rawfil, setup, side, /NOIMG, /NOWAV, /CLOBBER)
       ;; Reread the order structure
       ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
   endif

   if n_elements(indx) NE 0 then arcs = arcs[indx]
   ;; Process as necessary
   nindx = n_elements(arcs)
   for kk=0L,nindx-1 do begin
       idx = kk
       if total(mike[arcs[idx]].arc_xyoff) NE 0. and $
         not keyword_set( CLOBBER ) then begin
           print, 'mike_arcalign: xyoff already set.  Use /CLOBBER to override'
           continue
       endif
       arc_fil = mike_getfil('arc_fil', subfil=mike[arcs[idx]].img_root, $
                             /name)
       if x_chkfil( arc_fil+'*' ) EQ 0 then begin
           rawfil = mike[arcs[idx]].rootpth+mike[arcs[idx]].img_root
           rslt = mike_procarc_sngl( rawfil, setup, side )
       endif
       
       ;; Grab the images
       arc = xmrdfits(arc_fil, 0, /silent)
       arci = xmrdfits(arc_fil, 1, /silent)
       
       if total(abs(ordr_str.arc_m)) EQ 0 then begin
           print, "Need to have Arc slopes in arc_m"
           print, "Run mike_tracearc and mike_fittrcarc"
           return
       endif
       
       ;; Run arc align for a quick offset
;       dxdy = mike_arcalign_work(arc, arci, ordr_str, CHK=chk, NSMPL=1L)

       ;; Apply the offset and then perform a more extensive calculation
;       temp_ordr = ordr_str
;       temp_ordr.lhedg = temp_ordr.lhedg + dxdy[0]
;       temp_ordr.rhedg = temp_ordr.rhedg + dxdy[0]
       dxdy = mike_arcalign_work(arc, arci, ordr_str, $
                                 CHK=chk, NSMPL=nsmpl, FITPRM=fitprm)

       ;; Save
;       mike[arcs[idx]].arc_xyoff = dxdy + dxdy_fin
       mike[arcs[idx]].arc_xyoff = fitprm
       
;       if abs(dxdy[0]) GT 2.0 then begin
;         print, 'mike_arcalign: One more round to be safe'
;         temp_ordr = ordr_str
;         temp_ordr.lhedg = temp_ordr.lhedg + dxdy[0]
;         temp_ordr.rhedg = temp_ordr.rhedg + dxdy[0]
;         dxdy2 = mike_arcalign_work(arc, arci, temp_ordr, CHK=chk, NSMPL=nsmpl)
;         mike[arcs[idx]].arc_xyoff = mike[arcs[idx]].arc_xyoff + dxdy2
;       endif
   

       obj_match = where(mike.arc_fil EQ arc_fil)

;       if obj_match[0] NE -1 then begin
;          print, 'Also resetting these exposures ' , mike[obj_match].img_root
;          mike[obj_match].arc_xyoff =  mike[arcs[idx]].arc_xyoff
;       endif

       print, 'mike_arcalign: arc_xyoff = ', mike[arcs[idx]].arc_xyoff
       print, 'mike_arcalign: Do -- IDL> xatv, '''+arc_fil+'''
       print, 'mike_arcalign: Do -- IDL> mike_chktrcflat, mike, '+strtrim(setup,2)+', ' $
         +strtrim(side,2)+', /NOSTOP, /FIT, xoff='+string(dxdy[0],format='(f6.3)')
       
   endfor

   ;; All done
   print, 'mike_arcalign: All done'
   return

end


