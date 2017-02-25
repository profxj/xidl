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
pro fig_wavedisp, arcfil, ordr_str, arc_info, red=red, nycoeff=nycoeff, $
                nocoeff=nocoeff, CLOBBER=clobber, SZ=sz, $
                DEBUG=debug, CHKRES=chkres, out_str=out_str, $
                QAFIL=qafil, YVAL=yval, ORIG=orig, SIGREJ=sigrej
                   
;
;  if  N_params() LT 3  then begin 
;      print,'Syntax - ' + $
;        'rslt = x_fit2darc( arcfil, ordr_str, arc_info, OUT_STR=) [v1.1]'
;      return, -1
;  endif

   colors = GetColor(/Load, Start=1)

   if keyword_set(red) then begin 
       if not keyword_set(ARCFIL) then arcfil = $
         getenv('MIKE_PAP')+'Arcs/Arc_mr0005.fits'
       
       if not keyword_set(OSTRFIL) then ostrfil = $
         getenv('MIKE_PAP')+'Flats/OStr_R_01.fits'
       
       if not keyword_set(ARC_FIT) then arc_fit = $
         getenv('MIKE_PAP')+'Arcs/Fits/Arc_mr0005_fit2D.fits'
   endif else begin
       if not keyword_set(ARCFIL) then arcfil = $
         getenv('MIKE_PAP')+'Arcs/Arc_mb0005.fits'
       
       if not keyword_set(OSTRFIL) then ostrfil = $
         getenv('MIKE_PAP')+'Flats/OStr_B_01.fits'
       
       if not keyword_set(ARC_FIT) then arc_fit = $
         getenv('MIKE_PAP')+'Arcs/Fits/Arc_mb0005_fit2D.fits'
   endelse
   
   ordr_str = xmrdfits(ostrfil, 1, /silent)
   fit_str = xmrdfits(arc_fit,1)
       
   ;; OUTPUT Structure;   
   ;; out_str = { nrm: nrm, nrmt: nrmt, ny: nycoeff, no: nocoeff,  res: res }

  if not keyword_set(SZ) then begin
      head = xheadfits(arcfil, /silent)
      sz = lonarr(2)
      sz[0] = sxpar(head,'NAXIS1')
      sz[1] = sxpar(head,'NAXIS2')
  endif
  
  ;; NORMALIZE PIX
  all_pix = dindgen(sz[1])
  npix = n_elements(all_pix)
  pix_nrm = 2. * (all_pix - fit_str.nrm[0])/fit_str.nrm[1]

  worky = flegendre(pix_nrm[*], fit_str.ny)
       
  nordr = n_elements(ordr_str)
;  wvvec = dindgen(npix)
  wv = dblarr(npix,nordr)
  
  for jj=0L,nordr-1 do begin
      ;; NORMALIZE ORDER
      ii = ordr_str[jj].order
      tsub = replicate(float(ii), npix)
      t_nrm = 2. * (tsub - fit_str.nrmt[0])/fit_str.nrmt[1]
      
      ;; work2d and wv
      work2d = dblarr(npix,fit_str.ny*fit_str.no)
      workt = flegendre(t_nrm[*], fit_str.no)
       
      for i=0,fit_str.no-1 do begin
          for j=0,fit_str.ny-1 do begin
              work2d[*,j*fit_str.no+i] = worky[*, j] * workt[*,i]
          endfor
      endfor
           
;           wv[*] = work2d # out_str.res 

      wv[*,jj] = work2d # fit_str.res / ordr_str[jj].order

  endfor

  disp = wv / (shift(wv,1,0) - wv )  ; per 2x binned pixel  ~8.2"/pix
  
  disp = disp / 1000.   ;;;;;;;;;;;***********

  disp[0,*] = disp[1,*]   
;   disp_asec = disp/4.1

  cpix = 1020L
  cent = wv[cpix,*]
  lo = cent + 0.5 *( wv[cpix,*] - shift(wv[cpix,*],0,1)  )
  hi = cent + 0.5 *( wv[cpix,*] - shift(wv[cpix,*],0,-1) )
  lopix = fix(lo)
  hipix = fix(hi)
  junk = 0.
  for i =0L,nordr-1 do begin
      junk = min ( abs(wv[*,i]- lo[i]), ind)
      if keyword_set(red) then lopix[i] = fix(ind) else hipix[i] = fix(ind)
      junk = min ( abs(wv[*,i]- hi[i]), ind)
      if keyword_set(red) then hipix[i] = fix(ind) else lopix[i] = fix(ind)
  endfor
  
  set_plot, 'x'
  !p.background = colors.white
  !p.charsize = 1.5
  !x.thick = 1.5
  !y.thick = 1.5
  !x.margin = [10,10]

   if keyword_set(red) then begin
       x_psopen , 'wavedisp_red.ps', /maxs
       colors = GetColor(/Load, Start=1)
       color_side = colors.red
       hipix[nordr-1]=hipix[nordr-2]+2
       lopix[0]=lopix[1]-4
   endif else begin
       x_psopen , 'wavedisp_blue.ps', /maxs
       colors = GetColor(/Load, Start=1)
       color_side = colors.blue
       lopix[nordr-1]=lopix[nordr-2]+2
       hipix[0]=hipix[1]-4
   endelse
   
   i = 0L
   yran= ( max(disp)- min(disp) ) *0.05
   xran= ( max(wv)- min(wv) ) *0.05
   plot, wv[*,i], disp[*,i]*0 $
         , ystyle=9, xstyle=1, linesty=0 $
         , xrange=[min(wv)-xran, max(wv)+xran]  $
         , yrange=[min(disp)-yran, max(disp)+yran]  $
         , color=colors.black, charsize=1.5, thick=4 $
         , ytitle= 'Resolution (!9l/Dl!X) per pixel / 1000 '$
         , xtitle= 'Wavelength (Angstroms)' 
   axis, yaxis=1 $
         , yrange=3.e5/([min(disp)-yran, max(disp)+yran]*1000.) $
         , ytitle = 'km/s per pixel'
   for i=0L, nordr-1 do $
     oplot, wv[*,i], disp[*,i], color=color_side, thick=4,linest=1
   for i=0L,nordr-1 do $
     oplot, wv[lopix[i]:hipix[i],i], disp[lopix[i]:hipix[i],i], $
            color=colors.black, thick=4
   
   x_psclose
   ;if keyword_set(red) then spawn, 'gv -seascape wavedisp_red.ps' $
   ;  else spawn, 'gv -seascape wavedisp_blue.ps' 

   return
   
end

