;+ 
; NAME:
; mike_fntobj   
;     Version 1.1
;
; PURPOSE:
;    Identify the object within the slit (single object) and trace it
;    along each order.  The code rectifies each order, collapses it
;    and looks for the object.  IF it finds it in >7 orders then it
;    peaks up on the object.  Otherwise, it guesses the object is at
;    the center of the order.  It then uses xy2tracset (trace_crude)
;    to trace the object along each order individually.  Finally, it
;    performs a PCA analysis on the trace_crude coefficients to create
;    a smoothed (pseudo-2D) solution for the trace.  The code outputs
;    an object structure which includes the trace and other info.
;
;
; CALLING SEQUENCE:
;   
;    mike_fntobj, mike, setup, obj_id, side, [exp], /STD, /CHK [v1.0]'
;
; INPUTS:
;   mike    -  MIKE structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   side    -  Blue (1) or Red (2) side
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;  Object structure (MIKEOBJSTRCT) including the trace and other information
;  regarding the object.  This structure is then filled up with the 1D
;  spectrum, etc.
;
; OPTIONAL KEYWORDS:
;  /STD     - Find object for standard star
;  /CHK     - Show diagnostics of the code performance
;  /NOCLOB  - Do not clobber existing files (default: clobber)
;  /DEBUG   - Debug
;  OBJAPER  - Set aperture to mask for sky subtraction (default: 20
;             unbinned pixels for obj, 26 for STD)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_fntobj, mike, 1L, 1L, /CHK, 
;
; PROCEDURES/FUNCTIONS CALLED:
;  mike_rectify
;  find_npeaks
;  mike_fweight
;  mike_fntobj_off
;  trace_crude
;  mike_basis
;
; REVISION HISTORY:
;   23-Sep-2003 Written by JXP (combined aspects of fndobj with trcobj)
;   10-Nov-2003 Changed fntobj_off, SMB
;-
;------------------------------------------------------------------------------

function mike_fntobj_off, ordr_str, img, qq, DEBUG=debug, xerr=xerr, $
                          CHK=chk

  ;; Rectify the image
  rect_image = mike_rectify(img, ordr_str[qq].lhedg, ordr_str[qq].rhedg)

  ;; Smooth
  sz = size(img,/dimen)
  vertical_median_smooth = 63*sz[1]/4096L
  half_smooth = (vertical_median_smooth+1)/2
  nf = sz[1] / half_smooth*half_smooth
  rect_rebin = djs_median(reform(rect_image[*,0:nf-1], $
              (size(rect_image))[1], half_smooth , sz[1]/half_smooth), 2)
  rect_smooth = median_row(transpose(rect_image), vertical_median_smooth,$
                   reflect=1)

  ;; Trace down the image
  down1= trace_crude(rect_rebin, xstart=(size(rect_image))[1]/2, $
                maxshift0=5,maxshifte=0.3)
  down = mike_fweight(rect_rebin, down1, niter=10)

  smash_image= total(rect_smooth,1) / (size(rect_image))[2] 
  smash_sz = n_elements(smash_image)
  center = find_npeaks(smash_image, nfind=1, width=3L, ypeak=peak, $
                 xerr=peak_err)
  center= median(down)
  if keyword_set( CHK ) then begin
      plot, smash_image, /yno, psym=10
      wait, 0.2
  endif
  slitfrac = center/(smash_sz-1)

  ;; Fraction
  frac = center - fix(center)
  value = convol(djs_median(smash_image,width=3,boundary='reflect'), $
                 [1.0-frac,frac],/edge_wrap)
  if keyword_set( CHK ) then djs_oplot, value, color='red'
  pixl  = fix(center) - 2 > 0
  pixmid = fix(center) + 1 < (smash_sz-1)
  pixr = fix(center) + 4 < (smash_sz-1)
  
  yl    = value[pixl]
  yr    = value[pixr]
  ymid  = value[pixmid]
  
  xerr = 2.0*sqrt(((yl> 0) +  (yr > 0) + 1.))/ ((2.0*ymid - yl -yr) > 0.1) / $
    sqrt((size(img,/dimen))[1]) 
  

  if keyword_set(DEBUG) then print, center, slitfrac, peak, xerr, peak_err
  xerr = xerr * (peak_err GT 0)
  return, slitfrac
  
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_fntobj, mike, setup, obj_id, side, exp, $
                   CHK=chk, NOCLOB=noclob, OBJAPER=objaper, $
                   DEBUG=debug, STD=std, NONLIN=nonlin, _EXTRA=extra

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'mike_fntobj, mike, setup, obj_id, side, [exp], /STD, /CHK, ' 
      print, '       POLY_NCOEFF=, OBJAPER= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  ;; QA
  if keyword_set(CHK) then window, 0, title='mike_fntobj'
  if setup LT 10 then fqa = 'QA/Obj0'+strtrim(setup,2) $
  else fqa = 'QA/Obj'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Find all relevant obj images
  trimtype = strtrim(mike.type,2)

  indx = where(mike.flg_anly NE 0 AND mike.setup EQ setup AND $
               mike.side EQ side[0]  AND $
               mike.obj_id EQ obj_id AND $
               (trimtype EQ 'OBJ' OR trimtype EQ 'STD'), nindx)
  if nindx EQ 0 then begin
      print, 'mike_fntobj: No images to find obj for!', obj_id, ' and side ',$
             side
      stop
      return
  endif
  
; APERTURES
  if not keyword_set( OBJAPER ) then begin
      if not keyword_set(STD) then objaper = [0.75, 0.75] $
      else objaper = [0.9, 0.9]
  endif


;  Read in order structure
  ordr_str = mike_getfil('ordr_str', setup, SIDE=side)
  sv_ordrstr = ordr_str
  nordr = n_elements(ordr_str)

;  Exposures
  if n_elements(exp) EQ 0 then exp = lindgen(nindx)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = mike[indx[exp[q]]].obj_fil

      ;; xyoff
      ordr_str = sv_ordrstr
      shftprm = mike[indx[exp[q]]].arc_xyoff
      shft = mike_shifti(shftprm, OSTR=ordr_str)
;      xyoff = mike[indx[exp[q]]].arc_xyoff[0]
;      ordr_str.lhedg = sv_ordrstr.lhedg + xyoff[0]
;      ordr_str.rhedg = sv_ordrstr.rhedg + xyoff[0]

      ;; Cen ref
      cen_ref = (ordr_str.lhedg + ordr_str.rhedg)/2. 
      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
          print, 'mike_fntobj: Using Obj structure -- ', objfil
          if x_chkfil(objfil+'*') EQ 0 then begin
              print, 'mike_fntobj: File ', objfil, $
                'doesnt exist! Turn off NOCLOB'
              return
          endif
          objstr = xmrdfits(objfil, 1, STRUCTYP='mikeobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          ;; objfil = 'Extract/Obj_'+mike[indx[exp[q]]].img_root
          objfil = mike_getfil('obj_fil', /name, $
              subfil='Extract/Obj_'+mike[indx[exp[q]]].img_root)

          mike[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { mikeobjstrct }
          objstr = replicate(tmp, nordr)
          objstr.UT = ' '
          objstr.field = ' '
          objstr.img_fil = ' '
          objstr.arc_fil = ' '
          objstr.exp = mike[indx[exp[q]]].exp
          objstr.field = mike[indx[exp[q]]].Obj
          nobj = 0L
      endelse
      objstr.flg_anly = 1

      ;; Read IMG+VAR
      imgfil = mike[indx[exp[q]]].img_final
      a = findfile(imgfil+'*', count=na)
      if na EQ 0 then begin
          print, 'mike_fntobj: No Final file ', imgfil
          print, '... try running mike_proc'
          continue
      endif
      objstr.img_fil = imgfil
      print, 'mike_fntobj: Opening files... ', imgfil
      img = xmrdfits(imgfil, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Read ARC
      flg_arc = 0
      objstr.arc_fil = strtrim(mike[indx[exp[q]]].arc_img,2)
      if q EQ 0 then begin
          img_arc = xmrdfits(strtrim(mike[indx[exp[q]]].arc_img,2), /silent) 
      endif else begin
          if mike[indx[exp[q]]].arc_img NE mike[indx[exp[q-1]]].arc_img then $
            img_arc = xmrdfits(mike[indx[exp[q]]].arc_img, /silent) $
          else flg_arc=1
      endelse
            

      ;; QA
      qafil = mike_getfil('qa_fntobj', setup, SUBFIL=mike[indx[exp[q]]].img_root)

      x_fntobj, objstr, ordr_str, img, ivar, qafil, $
        CHK=chk, NOCLOB=noclob, OBJAPER=objaper, FCOEFF=fcoeff, $
        FULLTHR=fullthr, DEBUG=debug, STD=std, _EXTRA=extra
      

      ;; 
      ;; Write Obj structure
      print, 'mike_fntobj: Creating ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'mike_fntobj: All done! '
  return
end


;      ;; FIND OBJ
;      print, 'mike_fntobj: Finding the object...'
;      frac_array = fltarr(nordr)
;      xerr_array = fltarr(nordr)
;      for qq = 0,nordr-1 do begin & $
;       frac_array[qq] = mike_fntobj_off(ordr_str, img, qq, xerr=xerr_temp, $
;                                  debug=debug, CHK=chk) & $
;       xerr_array[qq] = xerr_temp & $
;      endfor
;      good = where(xerr_array GT 0 AND xerr_array LT 0.4, ngood)
;
;      if ngood GT 7 then begin
;        lad_coeff = ladfit(good, frac_array[good])
;        lad_fit = poly(findgen(nordr),lad_coeff)
;        ; frac_guess = (lad_fit + median(frac_array[good]))/2.
;        frac_guess = lad_fit 
;      endif else $
;        if ngood GE 1 then frac_guess=fltarr(nordr) + median(frac_array[good]) $
;        else               frac_guess=fltarr(nordr) + 0.5
;
;      ordr_width = (ordr_str.rhedg-ordr_str.lhedg)
;      trace_guess = ordr_width * (frac_guess##replicate(1,sz_img[1])) + $
;                    ordr_str.lhedg
;
;      if keyword_set(CHK) then begin
;          plot, frac_array, xtitle='Trace Number', ytitle='slit fraction', $
;            xmargin=[13,1], xstyle=1
;          if good[0] NE -1 then oplot, [good], [frac_array[good]],ps=1 
;          oplot, frac_guess
;      endif
;      frac_good = good          ; Save for QA
;          
;  
;      ;; TRACE
;      print, 'mike_fntobj: Tracing...'
;      yset = findgen(sz_img[1]) # replicate(1,nordr)
;      radius = fix(median(ordr_width)/4.) +1
;      x2 = mike_fweight(img, trace_guess, yset, invvar=ivar, $
;                 niter=10, xerr=xerr, radius=radius, sig=0.7)
;
;      shift = total((x2-trace_guess)/xerr^2,1)/total(1./xerr^2,1)
;      xivar = (xerr LT 900)/(xerr^2 + 0.01^2)
;      ncoeff = 7
;      xy2traceset, yset, x2, tset, invvar=xivar, ncoeff=ncoeff, yfit=x2fit, $
;           outmask=outmask, upper=10, lower=10
;
;
;      ;  Find eigenvectors.....
;      print, 'MIKE_FNTOBJ:  PCA...'
;      x3fit = x_basis(tset, x2, outmask, ncoeff=ncoeff, eigenvec=eigenvec,$
;                        msktrc=msktrc, outstr=pca_out, NONLIN=nonlin, $
;                     POLY_NCOEFF=poly_ncoeff)
;
;      if keyword_set( DEBUG ) then begin
;          tmp = img
;          for qq=0L,nordr-1 do begin
;              rnd_trc = round(x3fit[*,qq]) 
;              plot_this = where(rnd_trc GE 0 AND rnd_trc LT sz_img[0])
;              if plot_this[0] NE -1 then $
;                tmp[rnd_trc[plot_this], plot_this] = -10000
;          endfor
;;          print, objstr[mm].xcen, objstr[mm].ycen
;          xatv, tmp, /block, /histeq
;          stop
;      endif
;          
;
;      if n_elements(x3fit) EQ 1 then begin
;        ;; time for plan B  use order center
;        print, 'MIKE_FNTOBJ: Using order centers as guide for trace'
;
;        xy2traceset, yset, cen_ref, ordrcen_set, ncoeff=ncoeff, $
;            yfit=x2fit, outmask=outmask
;
;        x3fit = mike_basis(ordrcen_set, x2, outmask, $
;                           ncoeff=ncoeff, eigenvec=eigenvec, $
;                           msktrc=msktrc, outstr=pca_out)
;
;      endif
;;
;     final tweak to traces...
;;
;      if x3fit[0] NE -1L then begin
; 
;          x4 = mike_fweight(img, x3fit, yset, invvar=ivar, $
;                            niter=5, xerr=xerr4, radius=radius/1.5, sig=0.7)
;          
;          x_ivar = 1./(xerr4 + (xerr4 EQ 0))^2 * (xerr4 GT 0 AND xerr4 LT 1.5) * $
;                   (x4 GT radius) * (x4 LT sz_img[0]-radius)
;          ; rebin into 30 sections x nord
;          nbin = long(sz_img[1] / 31)
;          lo = long(0.5*nbin)
;          x_ivar_check = total(reform(x_ivar[lo:lo+30*nbin-1,*],nbin, 30, nordr),1)
;          fullfit_fib = where(total(x_ivar_check GT 1000.,1) EQ 30, nfullfit)
;
;          print, 'There are ', nfullfit, " fibers which don't require PCA interpolation"
;          
;
;          if keyword_set( DEBUG ) then begin
;              tmp = img
;              for qq=0L,nordr-1 do begin
;                  rnd_trc = round(x4[*,qq]) 
;                  plot_this = where(rnd_trc GE 0 AND rnd_trc LT sz_img[0])
;                  if plot_this[0] NE -1 then $
;                    tmp[rnd_trc[plot_this], plot_this] = -10000
;              endfor
;;          print, objstr[mm].xcen, objstr[mm].ycen
;              xatv, tmp, /block, /histeq
;              stop
;          endif
;          
;          xy2traceset, yset, 1.0d*(x4-x3fit), tset4, invvar=(xerr4 LT 1.5), $
;            ncoeff=5, yfit=x4fit, outmask=outmask
;          xfinal_fit = x3fit 
;          addin = where(total(x4fit^2,1) LT 4*sz_img[1])
;          xfinal_fit[*,addin] =  x3fit[*,addin] + x4fit[*,addin]
;
;          if fullfit_fib[0] NE -1 then begin
;               xy2traceset, yset, x4, tset4full, $
;                  ncoeff=15, yfit=x4full_fit, outmask=outmask
;               xfinal_fit[*,fullfit_fib] = x4full_fit[*,fullfit_fib]
;          endif
;
;      endif else begin
;          xfinal_fit = x2fit
;      endelse
;
;
;      ;; Fill up the Object Structure
;      nordr = n_elements(ordr_str)
;      objstr[0:nordr-1].trace[0:sz_img[1]-1] = xfinal_fit
;      objstr[0:nordr-1].aper = objaper # replicate(1,nordr)
;      objstr[0:nordr-1].obj_id = 'a'
;      objstr[0:nordr-1].order = ordr_str.order
;      medx = sz_img[1]/2.
;      objstr[0:nordr-1].xcen = medx
;
;      if x3fit[0] EQ -1L then begin
;          ;; No basis
;          objstr[0:nordr-1].ycen = (x2fit[medx,*] - trace_guess[medx,*])[*]
;      endif else begin
;          objstr[0:nordr-1].ycen = (x3fit[medx,*] - trace_guess[medx,*])[*]
;      endelse
