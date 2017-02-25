;+ 
; NAME:
; uves_fntobj   
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
;    uves_fntobj, uves, setup, obj_id, side, [exp], /STD, /CHK [v1.0]'
;
; INPUTS:
;   uves    -  MIKE structure
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
;   uves_fntobj, uves, 1L, 1L, /CHK, 
;
; PROCEDURES/FUNCTIONS CALLED:
;  uves_rectify
;  find_npeaks
;  uves_fweight
;  uves_fntobj_off
;  trace_crude
;  uves_basis
;
; REVISION HISTORY:
;   23-Sep-2003 Written by JXP (combined aspects of fndobj with trcobj)
;   10-Nov-2003 Changed fntobj_off, SMB
;-
;------------------------------------------------------------------------------

function uves_fntobj_off, ordr_str, img, qq, DEBUG=debug, xerr=xerr, $
                          CHK=chk

  ;; Rectify the image
  rect_image = uves_rectify(img, ordr_str[qq].lhedg, ordr_str[qq].rhedg)

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
  down = x_fweight(rect_rebin, down1, niter=10)

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

pro uves_fntobj, uves, setup, obj_id, side, exp, $
                   CHK=chk, NOCLOB=noclob, OBJAPER=objaper, $
                   DEBUG=debug, STD=std, NONLIN=nonlin, _EXTRA=extra

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'uves_fntobj, uves, setup, obj_id, side, [exp], /STD, /CHK, ' 
      print, '       POLY_NCOEFF=, OBJAPER= [v1.1]'
      return
  endif 
  
;  Optional Keywords


;  Find all relevant obj images
  trimtype = strtrim(uves.type,2)

  indx = where(uves.flg_anly NE 0 AND uves.setup EQ setup AND $
               uves.side EQ side[0]  AND $
               uves.obj_id EQ obj_id AND $
               (trimtype EQ 'OBJ' OR trimtype EQ 'STD'), nindx)
  if nindx EQ 0 then begin
      print, 'uves_fntobj: No images to find obj for!', obj_id, ' and side ',$
             side
      stop
      return
  endif

  ;; QA
  wcen = round(uves[indx[0]].xdangl)
  cwcen = strtrim(wcen,2)

  if keyword_set(CHK) then window, 0, title='uves_fntobj'
  if setup LT 10 then c_s = '0'+strtrim(setup,2) $
  else c_s = strtrim(setup,2)
  fqa = 'QA/Obj'+cwcen+'_'+c_s
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa
  
; APERTURES
  if not keyword_set( OBJAPER ) then begin
      if not keyword_set(STD) then objaper = [0.75, 0.75] $
      else objaper = [0.9, 0.9]
  endif


;  Read in order structure
  ordr_str = uves_getfil('ordr_str', setup, WCEN=wcen)
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
      objfil = uves[indx[exp[q]]].obj_fil

      ;; xyoff
      ordr_str = sv_ordrstr
      shftprm = uves[indx[exp[q]]].arc_xyoff
;      shft = uves_shifti(shftprm, OSTR=ordr_str)
;      xyoff = uves[indx[exp[q]]].arc_xyoff[0]
;      ordr_str.lhedg = sv_ordrstr.lhedg + xyoff[0]
;      ordr_str.rhedg = sv_ordrstr.rhedg + xyoff[0]

      ;; Cen ref
      cen_ref = (ordr_str.lhedg + ordr_str.rhedg)/2. 
      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
          print, 'uves_fntobj: Using Obj structure -- ', objfil
          if x_chkfil(objfil+'*') EQ 0 then begin
              print, 'uves_fntobj: File ', objfil, $
                'doesnt exist! Turn off NOCLOB'
              return
          endif
          objstr = xmrdfits(objfil, 1, STRUCTYP='uvesobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          ;; objfil = 'Extract/Obj_'+uves[indx[exp[q]]].img_root
          objfil = uves_getfil('obj_fil', /name, objn=uves[indx[exp[q]]].img_root)

          uves[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { uvesobjstrct }
          objstr = replicate(tmp, nordr)
          objstr.UT = ' '
          objstr.field = ' '
          objstr.img_fil = ' '
          objstr.arc_fil = ' '
          objstr.exp = uves[indx[exp[q]]].exp
          objstr.field = uves[indx[exp[q]]].Obj
          nobj = 0L
      endelse
      objstr.flg_anly = 1

      ;; Read IMG+VAR
      imgfil = uves[indx[exp[q]]].img_final
      a = findfile(imgfil+'*', count=na)
      if na EQ 0 then begin
          print, 'uves_fntobj: No Final file ', imgfil
          print, '... try running uves_proc'
          continue
      endif
      objstr.img_fil = imgfil
      print, 'uves_fntobj: Opening files... ', imgfil
      img = xmrdfits(imgfil, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; Read ARC
      flg_arc = 0
      objstr.arc_fil = strtrim(uves[indx[exp[q]]].arc_img,2)
      if q EQ 0 then begin
          img_arc = xmrdfits(strtrim(uves[indx[exp[q]]].arc_img,2), /silent) 
      endif else begin
          if uves[indx[exp[q]]].arc_img NE uves[indx[exp[q-1]]].arc_img then $
            img_arc = xmrdfits(uves[indx[exp[q]]].arc_img, /silent) $
          else flg_arc=1
      endelse
            

      ;; QA
      qafil = uves_getfil('qa_fntobj', setup, WCEN=wcen, $
                          OBJN=uves[indx[exp[q]]].img_root)

      x_fntobj, objstr, ordr_str, img, ivar, qafil, $
        CHK=chk, NOCLOB=noclob, OBJAPER=objaper, FCOEFF=fcoeff, $
        FULLTHR=fullthr, DEBUG=debug, STD=std, _EXTRA=extra
      

      ;; 
      ;; Write Obj structure
      print, 'uves_fntobj: Creating ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'uves_fntobj: All done! '
  return
end


