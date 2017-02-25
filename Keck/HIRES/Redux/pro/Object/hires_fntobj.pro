;+ 
; NAME:
; hires_fntobj   
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
;  The main call is to x_fntobj
;
;
; CALLING SEQUENCE:
;   
;    hires_fntobj, hires, setup, obj_id, chip, [exp], /STD, /CHK [v1.0]'
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
;   chip    -  Blue (1), Green (2) OR Red (3) chip
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
;  FNTOBJ KEY =   MIN_PCA, _EXTRA
;  /TRIMRED - Deals with only a few orders on the RED chip + bad
;             trace.  Probably obsolete.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_fntobj, hires, 1L, 1L, /CHK, 
;
; PROCEDURES/FUNCTIONS CALLED:
;  hires_rectify
;  find_npeaks
;  hires_fweight
;  hires_fntobj_off
;  trace_crude
;  hires_basis
;
; REVISION HISTORY:
;   23-Sep-2003 Written by JXP (combined aspects of fndobj with trcobj)
;   10-Nov-2003 Changed fntobj_off, SMB
;   04-Feb-2013   Modified from HIRES by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_fntobj, hires, setup, obj_id, ichip, exp, $
                  CHK=chk, NOCLOB=noclob, OBJAPER=objaper, $
                  MINPCA=minpca, TRIMRED=trimred, $
                  DEBUG=debug, STD=std, _EXTRA=extra

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'hires_fntobj, hires, setup, obj_id, chip, [exp], /STD, /CHK, OBJAPER= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  ;; QA
  if keyword_set(CHK) then window, 0, title='hires_fntobj'
  if not keyword_set(FCOEFF) then fcoeff = 5L
  if not keyword_set(FULLTHR) then fullthr = 500.
  if setup LT 10 then fqa = 'QA/Obj0'+strtrim(setup,2) $
  else fqa = 'QA/Obj'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

  ;;  Find all relevant obj images
  for nn=0L,n_elements(ichip)-1 do begin
      chip = ichip[nn]
      
      trimtype = strtrim(hires.type,2)
      
      indx = where(hires.flg_anly NE 0 AND hires.setup EQ setup AND $
                   hires.chip EQ chip[0]  AND $
                   hires.obj_id EQ obj_id AND $
                   (trimtype EQ 'OBJ' OR trimtype EQ 'STD'), nindx)
      if nindx EQ 0 then begin
          print, 'hires_fntobj: No images to find obj for!', obj_id, ' ' + $
                 'and chip ',$
                 chip
          stop
          return
      endif
  
; APERTURES
      ;; Setting this automatically;  Not optimal...
      if not keyword_set( OBJAPER ) then begin
          if not keyword_set( STD ) then objaper = [0.60,0.60] $
          else objaper = [0.85,0.85] ; STD
      endif
      print, 'hires_fntobj: Using aperture ', objaper
      
      ;;  Read in order structure
      ordr_str = hires_getfil('ordr_str', setup, CHIP=chip)
      sv_ordrstr = ordr_str
      nordr = n_elements(ordr_str)
      
      ;;  Exposures
      if not keyword_set(exp) then exp = lindgen(nindx)

      ;; PCA
      case chip of 
          3: if not keyword_set(MINPCA) then MINPCA=7
          else: 
      endcase
      

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
      for q=0L,n_elements(exp)-1 do begin
          ;; Flag
          flg_objstr = 0
          ;; Look for obj file
          objfil = hires[indx[exp[q]]].obj_fil
          
          ;; xyoff
          ordr_str = sv_ordrstr
          ordr_str.lhedg = sv_ordrstr.lhedg + $
                           hires[indx[exp[q]]].arc_xyoff[0]
          ordr_str.rhedg = sv_ordrstr.rhedg + $
                           hires[indx[exp[q]]].arc_xyoff[0]
          
          ;; Cen ref
          cen_ref = (ordr_str.lhedg + ordr_str.rhedg)/2. 
          if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
              print, 'hires_fntobj: Using Obj structure -- ', objfil
              if x_chkfil(objfil+'*') EQ 0 then begin
                  print, 'hires_fntobj: File ', objfil, $
                         'doesnt exist! Turn off NOCLOB'
                  return
              endif
              objstr = xmrdfits(objfil, 1, STRUCTYP='hiresobjstrct', /silent)
              nobj = n_elements(objstr)
              flg_objstr = 1
          endif else begin
              ;; Set obj_fil
              ;; objfil = 'Extract/Obj_'+hires[indx[exp[q]]].img_root
              objfil = hires_getfil('obj_fil', /name, $
                                    FRAME=hires[indx[exp[q]]].frame,$
                                    CHIP=hires[indx[exp[q]]].chip)
              
              hires[indx[exp[q]]].obj_fil = objfil
              ;; Create objects
              tmp = { hiresobjstrct }
              objstr = replicate(tmp, nordr)
              objstr.UT = ' '
              objstr.field = ' '
              objstr.img_fil = ' '
              objstr.arc_fil = ' '
              objstr.exp = hires[indx[exp[q]]].exp
              objstr.field = hires[indx[exp[q]]].Obj
              nobj = 0L
          endelse
          objstr.flg_anly = 1
          
          ;; Read IMG+VAR
          imgfil = hires[indx[exp[q]]].img_final
          a = findfile(imgfil+'*', count=na)
          if na EQ 0 then begin
              print, 'hires_fntobj: No Final file ', imgfil
              print, '... try running hires_proc'
              continue
          endif
          objstr.img_fil = imgfil
          print, 'hires_fntobj: Opening files... ', imgfil
          img = xmrdfits(imgfil, /silent)
          ivar = xmrdfits(imgfil, 1, /silent)
          sz_img = size(img, /dimensions)
          
          ;; Read ARC
;      flg_arc = 0
;      objstr.arc_fil = strtrim(hires[indx[exp[q]]].arc_img,2)
;      if q EQ 0 then begin
;;          img_arc = xmrdfits(strtrim(hires[indx[exp[q]]].arc_img,2), /silent) 
;      endif else begin
;          if hires[indx[exp[q]]].arc_img NE hires[indx[exp[q-1]]].arc_img then $
;            img_arc = xmrdfits(hires[indx[exp[q]]].arc_img, /silent) $
;          else flg_arc=1
;      endelse
            
          qafil = hires_getfil('qa_fntobj', setup, CHIP=hires[indx[exp[q]]].chip,$
                               FRAME=hires[indx[exp[q]]].frame)
          ;; FIND+TRACE OBJ
          x_fntobj, objstr, ordr_str, img, ivar, qafil, $
                    CHK=chk, NOCLOB=noclob, OBJAPER=objaper, FCOEFF=fcoeff, $
                    MINPCA=minpca, FULLTHR=fullthr, DEBUG=debug, STD=std, _EXTRA=extra
          
          ;; TRIMRED
          if keyword_set(TRIMRED) then begin
              nobjs = n_elements(objstr)
              dff = objstr[nobjs-1].trace - shift(objstr[nobjs-1].trace,1)
              gd = where(objstr[nobjs-1].trace GT 0. and dff LT 0., ngd)
              npix = n_elements(objstr[nobjs-1].trace)
              ;; Spline
;              splin = spl_init(findgen(ngd)+gd[0], $
;                               objstr[nobjs-1].trace[gd], /double)
;              objstr[nobjs-1].trace[0:gd[0]-1] = $
;                spl_interp(findgen(ngd)+gd[0], objstr[nobjs-1].trace[gd], $
;                           splin, findgen(gd[0]), /doub)
              ;; Fit
              fitstr = x_setfitstrct(NORD=3, FUNC='POLY')
              fit = x_fit(findgen(ngd)+gd[0], $
                          objstr[nobjs-1].trace[gd], FITSTR=fitstr)
              objstr[nobjs-1].trace[0:gd[0]+49] = $
                x_calcfit(findgen(gd[0]+50), FITSTR=fitstr)
              
              print, 'hires_fntobj: Trimming the reddest order'
              print, 'hires_fntobj: Note the QA plot is NOT updated!!'
          endif
          
          ;; Write Obj structure
          print, 'hires_fntobj: Creating ', objfil
          mwrfits, objstr, objfil, /create, /silent
          spawn, 'gzip -f '+objfil
      endfor
  endfor
  
;  DONE
  print, 'hires_fntobj: All done! '
  return
end


