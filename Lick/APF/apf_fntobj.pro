;+ 
; NAME:
; apf_fntobj   
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
;    apf_fntobj, apf, setup, obj_id, [exp], /STD, /CHK [v1.0]'
;
; INPUTS:
;   apf   -  HIRES structure
;   setup   -  Setup ID
;   obj_id  -  Object ID  (e.g. 1L)  (or STD index if /STD set)
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
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   apf_fntobj, apf, 1L, 1L, /CHK, 
;
; PROCEDURES/FUNCTIONS CALLED:
;  apf_rectify
;  find_npeaks
;  apf_fweight
;  apf_fntobj_off
;  trace_crude
;  apf_basis
;
; REVISION HISTORY:
;   23-Sep-2003 Written by JXP (combined aspects of fndobj with trcobj)
;   10-Nov-2003 Changed fntobj_off, SMB
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro apf_fntobj, apf, setup, obj_id, exp, $
                  CHK=chk, NOCLOB=noclob, OBJAPER=objaper, $
                  MINPCA=minpca, TRIMRED=trimred, $
                  DEBUG=debug, STD=std, _EXTRA=extra

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'apf_fntobj, apf, setup, obj_id, [exp], /STD, /CHK, OBJAPER= [v1.1]'
      return
  endif 
  
;  Optional Keywords

  ;; QA
  if keyword_set(CHK) then window, 0, title='apf_fntobj'
  if not keyword_set(FCOEFF) then fcoeff = 5L
  if not keyword_set(FULLTHR) then fullthr = 500.
  if setup LT 10 then fqa = 'QA/Obj0'+strtrim(setup,2) $
  else fqa = 'QA/Obj'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

  ;;  Find all relevant obj images
  trimtype = strtrim(apf.type,2)
  indx = where(apf.flg_anly NE 0 AND apf.setup EQ setup AND $
               apf.obj_id EQ obj_id AND $
               (trimtype EQ 'OBJ' OR trimtype EQ 'STD'), nindx)
  if nindx EQ 0 then begin
     print, 'apf_fntobj: No images to find obj for!', obj_id
     stop
     return
  endif
  
; APERTURES
      ;; Setting this automatically;  Not optimal...
  if not keyword_set( OBJAPER ) then begin
     if not keyword_set( STD ) then objaper = [0.60,0.60] $
     else objaper = [0.85,0.85] ; STD
  endif
  print, 'apf_fntobj: Using aperture ', objaper
  
  ;;  Read in order structure
  ordr_str = apf_getfil('ordr_str', setup)
  sv_ordrstr = ordr_str
  nordr = n_elements(ordr_str)
  
  ;;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)
  

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
  for q=0L,n_elements(exp)-1 do begin
     ;; Flag
     flg_objstr = 0
     ;; Look for obj file
     objfil = apf[indx[exp[q]]].obj_fil
     
     ;; xyoff
     ordr_str = sv_ordrstr
     ordr_str.lhedg = sv_ordrstr.lhedg + $
                           apf[indx[exp[q]]].arc_xyoff[0]
     ordr_str.rhedg = sv_ordrstr.rhedg + $
                      apf[indx[exp[q]]].arc_xyoff[0]
     
     ;; Cen ref
     cen_ref = (ordr_str.lhedg + ordr_str.rhedg)/2. 
     if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
        print, 'apf_fntobj: Using Obj structure -- ', objfil
        if x_chkfil(objfil+'*') EQ 0 then begin
           print, 'apf_fntobj: File ', objfil, $
                  'doesnt exist! Turn off NOCLOB'
           return
        endif
        objstr = xmrdfits(objfil, 1, STRUCTYP='apfobjstrct', /silent)
        nobj = n_elements(objstr)
        flg_objstr = 1
     endif else begin
        ;; Set obj_fil
        ;; objfil = 'Extract/Obj_'+apf[indx[exp[q]]].img_root
        objfil = apf_getfil('obj_fil', /name, FRAME=apf[indx[exp[q]]].frame)
        
        apf[indx[exp[q]]].obj_fil = objfil
        ;; Create objects
        tmp = { apfobjstrct }
        objstr = replicate(tmp, nordr)
        objstr.UT = ' '
        objstr.field = ' '
        objstr.img_fil = ' '
        objstr.arc_fil = ' '
        objstr.exp = apf[indx[exp[q]]].exp
        objstr.field = apf[indx[exp[q]]].Obj
        nobj = 0L
     endelse
     objstr.flg_anly = 1
     
     ;; Read IMG+VAR
     imgfil = apf[indx[exp[q]]].img_final
     a = findfile(imgfil+'*', count=na)
     if na EQ 0 then begin
        print, 'apf_fntobj: No Final file ', imgfil
        print, '... try running apf_proc'
        continue
     endif
     objstr.img_fil = imgfil
     print, 'apf_fntobj: Opening files... ', imgfil
     img = xmrdfits(imgfil, /silent)
     ivar = xmrdfits(imgfil, 1, /silent)
     sz_img = size(img, /dimensions)
     
     qafil = apf_getfil('qa_fntobj', setup, FRAME=apf[indx[exp[q]]].frame)
     ;; FIND+TRACE OBJ
     x_fntobj, objstr, ordr_str, img, ivar, qafil, $
               CHK=chk, NOCLOB=noclob, OBJAPER=objaper, FCOEFF=fcoeff, $
               MINPCA=minpca, FULLTHR=fullthr, DEBUG=debug, STD=std, _EXTRA=extra
     
     ;; Write Obj structure
     print, 'apf_fntobj: Creating ', objfil
     mwrfits, objstr, objfil, /create, /silent
     spawn, 'gzip -f '+objfil
  endfor
  
;  DONE
  print, 'apf_fntobj: All done! '
  return
end


