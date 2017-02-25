; NAME:
; imacsls_fndobj   
;     Version 1.1
;
; PURPOSE:
;    Finds objects in the silt and traces them.  I recommend using
;  the interactive mode to identify (and mask) serendipitous objects.
;
; CALLING SEQUENCE:
;   
;  imacsls_fndobj, imacsls, setup, obj_id, [exp], REFWV=, SCICLM=,
;    /INTER, /STD, /CHK, MAXOFF=, /NOCLOB, NSIG=, /DEBUG
;
; INPUTS:
;   imacsls -  IMACS long slit structure
;   setup   -  Setupe ID
;   obj_id  -  Object ID  (e.g. 0L, 1L, etc)
;   [exp]   -  Exposure frames (e.g. [0L, 1L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /STD     - Find object for standard star
;  /CHK     - Show overall trace
;  REFWR=   - Reference row for identifiying the obj 
;  /INTER   - Interactively identify objects and set apertures
;             (recommended)
;  MAXOFF=  - Minimum offset between max obj and center of slit
;             (default: 5.)
;  NSIG=    - Number of sigma significance of the object
;  SCICLM=  - Expected column for the object
;  
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_fndobj, imacsls, 1L, 0L, /CHK, /INTER 
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP 
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_fndobj, imacsls, setup, obj_id, side, exp, REFRW=refrw, $
                    SCICLM=sciclm, INTER=inter, $
                    CHK=chk, MAXOFF=maxoff, NOCLOB=noclob, $
                    NSIG=nsig, DEBUG=debug, STD=std

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'imacsls_fndobj, imacsls, setup, obj_id, side, [exp], /INTER, /STD, /CHK [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(MAXOFF) then maxoff=5.
  if not keyword_set( NSIG ) then nsig = 3.
  if not keyword_set( STD ) then type = 'OBJ' else type = 'STD'

;  Find all relevant obj images
  indx = where(imacsls.flg_anly NE 0 AND imacsls.setup EQ setup AND $
               imacsls.side EQ side  AND $
               imacsls.obj_id EQ obj_id AND strtrim(imacsls.type,2) EQ type, $
               nindx)
  if nindx EQ 0 then begin
      print, 'imacsls_fndobj: No images to find obj for!', obj_id, ' and side ',$
        side
      return
  endif

  if not keyword_set( SCICLM ) then begin
      ;; Key column
      case side of
          1: begin              ; blue
              case imacsls[indx[0]].slit of
                  0.75: sciclm = round(1175. / imacsls[indx[0]].cbin)
                  else: stop
              endcase
          end
          2: begin ; red
              case imacsls[indx[0]].slit of
                  0.75: sciclm = round(873. / imacsls[indx[0]].cbin)
                  else: stop
              endcase
          end
          else: stop
      endcase
  endif

  ;; REFRW
  if not keyword_set(REFRW) then refrw = 3450L/imacsls[indx[0]].rbin
  
; APERTURES
  if not keyword_set( OBJAPER ) then objaper = [12.,12.]/imacsls[indx[0]].cbin

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  Loop
  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = imacsls[indx[exp[q]]].obj_fil

      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
          print, 'imacsls_fndobj: Using Obj structure -- ', objfil
          if x_chkfil(objfil+'*') EQ 0 then begin
              print, 'imacsls_fndobj: File ', objfil, $
                'doesnt exist! Turn off NOCLOB'
              return
          endif
          objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          objfil = 'Extract/Obj_'+imacsls[indx[exp[q]]].img_root
          imacsls[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          objstr = { dblsobjstrct }
          objstr.slit_fil = ' '
          objstr.spec2d_fil = ' '
          objstr.img_fil = ' '
          objstr.UT = ' '
          objstr.instr_strct = ' '
          objstr.field = ' '
          objstr.exp = imacsls[indx[exp[q]]].exp
          objstr.obj_id = 'a'
          nobj = 1L
      endelse

      ;; Read IMG+IVAR
      imgfil = 'Final/f_'+imacsls[indx[exp[q]]].img_root
      a = findfile(imgfil+'*', count=na)
      if na EQ 0 then begin
          print, 'imacsls_fndobj: No Final file ', imgfil
          stop
          continue
      endif
      objstr.spec2d_fil = imgfil
      objstr.aper = objaper
      print, 'imacsls_fndobj: Opening files... ', imgfil
      img = xmrdfits(imgfil, /silent)
      ivar = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)

      ;; FIND OBJ
 
      ;; Smash
      print, 'imacsls_fndobj: Smashing...'
      smsh = djs_median(img[sciclm-50:sciclm+50,refrw-7:refrw+7],2)

      ;; Interactive
      if keyword_set( INTER ) then begin
;          tmp = { specobjstrct }
;          tmp.obj_id
          if keyword_set( tmpobj ) then delvarx, tmpobj
          dum = x_setapergui( smsh, OBJSTR=tmpobj )
          npk = n_elements(tmpobj)
          sciobj = where(tmpobj.obj_id EQ 'a', nsci)
          sci = sciobj[0]
          if nsci EQ 0 then begin
              print, 'imacsls_fndobj: You failed to identify the sci obj!'+$
                '  Returning!'
              return
          endif
          center = tmpobj[sci].ycen
          off_sci = center-50+sciclm
      endif else begin ;; AUTO
          x_fndobj, smsh, center, PEAK=peak, NSIG=nsig, PKWDTH=3L, $
            EDGES=edges, NORDB=2L
          ;; Check for center
          if center[0] EQ -1. AND n_elements(center) EQ 1 then begin
              print, 'imacsls_fndobj: No obj found!!'
              print, "giving up on X's crappy fndobj, let's use find_npeaks"
              center = find_npeaks(smsh, nfind=1, width=3L, ypeak=peak)
              print, 'found a peak with height', peak[0], ' at pixel ', center[0]
          end
          ;; Set Science obj as max within 20 pixels of center
          nearobj = where(abs(center-50) LT maxoff, nnear)
          if nnear EQ 0 then begin
              print, 'imacsls_fndobj: No obj near center of order!!  Punting'
              print, 'imacsls_fndobj: Try running with interactive or choosing' + $
                ' another region of the chip for centroiding.' 

              stop
              return
          endif
          ;; Find Max near the center
          mx = max(smsh[center[nearobj]], sci)
          sci = nearobj[sci]
          off_sci = center[sci]-50+sciclm
      endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; TRACING

      ;; xcen, ycen
      objstr.xcen = off_sci
      objstr.ycen = refrw

      xcen_pos = trace_crude(img, ivar, yset=ycen_pos,XSTART=[off_sci],$
                             radius=2.5, ystart=[float(refrw)], xerr=xerr_pos, $
                             maxshifte=0.1)
      objstr.trace[0:sz_img[1]-1] = xcen_pos
                  
      ;; CHK
      if keyword_set( CHK ) then begin
          tmp = img
          rnd_trc = round(objstr.trace[0:sz_img[1]-1])
          trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
          tmp[trc_msk] = -10000
          xatv, tmp, /block, min=-70, max=1000
      endif
      
      ;; Write Obj structure
      print, 'imacsls_fndobj: Creating ', objfil
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'imacsls_fndobj: All done! '
  return
end
