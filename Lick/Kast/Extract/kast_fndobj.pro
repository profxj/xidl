;+ 
; NAME:
; kast_fndobj   
;     Version 1.0
;
; PURPOSE:
;    Find object on the image and then trace using x_trace.
;  The code can look automatically or the user can look interactively.
;
; CALLING SEQUENCE:
; kast_fndobj, kast, setup, side, obj_id, [exp], $
;                /CHK, /NOCLOBBER, SCICLM=, /STD, /AUTO, SCIROW=
; 
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;   side  --  Specific camera [blue (1) vs. red (2)]
; obj_id  --  Object value
;  [exp]  --  Exposure indices
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  SCICLM     -- Column to look for object in
;  SCIROW     -- Approximate row expected for the object
;  /AUTO      -- Look for the object automatically
;  /NOCLOBBER -- Do not overwrite the object structure
;  /STD       -- Standard star
;  /CHK       -- Plot the image and the trace
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_fndobj, kast 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;   12-Nov-2005 Added SCICLM defaults for d55, blue G2 (600/4310) and
;               red 600/7500, KLC
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_fndobj, kast, setup, side, obj_id, exp, $
                 CHK=chk, NOCLOBBER=noclobber, SCICLM=sciclm, $
                 STD=std, AUTO=auto, SCIROW=scirow

;
  if  N_params() LT 4  then begin 
      print,'Syntax - ' + $
        'kast_fndobj, kast, setup, side, obj_id, [exp], SCICLM=  [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( SCIROW ) then scirow = 89L

;  Find all relevant obj
  if not keyword_set( STD ) then begin
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.obj_id EQ obj_id AND kast.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'kast_fndobj: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'kast_fndobj, kast, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
      indx = indx[obj_id]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

;  Loop

  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = strtrim(kast[indx[exp[q]]].obj_fil,2)
      if x_chkfil(objfil+'*') NE 0 and keyword_set(NOCLOBBER) then begin
          print, 'kast_fndobj: Using Obj structure -- ', objfil
          objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          print, 'kast_fndobj: Creating Obj structure -- ', objfil
          objfil = 'Extract/Obj_'+kast[indx[exp[q]]].img_root
          kast[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { specobjstrct }
          objstr = replicate(tmp, 10)
          objstr.slit_fil = ' '
          objstr.spec2d_fil = ' '
          objstr.img_fil = ' '
          objstr.UT = ' '
          objstr.instr_strct = ' '
          objstr.field = ' '
          objstr.exp = kast[indx[exp[q]]].exp
          nobj = 0L
          ;; SET SCICLM (XCEN)
          if not keyword_set( SCICLM ) then begin
              case strtrim(kast[indx[exp[q]]].splitter,2) of
                  'd46': begin
                      case side of
                          1: begin ;; blue
                              case strtrim(kast[indx[exp[q]]].grising,2) of
                                  '452/3306': SCICLM = 700L
                                  else: stop
                              endcase
                          end
                          2: begin ;; red
                              case strtrim(kast[indx[exp[q]]].grising,2) of
                                  '300/7500': SCICLM = 515L
                                  else: stop
                              endcase
                          end
                          else: stop
                      endcase
                  end
                  'mirror': begin
                      case side of
                          1: begin ;; blue
                              case strtrim(kast[indx[exp[q]]].grising,2) of
                                  '452/3306': SCICLM = 700L
                                  else: stop
                              endcase
                          end
                          2: stop
                          else: stop
                      endcase
                  end
                  'd55': begin
                      case side of
                          1: begin ;; blue
                              case strtrim(kast[indx[exp[q]]].grising,2) of
                                  '600/4310': SCICLM = 720L
                                  else: stop
                              endcase
                          end
                          2: begin ;; red
                              case strtrim(kast[indx[exp[q]]].grising,2) of
                                  '600/7500': SCICLM = 100L
                                  else: stop
                              endcase
                          end
                          else: stop
                      endcase 
                  end 
                  else: stop
              endcase
          endif
          objstr.xcen = round(sciclm)
      endelse

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+kast[indx[exp[q]]].img_root
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'kast_fndobj: No Final file ', imgfil
          stop
      endif
      objstr.spec2d_fil = imgfil
      img = xmrdfits(imgfil, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz_img = size(img, /dimensions)


      ;; Smash
      smsh = fltarr(sz_img[1])
      print, 'kast_fndobj: Smashing'
      smsh = djs_median(img[sciclm-20L:sciclm+20L,*],1)

      ;; Peaks
      print, 'kast_fndobj: Finding obj'
      gdsmsh = where(smsh GT 0.)
      peak = 1.
      x_fndobj, smsh[gdsmsh], center, NSIG=3., PEAK=peak, $
        PKWDTH=4L, EDGES=edges 
      gdpeak = gdsmsh[peak]
      npk = n_elements(peak)

      ;; Chk
;      if keyword_set( CHK ) then $
;        x_splot, smsh, /block, XTWO=gdpeak, $
;        YTWO=fltarr(n_elements(center)), psym_y2=1

      ;; Identify SCI
      mn = min(abs(gdpeak-scirow), sci)
      off = gdpeak[sci]-center[sci]

      ;; Setup Obj
      if flg_objstr EQ 0 then begin
          ;; Science obj
          objstr[0].xcen = sciclm
          objstr[0].ycen = center[sci]+off
          objstr[0].obj_id = 'a'
          objstr[0].aper[*] = 8.
          nobj = nobj + 1

          if not keyword_set( STD ) then begin
              ;; Other obj
              for jj=0L,npk-1 do begin
                  ;; Skip sci
                  if jj EQ sci then continue
                  ;; Set
                  off = float(gdpeak[jj]-peak[jj])
                  objstr[nobj].xcen = sciclm
                  objstr[nobj].ycen = center[jj]+off
                  objstr[nobj].obj_id = x_objnumid(nobj)
                  objstr[nobj].aper[*] = 8.
                  nobj = nobj + 1
              endfor
          endif
          ;; Reset
          objstr = objstr[0:nobj-1]
      endif
      
      ;; AUTO
      if not keyword_set( AUTO ) then begin
          dum = x_setapergui( smsh, OBJSTR=objstr )
          nobj = n_elements(objstr)
          ;; Set xcen for new obj
          newobj = where(objstr.xcen EQ 0., nnew)
          if nnew NE 0 then begin
              for i=0L,nnew-1 do begin
                  peak = round(objstr[newobj[i]].ycen)
                  objstr[newobj[i]].xcen = sciclm
                  objstr[newobj[i]].aper[*] = 8.
              endfor
          endif
      endif

      ;; TRACE
      for jj=0L,nobj-1 do begin

          trace = x_trace(img, objstr[jj].ycen, objstr[jj].xcen, $
                          VAR=var, TNORD=3, /CRUDE, /ROT, MXSHFT=0.02)
          ;; Test
          if keyword_set( CHK ) then begin
              rnd_trc2 = lindgen(sz_img[0])
              trc_msk = rnd_trc2 + round(trace)*sz_img[0]
              tmp = img
              tmp[trc_msk] = -10000
              xatv, tmp, /block
          endif
          
          ;; Save
          objstr[jj].trace[0:sz_img[0]-1] = trace

      endfor
          

      ;; Write Obj structure
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'kast_fndobj: All done! '
  return
end
