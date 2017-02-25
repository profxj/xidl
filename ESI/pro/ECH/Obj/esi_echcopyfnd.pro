;+ 
; NAME:
; esi_echcopyfnd   
;     Version 1.0
;
; PURPOSE:
;    Finds all objects and traces them
;
; CALLING SEQUENCE:
;   
;  esi_echcopyfnd, esi, obj_id
;
; INPUTS:
;   esi   -  ESI structure
;   indx  -  Indices of objects to process
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   FLAT  - Flat file
;   BIAS  - Bias frame
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_echcopyfnd, esi 
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   22-Jul-2002 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_echcopyfnd, esi, obj_id, exp, NOCLOB=noclob, OFFSET=offset, $
                   APER=aper, CHK=chk, CPYIDX=cpyidx

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'esi_echcopyfnd, esi, obj_id, [exp], /NOCLOB, OFFSET=, APER=, /CHK [v1.0]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set(OFFSET) then offset = 0.
  if not keyword_set(APER) then aper = 11.
  if not keyword_set(OBJ_ID) then obj_id = 0L
  if not keyword_set(CPYIDX) then begin
      print, 'esi_echcopyfnd: You must set CPYOBJ!'
      return
  endif

;  Find all relevant obj images
  indx = where(esi.flg_anly NE 0 AND esi.mode EQ 2 AND $
               esi.obj_id EQ obj_id AND strtrim(esi.type,2) EQ 'OBJ', nindx)
  if nindx EQ 0 then begin
      print, 'esi_echcopyfnd: No images to find obj for!', obj_id
      return
  endif
  
;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

; Opening template
  cpyobj = xmrdfits(esi[cpyidx].obj_fil, 1, STRUCTYP='esiobjstrct', /silent)

;  Loop

  for q=0L,n_elements(exp)-1 do begin
      ;; Flag
      flg_objstr = 0
      ;; Look for obj file
      objfil = esi[indx[exp[q]]].obj_fil
      ;; REFORDR
      refordr = esi[cpyidx].refordr
      esi[indx[exp[q]]].refordr = refordr
      if strlen(strtrim(objfil,2)) NE 0 and keyword_set(NOCLOB) then begin
          print, 'esi_echcopyfnd: Using Obj structure -- ', objfil
          objstr = xmrdfits(objfil, 1, STRUCTYP='esiobjstrct', /silent)
          nobj = n_elements(objstr)
          flg_objstr = 1
      endif else begin
          ;; Set obj_fil
          objfil = 'Extract/Obj_'+esi[indx[exp[q]]].img_root
          esi[indx[exp[q]]].obj_fil = objfil
          ;; Create objects
          tmp = { esiobjstrct }
          objstr = replicate(tmp, 50)
;          objstr.slit_fil = ' '
;          objstr.spec2d_fil = ' '
          objstr.img_fil = ' '
          objstr.UT = ' '
;          objstr.instr_strct = ' '
          objstr.field = ' '
          objstr.exp = esi[indx[exp[q]]].exp
          nobj = 0L
      endelse

      ;; Set IMG
      imgfil = 'Final/f_'+esi[indx[exp[q]]].img_root
      a = findfile(imgfil+'*', count=na)
      if na EQ 0 then begin
          print, 'esi_echfndobj: No Final file ', imgfil
          stop
          continue
      endif
      objstr.spec2d_fil = imgfil

      ;; Setup objstr (default aper to 11)
      if flg_objstr EQ 0 then begin
          print, 'esi_echcopyfnd: Tracing..'
          ;; ONLY SCIENCE

          for qq=0L,9 do begin
              objstr[qq].aper[*] = aper
              objstr[qq].obj_id = 'a'
              objstr[qq].order = qq

              ;; xcen
              objstr[qq].xcen = cpyobj[qq].xcen
              ;; ycen
              objstr[qq].ycen = cpyobj[qq].ycen + offset
              ;; TRACE
              objstr[qq].trace = cpyobj[qq].trace + offset
          endfor
          nobj = nobj + 1

          ;; CHK
          if keyword_set( CHK ) then begin
              tmp = xmrdfits(esi[indx[exp[q]]].img_final,0,/silent)
	      sz_img = size(tmp, /dimensions)
              for qq=0L,9 do begin
                  rnd_trc = round(objstr[qq].trace[0:sz_img[1]-1])
                  trc_msk = rnd_trc + lindgen(sz_img[1])*sz_img[0]
                  tmp[trc_msk] = -10000
              endfor
              print, objstr[refordr].xcen, objstr[refordr].ycen
              xatv, tmp, /block, min=-70, max=770
          endif

          ;; Take only the relevant structures
          objstr = objstr[0:nobj*10L-1]
      endif
      

      ;; Write Obj structure
      mwrfits, objstr, objfil, /create, /silent
      spawn, 'gzip -f '+objfil

  endfor
  
;  DONE
  print, 'esi_echcopyfnd: All done! '
  return
end
