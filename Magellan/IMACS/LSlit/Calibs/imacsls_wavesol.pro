;+ 
; NAME:
; imacsls_wavesol   
;   Version 1.1
;
; PURPOSE:
;  Calibrate an arc line spectrum.  This can be done automatically.
;
; CALLING SEQUENCE:
;  imacsls_wavesol, imacsls, setup, obj_id, side, exp, /STD, $
;                  /SILENT, /AINTER, /AUTO, /CALIB, CALIBFIL=, NFIN=
;   
; INPUTS:
;   imacsls  - IMACS structure
;   setup    - Setup ID
;   obj_id   - Object ID
;   side     - Side (CCD)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /AUTO
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   spec = imacsls_wavesol('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro imacsls_wavesol, imacsls, setup, obj_id, side, exp, STD=std, $
                  SILENT=silent, AINTER=ainter, AUTO=auto, CALIB=calib, $
                  CALIBFIL=calibfil, NFIN=nfin

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'imacsls_extract, imacsls, setup, obj_id, side, [exp], /AUTO, /SILENT'
	print, '   /CALIB, /AINTER, NFIN=, CALIBFIL=  [v1.1]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( NFIN ) then nfin = 5L
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0

; Find objects
  if not keyword_set( STD ) then begin
      indx = where(imacsls.flg_anly NE 0 AND imacsls.mode EQ 1 AND $
                   imacsls.side EQ side AND imacsls.setup EQ setup AND $
                   imacsls.obj_id EQ obj_id AND imacsls.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'imacsls_wavesol: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(imacsls.flg_anly NE 0 AND imacsls.mode EQ 1 AND $
                   imacsls.side EQ side AND imacsls.setup EQ setup AND $
                   imacsls.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'imacsls_extract, imacsls, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
      indx = indx[obj_id]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

  if side EQ 1 then begin
      case strtrim(imacsls[indx[0]].grising,2) of 
          '300l': begin
              templt = getenv('XIDL_DIR')+'/Magellan/IMACS/LSlit/Calibs/imacslsfit_300lB.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/imacs_lslit.lst'
          end
          else: stop
      endcase
  endif else begin
      case strtrim(imacsls[indx[0]].grising,2) of 
          '300l': begin
              templt = getenv('XIDL_DIR')+'/Magellan/IMACS/LSlit/Calibs/imacslsfit_300lR.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/imacs_lslit.lst'
          end
          else: stop
      endcase
  endelse

;  Loop
  for q=0L,n_elements(exp)-1 do begin

      ;; Open objfil
      objfil = imacsls[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'imacsls_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='dblsobjstrct', /silent)

      ;; Nobj
      nobj = n_elements(objstr)
      if keyword_set( STD) then nobj = 1

      ;; Read arc
      arc = imacsls_getfil('arc_fil', 1, SIDE=side)
      sz = size(arc, /dimensions)

      for kk=0L,nobj-1 do begin
          ;; Extract
          arcsp = x_extract(transpose(arc), $
                            [objstr.xcen-7.,objstr.xcen+7.], $
                            objstr.trace[0:sz[1]-1], CAPER=objstr.ycen)

          ;; CALIBRATE
          if keyword_set( AUTO ) then begin
              ;; Grab template
              restore, templt
              ;; Cross-correlate
              step = lindgen(400L) - 200L
              corr = c_correlate((aspec<2000.), (arcsp<2000.), step, /double)
              mx = max(corr, imx)
              imx = step[imx]
              print, 'imacsls_wavesol: Offseting ', strtrim(imx,2), ' pixels'
              ;; Line list
              x_arclist, linlist, lines
              ;; ID
              x_templarc, arcsp, lines, calib, SHFT=imx
              ;; Peak up
              x_arcpeakup, arcsp, lines, FFIT=ffit, WV=wave, INTER=ainter, $
                NFIN=nfin
          endif else begin
              x_identify, arcsp, calib, LINELIST=linlist
              wave = x_calcfit(n_elements(arcsp), FITSTR=calib)
              if keyword_set( CALIBFIL ) then begin
                  aspec = arcsp
                  save, aspec, calib, filename=calibfil
              endif
          endelse

          ;; VACUUM
          if not keyword_set( NOVAC ) then airtovac, wave
          ;; Save
          objstr.wave[0:sz[1]-1] = float(wave)
      endfor
      ;; Write
      print, 'imacsls_wavesol: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

  endfor


  print, 'imacsls_extract: All Done!'
  return
end
