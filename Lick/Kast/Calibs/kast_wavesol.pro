;+ 
; NAME:
; kast_wavesol   
;   Version 1.1
;
; PURPOSE:
;    
;
; CALLING SEQUENCE:
;   
;   spec = x_apall(ydat, [head])
;
; INPUTS:
;   ydat       - Values 
;   [head]     - Header
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   wave       - wavelength array
;   OVR        - String array for ov region:  '[2050:2100, *]'
;   ERROR      - Variance array
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   spec = kast_wavesol('spec.fits')
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_wavesol, kast, setup, side, obj_id, exp, STD=std, $
                  SILENT=silent, AINTER=ainter, AUTO=auto, CALIB=calib, $
                  CALIBFIL=calibfil, NFIN=nfin

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_extract, kast, setup, side, obj_id, [exp]  [v1.0]'
    return
  endif 

;  Optional Keywords
  if not keyword_set( NFIN ) then nfin = 5L
  if not keyword_set( SETUP ) then setup = 0
  if not keyword_set( SIDE ) then side = 1
  if not keyword_set( OBJ_ID ) then obj_id = 0

; Find objects
  if not keyword_set( STD ) then begin
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.obj_id EQ obj_id AND kast.type EQ 'OBJ', nindx)
      if nindx EQ 0 then begin
          print, 'kast_wavesol: No images to find obj for!', obj_id
          return
      endif
  endif else begin  ; STANDARD STAR
      indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                   kast.side EQ side AND kast.setup EQ setup AND $
                   kast.type EQ 'STD', nindx)
      if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
          print,'Syntax - ' + $
            'kast_extract, kast, setup, side, EXP, /STD   [v1.0]'
          return
      endif 
      indx = indx[obj_id]
      nindx = 1L
  endelse

;  Exposures
  if not keyword_set(exp) then exp = lindgen(nindx)

  if side EQ 1 then begin
      case strtrim(kast[indx[0]].grising,2) of 
          '452/3306': begin
              templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g1.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue3.lst'
          end
          else: stop
      endcase
  endif else begin
      case strtrim(kast[indx[0]].grising,2) of 
          '300/7500': begin
              templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_3007500.idl'
              linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_red.lst'
          end
          else: stop
      endcase
  endelse

;  Loop
  for q=0L,n_elements(exp)-1 do begin

      ;; Open objfil
      objfil = kast[indx[exp[q]]].obj_fil 
      if x_chkfil(objfil+'*') EQ 0 then begin
          print, 'kast_extract: No Obj file ', objfil
          continue
      endif
      objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)

      ;; Nobj
      nobj = n_elements(objstr)
      if keyword_set( STD) then nobj = 1

      ;; Read IMG+VAR
      imgfil = 'Final/f_'+kast[indx[exp[q]]].img_root
      if x_chkfil(imgfil+'*') EQ 0 then begin
          print, 'kast_extract: No Final file ', imgfil, '  Continuing...'
          continue
      endif
      img = xmrdfits(imgfil, /silent)
      var = xmrdfits(imgfil, 1, /silent)
      sz = size(img, /dimensions)

      ;; Read arc
      arc_fil = kast[indx[exp[q]]].arc_fil
      if x_chkfil(arc_fil+'*') EQ 0 then begin
          print, 'kast_extract: No Arc file ', arc_fil, '  Continuing...'
          continue
      endif
      arc = x_readimg(arc_fil, /fscale)

      for kk=0L,nobj-1 do begin
          ;; Extract
          arcsp = x_extract(arc, [objstr[kk].ycen-7.,objstr[kk].ycen+7.], $
                            objstr[kk].trace[0:sz[0]-1], CAPER=objstr[kk].xcen)

          ;; CALIBRATE
          if keyword_set( AUTO ) then begin
              ;; Grab template
              restore, templt
              ;; Cross-correlate
              step = lindgen(400L) - 200L
              corr = c_correlate((aspec<2000.), (arcsp<2000.), step, /double)
              mx = max(corr, imx)
              imx = step[imx]
              print, 'kast_wavesol: Offseting ', strtrim(imx,2), ' pixels'
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
          objstr[kk].wave[0:sz[0]-1] = float(wave)
      endfor
      ;; Write
      print, 'kast_wavesol: Updating ', objfil
      mwrfits, objstr, objfil, /create
      spawn, 'gzip -f '+objfil

  endfor


  print, 'kast_extract: All Done!'
  return
end
