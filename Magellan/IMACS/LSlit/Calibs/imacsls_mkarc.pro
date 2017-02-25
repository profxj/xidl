;+ 
; NAME:
; imacsls_mkarc   
;     Version 1.0
;
; PURPOSE:
;    Process arc file and combine to produce one output file
;
; CALLING SEQUENCE:
; imacsls_mkarc, imacsls, setup, /CLOBBER
;
; INPUTS:
;   imacsls  -  IMACS structure
;   setup    -  Setup ID value
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /CLOBBER - Overwrite exisiting file if it exists
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_mkarc, imacsls, setup
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   09-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro imacsls_mkarc, imacsls, setup, CLOBBER=clobber

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'imacsls_mkarc, imacsls, setup, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OVLBL ) then ovlbl = 'OV/ov_'

  c_set = strtrim(setup,2)

; BLUE
  for qq=1,2 do begin

      ;; Check for prior image
      outfil = imacsls_getfil('arc_fil', 1, SIDE=qq, /name)
      a = findfile(outfil+'*', count=na)
      if na NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'imacsls_mkarc: Arc ', outfil, ' exists.  Continuing'
          continue
      endif

      ;; Grab all Arc files
      gdarc = where(imacsls.mode EQ 1 AND imacsls.flg_anly NE 0 AND $
                    imacsls.side EQ qq AND $
                    imacsls.setup EQ setup AND $
                    strtrim(imacsls.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'imacsls_mkarc: No arc files! Continuing'
          continue
      endif

      ;; OV subtract
      imacsls_subbias, imacsls, gdarc

      ;; Median Combine
      if narc GT 1 then stop $
      else fin_arc = xmrdfits(ovlbl+imacsls[gdarc].img_root, /silent)
;          xcombine, ovlbl+imacsls[gdarc].img_root, fin_arc, head, $
;;            FCOMB=2, SCALE=imacsls[gdarc].exp, GAIN=imacsls[gdarc[0]].gain, $
;            RN=imacsls[gdarc[0]].readno

      ;; Flatten
      print, 'imacsls_mkarc: Flattening'
      flat = imacsls_getfil('flat_fil', setup, side=qq)
      fin_arc = fin_arc / flat

      ;; Output
      mwrfits, fin_arc, outfil, head, /create, /silent
      spawn, 'gzip -f '+outfil
      print, 'imacsls_mkarc: Arc created ', outfil+'.gz'
      ;; OV
      if not keyword_set( SVOV ) then imacsls_delov, imacsls, gdarc

  endfor
  print, 'imacsls_mkarc: All Done! '
  return
end
