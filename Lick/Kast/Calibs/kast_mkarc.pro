;+ 
; NAME:
; kast_mkarc   
;     Version 1.0
;
; PURPOSE:
;    Process arc files (flatten, combine)
;
; CALLING SEQUENCE:
;  kast_mkarc, kast, setup, /CLOBBER
;
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
; /CLOBBER -- Clobber any previous processed image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_mkarc, kast
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-Aug-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_mkarc, kast, setup, CLOBBER=clobber

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'kast_mkarc, kast, setup, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
;  if not keyword_set( SIDE ) then side = 3L
  if not keyword_set( OVLBL ) then ovlbl = 'Raw/'

  c_set = strtrim(setup,2)

; Loop on side
  for qq=1,2 do begin
      ;; LBL
      if qq EQ 1 then lbl = 'B' else lbl = 'R'
      ;; Check for prior image
      outfil = 'Arcs/Arc'+lbl+c_set+'.fits'
      a = findfile(outfil+'*', count=na)
      if na NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'kast_mkarc: Arc ', outfil, ' exists.  Continuing'
          continue
      endif
      ;; Grab all Arc files
      gdarc = where(kast.mode EQ 1 AND kast.flg_anly NE 0 AND kast.side EQ qq AND $
                   kast.setup EQ setup AND strtrim(kast.type,2) EQ 'ARC', narc)
      if narc EQ 0 then begin
          print, 'kast_mkarc: No arc files! Continuing'
          continue
      endif

      ;; Median Combine
      if narc GT 1 then begin
          xcombine, ovlbl+kast[gdarc].img_root, fin_arc, head, $
            FCOMB=2, SCALE=kast[gdarc].exp, GAIN=kast[gdarc[0]].gain, $
            RN=kast[gdarc[0]].readno
      endif else fin_arc = mrdfits(ovlbl+kast[gdarc].img_root, /silent)

      ;; Flatten
      print, 'kast_mkarc: Flattening'
      ffil = 'Flats/Flat'+lbl+c_set+'.fits'
      if x_chkfil(ffil+'*') EQ 0 then begin
          print, 'kast_mkarc: Make Flat first!!'
          continue
      endif
      flat = xmrdfits(ffil, /silent)
      fin_arc = fin_arc / flat

      ;; Output
      mwrfits, fin_arc, outfil, head, /create, /silent
      spawn, 'gzip -f '+outfil
      print, 'kast_mkarc: Arc created ', outfil+'.gz'

  endfor
  print, 'kast_mkarc: All Done! '
  return
end
