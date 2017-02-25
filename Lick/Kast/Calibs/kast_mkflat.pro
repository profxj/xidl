;+ 
; NAME:
; kast_mkflat   
;     Version 1.1
;
; PURPOSE:
;    Process flat file
;
; CALLING SEQUENCE:
;   
;  kast_mkflat, kast, slit
;
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
; /CLOBBER -- Clobber any previous processed image
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_mkflat, kast, 1L
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro kast_mkflat, kast, setup, CLOBBER=clobber
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'kast_mkflat, kast, setup, /CLOBBER [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( OVLBL ) then ovlbl = 'Raw/'

  c_set = strtrim(setup,2)

; BLUE
  for qq=1,2 do begin
      ;; LBL
      if qq EQ 1 then lbl = 'B' else lbl = 'R'
      ;; Check for prior image
      outfil = 'Flats/Flat'+lbl+c_set+'.fits'
      a = findfile(outfil+'*', count=na)
      if na NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'kast_mkflat: Flat ', outfil, ' exists.  Continuing'
          continue
      endif
      ;; Grab all QTZ files
      gdflat = where(kast.mode EQ 1 AND kast.flg_anly NE 0 AND kast.side EQ qq AND $
                   kast.setup EQ setup AND strtrim(kast.type,2) EQ 'QTZ', nflat)
      if nflat EQ 0 then begin
          print, 'kast_mkflat: No flat files! Continuing'
          continue
      endif

      ;; Median Combine
      if nflat GT 1 then begin
          xcombine, ovlbl+kast[gdflat].img_root, fin_flat, head, $
            FCOMB=2, SCALE='MED', GAIN=kast[gdflat[0]].gain, $
            RN=kast[gdflat[0]].readno
      endif else fin_flat = mrdfits(ovlbl+kast[gdflat].img_root, /silent)

      ;; Normalize
      print, 'kast_mkflat: Normalizing'
      sz = size(fin_flat, /dimensions)
      medflt = djs_median(fin_flat[*,(sz[1]/2)-20L:(sz[1]/2)+20L],2)
      if kast[gdflat[0]].splitter EQ 'open' OR $
        kast[gdflat[0]].splitter EQ 'mirror' then en = 100 else en = 15
      bset = bspline_iterfit(findgen(n_elements(medflt)),medflt, yfit=yfit,$
                             everyn=en) 
      fin_flat = fin_flat / ( yfit # replicate(1., sz[1]))

      ;; Output
      mwrfits, fin_flat, outfil, head, /create, /silent
      spawn, 'gzip -f ', outfil
      print, 'kast_mkflat: Flat created ', outfil+'.gz'

  endfor
  print, 'kast_mkflat: All Done! '
  return
end
