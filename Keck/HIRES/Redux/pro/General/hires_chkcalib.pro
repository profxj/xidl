;+ 
; NAME:
; hires_chkcalib
;     Version 1.1
;
; PURPOSE:
;    Verifies whether every science exposure has the requisite
;    calibration files.  If not, it squawks
;
; CALLING SEQUENCE:
;  hires_chkcalib
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  IMROOT=  -- Path to images
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_chkcalib
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   27-Oct-2005 Written by JXP
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_chkcalib, IMROOT=imroot

;
;  if  N_params() LT 1  then begin 
;      print,'Syntax - ' + $
;        'hires_chkcalib, IMROOT= [v1.0]'
;      return
;  endif 
  
;  Optional Keywords
  if not keyword_set( ETOLER ) then etoler = 0.00101
  if not keyword_set( XTOLER ) then xtoler = 0.00101
  if not keyword_set( IMROOT ) then imroot = 'hires'
  if not keyword_set( BADOUT ) then badout = 'chk_calib.txt'

  print, 'hires_chkcalib: Using tolerances ', etoler, xtoler

  ;; Create the structure
  hires_file = '[0-9][0-9][0-9][0-9].fits'
  img_archive = findfile('HI*'+hires_file+'*', count=narc) 
  img_tele = findfile(imroot+hires_file+'*',count=ntele) 
  img = [img_archive, img_tele]

  notempty = where(img NE '',nimg)
  
  if nimg EQ 0 then begin
      print, 'hires_chkcalib: No images !  Consider setting IMROOT'
      return
  endif
  img=img[notempty]

  ;; Create the structure
  hires_strct, hires, FILE_LIST=img, /NOFILE, /SILENT, /NOMKDIR

  ;; Run hires_setup
  hires_setup, hires
  uset = hires[uniq(hires.setup, sort(hires.setup))].setup
  gdset = where(uset LT 50)
  uset = uset[gdset]
  nset = n_elements(uset)

  close, /all
  openw, 19, badout
  printf, 19, 'Setup Deck  BLKF  ECHANG   XDANG  RBIN CBIN'
  for ii=0L,nset-1 do begin
      print, '------------------------------------------'
      set = uset[ii]
      xall = where(hires.setup EQ set)
      print, 'Setup = ', set
      printf, 19, '------------------------------------------'
      printf, 19, set, hires[xall[0]].decker, $
        hires[xall[0]].block, $
        hires[xall[0]].echangl, $
        hires[xall[0]].xdangl,$
        hires[xall[0]].rowbin, $
        hires[xall[0]].colbin, $
        FORMAT='(i2,5x,a2,1x,a6,1x,f7.4,1x,f7.4,2x,i2,3x,i2)'
      print, set, hires[xall[0]].decker, $
        hires[xall[0]].block, $
        hires[xall[0]].echangl, $
        hires[xall[0]].xdangl,$
        hires[xall[0]].rowbin, $
        hires[xall[0]].colbin, $
        FORMAT='(i2,5x,a2,1x,a6,1x,f7.4,1x,f7.4,2x,i2,3x,i2)'
      ;; Arcs?
      arcs = where(hires.setup EQ set and hires.type EQ 'ARC',narc)
      if narc NE 0 then begin
          printf, 19, 'Arcs: ', hires[arcs].img_root
      endif else begin
          print, 'hires_chkcalib: WARNING!!  No arcs for setup ', set
          printf, 19, 'hires_chkcalib: WARNING!!  No arcs for setup ', set
      endelse
      ;; Flats?
      flats = where(hires.setup EQ set and hires.type EQ 'TFLT',narc)
      if narc NE 0 then begin
          printf, 19, 'Flats: ', hires[flats].img_root
      endif else begin
          print, 'hires_chkcalib: WARNING!!  No flats for setup ', set
          printf, 19, 'hires_chkcalib: WARNING!!  No flats for setup ', set
      endelse
  endfor

  close, /all

  return
end
