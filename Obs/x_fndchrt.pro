;+ 
; NAME:
; x_fndchrt   
;    Version 1.0
;
; PURPOSE:
;    Given an array of strings (QSO names), return info and plot 
;
; CALLING SEQUENCE:
;   
; 
;
; INPUTS:
;   strings - Array of strings
;
; RETURNS:
;   uniq  - Array of unique members
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   imsize - Arcmin of image (default is 5')
;    
;
; OPTIONAL OUTPUTS:
;  COUNT - number of unique strings
;
; COMMENTS:
;
; EXAMPLES:
;   x_fndchrt, 'targets.list'
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   21-Nov-2003 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_fndchrt, targlist, OUTDIR=outdir, imsize=imsize, survey=survey

  if  N_params() LT 1 then begin 
      print,'Syntax - ' + $
        'x_fndchrt, targlist, OUTDIR=, IMSIZE= (v1.0)'
      return
  endif 

  if not keyword_set( OUTDIR ) then outdir='./'
  if not keyword_set( IMSIZE ) then imsize = 5.
  if not keyword_set( SURVEY ) then survey = '2r'

  ;; Readlist
  readcol, targlist, nam, ra, dec, FORMAT='A,A,A'

  ;; Loop

  nobj = n_elements(nam)

  for q=0L,nobj-1 do begin
      ;; Grab ra, dec
      x_radec, ra[q], dec[q], rad, decd
      ;; Grab dss image
      querydss, [rad,decd], img, hdr, imsize=imsize, survey=survey
      ;; Write to fits
      flg = 1
      if n_elements(img) GT 1 then begin
          if keyword_set(SVFITS) then imnm=nam[q]+'.fits' else imnm='tmp.fits'
          mwrfits, img, imnm, hdr, /create 
      endif else begin
          print, 'Image not found!  Try survey=''1'''
          flg=0
      endelse
      ;; Showfits
      if flg NE 1 OR keyword_set( NOPS ) then continue
      spwncmd = 'showfits -objnm='+strtrim(nam[q],2)+' -ps -fi=' $
        +outdir+strtrim(nam[q],2)+'.ps '+imnm 
      spawn, spwncmd
  endfor

  return
end


