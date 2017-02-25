;+
; NAME:
;   automated_fitwcs
;
; PURPOSE:
; Procedure that fits wcs solution to the image. 
; Not optimal, and fails sometimes. Need to improve.
;    
; CALLING SEQUENCE:
;
; INPUTS:
;   
; 
; 
; OPTIONAL INPUTS:
;  
; 
; OUTPUTS: 
;
; Generates wcs solution
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;    05-Jun-2012  Written by MF
;                 Revised by KK 
;                 Revised by MF (Can be significantly improved)
;
;-
;------------------------------------------------------------------------------
;; 


pro automated_fitwcs, filein

  ;;create temp file
  spawn, 'cp '+filein+' astro_tmp.fits' 

  ;;read first guess n ra-dec
  head=headfits('astro_tmp.fits',exten=0)
  ra=fxpar(head,"RA")
  dec=fxpar(head,"DEC")


  ;;Extra code to deal with the systematic offset so we can hopefully get better fits
  ;;This is quite empirical and not necesserely the best.

  x_radec, ra, dec, ra, dec 
  offset = xyoffset(head) 
  ;;Convert offset in pixels to arcseconds of ra and dec (roughly)
  offset_arcsec = offset*0.5/(60*60) 
  ra = ra-offset_arcsec[0]         ;try to set ra to center of image
  dec = dec-offset_arcsec[1]       ;try to set dec to center of image
  x_radec, ra, dec, ra, dec, /flip ;Converts ra and dec back to sexigasimal



  ;;fit new solution using astronomy.net
  splog, 'Fit new solution for ', filein
  spawn, ['solve-field', '-g','--no-plots','--overwrite','--ra',ra,'--dec',dec,'--radius','20','-d','50','astro_tmp.fits'], /noshell

  ;;accept new solution
  splog, 'Accept new solution'
  spawn, 'mv astro_tmp.new astro_tmp.fits', EXIT_STATUS=exit
  ;;check if solved or not
  if(exit EQ 1) then begin
     splog, 'WCS solutions not found!'
     spawn, 'touch '+filein+'NOWCSFOUND'
  endif

  newHeader = HEADFITS('astro_tmp.fits',EXTEN=0)
  MODFITS,fileIn,0,newHeader,EXTEN_NO=0
  
 ; spawn, 'cphead -w astro_tmp.fits '+filein
  
  
 ; fits=mrdfits(filein,0,head,/silent)
 ; sxaddpar, head, "HISTORY", 'ASTROMETRY FIT DONE WITH astronomy.net', AFTER='COMMENT'
  
  ;;clean
  spawn, 'rm -fr astro_tmp*'

  ;;clean /tmp (quite aggressive)
  spawn, 'rm -f /tmp/tmp.sanitized.*'
  spawn, 'rm -f /tmp/tmp.ppm.*'


end
