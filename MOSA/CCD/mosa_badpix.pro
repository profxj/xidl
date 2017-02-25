;+ 
; NAME:
; mosa_badpix
;    Version 1.1
;
; PURPOSE:
;  Create a set of IRAF scripts to identify bad pixels in Image 3
;  of the MOSA images.
;
; CALLING SEQUENCE:
;   mosa_badpix, finlist, CLFIL=, /IM3, /IBAND, /IM4
;
; INPUTS:
;   finlst  - List of images
;
; RETURNS:
;
; OUTPUTS:
;  CLFIL=  -- Name of IRAF script file [default: 'badpix.cl']
;
; OPTIONAL KEYWORDS:
;  /IM3 -- Add extra lines for CCD 3
;  /IM4 -- Add extra lines for CCD 4
;  /IBAND -- I do not think this recommended
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   Apr-2003 Written by JXP
;------------------------------------------------------------------------------
pro mosa_badpix, finlst, CLFIL=clfil, IM3=im3, IBAND=iband, IM4=im4

  if not keyword_set( clfil ) then clfil = 'badpix.cl'

  ;; Read list
  readcol, finlst, fil, format='A'

  close, /all
  openw, 1, clfil
  ;; Loop
  for q=0L,n_elements(fil)-1 do begin
      ;; IM1
      if keyword_set( IBAND ) then $
      printf, 1, 'imreplace bpm'+strtrim(strmid(fil[q],3),2)+ $
        '/bpm_im1.pl[82:163,2020:2108] 1'
      ;; IM3
      if keyword_set( IM3) then $
        printf, 1, 'imreplace bpm'+strtrim(strmid(fil[q],3),2)+ $
        '/bpm_im3.pl[1353:1354,1:1986] 1'
      ;; IM4
      if keyword_set( IM4) then $
        printf, 1, 'imreplace bpm'+strtrim(strmid(fil[q],3),2)+ $
        '/bpm_im4.pl[*,2580:2586] 1'
  endfor

  close, /all
  return
end

  
