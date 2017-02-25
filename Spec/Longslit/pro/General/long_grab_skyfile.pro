;+ 
; NAME:
; long_grab_skyfile   
;     Version 1.1
;
; PURPOSE:
;  Identify the Skyfile closest in wavelength
;
; CALLING SEQUENCE:
;   
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
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
;   23-May-2010 Written by JXP
;-
;------------------------------------------------------------------------------

function long_grab_skyfile, root, wv_mnx


  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'skyfile = long_grab_skyfile(root, wv_mnx) [v1.1]'
      return, -1
  endif 

  ;; Grab the files
  all_fil = findfile(root+'*.sav', count=nfil)
  if nfil EQ 0 then begin
     print, 'long_grab_skyfile: There are no sky files with this root!', iroot
     stop
  endif

  ;;  Parse for wavelengths
  wv_over = fltarr(nfil)
  for qq=0L,nfil-1 do begin
     prs = strsplit(all_fil[qq], '_', /extract)
     nprs = n_elements(prs)
     wv0 = float(prs[nprs-2])
     ipos = strpos(prs[nprs-1], '.sav')
     wv1 = float(strmid(prs[nprs-1],0,ipos))
     ;;
     wv_over[qq] = (wv1 < wv_mnx[1]) - max(wv0 > wv_mnx[0])
  endfor

  mx = max(wv_over,imx)

  return, all_fil[imx]
end

