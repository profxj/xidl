;+ 
; NAME:
; wfccd_cleanarc
;    Version 1.0
;
; PURPOSE:
;    Finds all bad wavelength values in the Arc image
;
; CALLING SEQUENCE:
;   
;   wfccd_cleanarc, wfccd, maskid, exp_id
;
; INPUTS:
;   wfccdstr    - WFCCD structure
;   mskid     - Long defining the mask to process
;
; RETURNS:
;
; OUTPUTS:
;   wfaimg      -  Cleaned WFCCD arc image (fits file)
;
; OPTIONAL KEYWORDS:
;  SLITSTR   - slit structure
;  WFASTR    - WFCCD arc structure
;  NOFITS   - Suppress fits output
;  OUTFIL   - Suppress fits output
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_cleanarc, wfccdstr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   02-May-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;  MAIN DRIVER
pro wfccd_cleanarc, wfccd, mask_id, exp_id, PER=per

;
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'wfccd_cleanarc, wfccd, mask_id, [exp_id], PER= [v1.0]'
    return
  endif 

; Set exp
  allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
              wfccd.mask_id EQ mask_id, nexp)
  if keyword_set(exp_id) then exp = allexp[exp_id] else exp=allexp[0]

;  Optional Keywords

  if not keyword_set( PER ) then per = 1.e-4

;  ARCFIL
  i = strpos(wfccd[exp].arc_fil, '_')
  arcfil = 'Arcs/ArcI'+strmid(wfccd[exp].arc_fil,i)
  arcfil = strtrim(arcfil,2)

; Read in Slit structure
  if not keyword_set( WFSLIT ) then $
    wfslit = xmrdfits(wfccd[exp].slit_fil,1,STRUCTYP='mslitstrct', /silent)
  
;;;;;;;;;;;;;;;;;;;;;;;
;;  Check to see if the Image exists

  a = findfile(arcfil, count=count)
  if count EQ 0 then begin
      print, 'wfccd_cleanarc: Arc image (', arcfil, ') does not exist!'
      return
  endif

;;;;
; INPUT

  wfaimg = xmrdfits(arcfil, /silent)
  sz = size(wfaimg, /dimensions)

;;;
; LOOP ON ROWS

  for q=0L,sz[1]-1 do begin
     case q of
         0L: begin  ; bottom row
             bad = where( abs(wfaimg[*,0]-wfaimg[*,1]) GT PER*wfaimg[*,1],nbad)
             if nbad NE 0 then wfaimg[bad,0L] = 0.
         end
         (sz[1]-1): begin ; top row
             bad = where( abs(wfaimg[*,sz[1]-1]-wfaimg[*,sz[1]-2]) $
                          GT PER*wfaimg[*,sz[1]-2],nbad)
             if nbad NE 0 then wfaimg[bad,sz[1]-1] = 0.
         end
         else: begin
             bad = where( abs(wfaimg[*,q]-wfaimg[*,q-1]) GT PER*wfaimg[*,q] AND $
                          abs(wfaimg[*,q]-wfaimg[*,q+1]) GT PER*wfaimg[*,q],nbad)
             if nbad NE 0 then wfaimg[bad,q] = 0.
         end
     endcase
 endfor

;;;
;  Zero out right edge of CCD where wavelength solution is fake

; CREATE 1D WAVE ARRAY
 nslit = n_elements(wfslit)
 for qq=0L,nslit-1 do begin
     wvoned = fltarr(sz[0])
     ycen = round(total(wfslit[qq].yedg_orig[*,*],2)/2.)
     for i=0L,sz[0]-1 do wvoned[i] = wfaimg[i,ycen[i]]
     dwv = wvoned - shift(wvoned,1)
; ??? MRB hack --- ignore very SMALL backwards steps in wavelngth
;    also, only check for bad stuff far enough to the right
     leftlimit=400L
     a = where(dwv[leftlimit:sz[0]-1] GE PER*wvoned,na)
     if na NE 0 then begin
         bdpx = a[0]+leftlimit 
         for ii=bdpx,sz[0]-1 do $
           wfaimg[ii,wfslit[qq].yedg_orig[ii,0]:wfslit[qq].yedg_orig[ii,1]]=0.
     endif 
 endfor
     
; OUTPUT

 mwrfits, wfaimg, arcfil, /create, /silent
 print, 'wfccd_cleanarc: All done! '

  return
end
