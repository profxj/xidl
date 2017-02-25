;+ 
; NAME:
; x_normflat   
;    Version 1.0
;
; PURPOSE:
;    Normalizes a multi-slit flat given the slit structure
;
; CALLING SEQUENCE:
;   
;   nrmflat = x_normflat( flat, slitstr, [var, nrmvar] )
;
; INPUTS:
;   flat       - Flat image or fits file
;   slistr     - Slit structure
;   [var]      - Variance array   
;
; RETURNS:
;   nrmflat    - Normalized flat   
;   nrmvar     - Normalized variance
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  YMED      -  Number of pixel offset from yedg to take median of
;  YNRM      -  Number of pixels outside yedg to normalize
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   nrmflat = x_normflat( flat, slitstr )
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Feb-2002 Written by JXP
;   15-Feb-2002 Revised to allow normalization of original
;   27-Mar-2004 Revised by MRB to not analyze flg_anly==0 slits
;-
;------------------------------------------------------------------------------

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function x_normflat, flat, slitstr, var, nrmvar, $
                     YMED=ymed, YNRM=ynrm, ORIG=orig, $
                     SILENT=silent, OUTNRM=outnrm


;  Error catching
  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'nrmflat = x_normflat(flat, slitstr, [var,, nrmvar], YMED=, YNRM=,'
    print, '       /ORIG, /SILENT, OUTNRM=) [v1.0]'
    return, -1
  endif 


;  Optional Keywords

  if not keyword_set( YMED ) then ymed = 2L
  if not keyword_set( YNRM ) then ynrm = 0L
  if keyword_set( VAR ) and arg_present(NRMVAR) then flg_var = 1 else flg_var=0

; Allow flat to be fits file

  dat = x_readimg(flat, /fscale)
  sz_img = size(dat, /dimensions)

  nslit = n_elements(slitstr)

; Sort to require descending order
  isrt = sort(slitstr.yedg_flt[1])
  isrt = reverse(isrt)

  if not keyword_set( ORIG ) then begin
      srt_yedg = slitstr[isrt].yedg_flt
      ; Detrermine median regions
      med_yedg = lonarr(2,nslit)
      for i=0L,nslit-1 do begin
          if(slitstr[i].flg_anly) then begin
; TOP
              if i GT 0 then $
                med_yedg[1,i] = (round(srt_yedg[1,i])-ymed) < $
                round(srt_yedg[0,i-1]) $
              else $
                med_yedg[1,i] = round(srt_yedg[1,i]) - ymed
; BOTTOM
              if i LT nslit-1 then $
                med_yedg[0,i] = (round(srt_yedg[0,i])+ymed) > $
                round(srt_yedg[1,i+1]) $
              else $
                med_yedg[0,i] = round(srt_yedg[0,i]) + ymed
          endif
      endfor

; Create normalization array
      tmp = fltarr(sz_img[0], nslit)
      for i=0L,nslit-1 do $
        if(slitstr[i].flg_anly) then $
        tmp[*,i] = djs_median(dat[*, med_yedg[0,i]:med_yedg[1,i]],2)
      ; Fit the norm array
      nrm = fltarr(sz_img[0], nslit)
      xx = findgen(sz_img[0])
      for i=0L,nslit-1 do $
        if(slitstr[i].flg_anly) then $
        nrm[*,i] = x_fitrej(xx, tmp[*,i], 'BSPLIN', 61, LSIGMA=5., HSIGMA=5.)

; NORMALIZE

  ; Find Edges of slit restricting against overlap
      nrm_yedg = temporary(med_yedg)
      for i=0L,nslit-1 do begin
        if(slitstr[i].flg_anly) then begin
; TOP
            if i GT 0 then $
              nrm_yedg[1,i] = (round(srt_yedg[1,i])+ynrm) < $
              (round(srt_yedg[0,i-1])-ynrm) $
            else $
              nrm_yedg[1,i] = round(srt_yedg[1,i]) + ynrm
; BOTTOM
            if i LT nslit-1 then $
              nrm_yedg[0,i] = (round(srt_yedg[0,i])-ynrm) > $
              (round(srt_yedg[1,i+1])+ynrm) $
            else $
              nrm_yedg[0,i] = (round(srt_yedg[0,i]) - ynrm) > 0
        endif
    endfor
    
; Create normalization image
    nrmimg = fltarr(sz_img[0], sz_img[1])
    for i=0L,nslit-1 do begin
        if(slitstr[i].flg_anly) then begin
            nrmimg[*,nrm_yedg[0,i]:nrm_yedg[1,i]] = $
              (nrm[*,i] # replicate(1.,nrm_yedg[1,i]-nrm_yedg[0,i]+1) )
        endif
    endfor
; Elminiate unused regions!
    a = where(nrmimg EQ 0)
    nrmflat = dat
    if a[0] NE -1 then begin
        nrmflat[a] = 0.
        nrmimg[a] = 1.
    endif
    
    if not keyword_set( SILENT) then print, 'x_normflat: Normalizing the flat'
    nrmflat = nrmflat / nrmimg
; VAR image
    if flg_var EQ 1 then nrmvar = var / (nrmimg)^2
    
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
endif else begin                ; ORIGINAL IMAGE
    srt_yedg = slitstr[isrt].yedg_orig
; Detrermine median regions
    med_yedg = lonarr(sz_img[0],2,nslit)
    for jj=0L,sz_img[0]-1 do begin
        for i=0L,nslit-1 do begin
            if(slitstr[i].flg_anly) then begin
                ;; TOP
                if i GT 0 then $
                  med_yedg[jj,1,i] = (round(srt_yedg[jj,1,i])-ymed) $
                  < round(srt_yedg[jj,0,i-1]) $
                else $
                  med_yedg[jj,1,i] = round(srt_yedg[jj,1,i]) - ymed
                ;; BOTTOM
                if i LT nslit-1 then $
                  med_yedg[jj,0,i] = (round(srt_yedg[jj,0,i])+ymed) $
                  > round(srt_yedg[jj,1,i+1]) $
                else $
                  med_yedg[jj,0,i] = round(srt_yedg[jj,0,i]) + ymed
            endif
        endfor
    endfor
    
    ;; Create normalization array
    tmp = fltarr(sz_img[0], nslit)
    for jj=0L, sz_img[0]-1 do begin
        for i=0L,nslit-1 do $
          if(slitstr[i].flg_anly) then $
          tmp[jj,i] = median(dat[jj, med_yedg[jj,0,i]:med_yedg[jj,1,i]])
    endfor
                                ; Fit the norm array
    nrm = fltarr(sz_img[0], nslit)
    xx = findgen(sz_img[0])
    for i=0L,nslit-1 do $
      if(slitstr[i].flg_anly) then $
      nrm[*,i] = x_fitrej(xx, tmp[*,i], 'BSPLIN', 61, LSIGMA=5., HSIGMA=5.)
    
    ;; NORMALIZE
    
    ;; Find Edges of slit restricting against overlap
    nrm_yedg = temporary(med_yedg)
    for jj=0L,sz_img[0]-1 do begin
        for i=0L,nslit-1 do begin
            if(slitstr[i].flg_anly) then begin
; TOP
                if i GT 0 then $
                  nrm_yedg[jj,1,i] = (round(srt_yedg[jj,1,i])+ynrm) < $
                  (round(srt_yedg[jj,0,i-1])-ynrm) $
                else $
                  nrm_yedg[jj,1,i] = round(srt_yedg[jj,1,i]) + ynrm
; BOTTOM
                if i LT nslit-1 then $
                  nrm_yedg[jj,0,i] = (round(srt_yedg[jj,0,i])-ynrm) > $
                  (round(srt_yedg[jj,1,i+1])+ynrm) $
                else $
                  nrm_yedg[jj,0,i] = (round(srt_yedg[jj,0,i]) - ynrm) > 0
            endif
        endfor
    endfor
    
; Create normalization image
    nrmimg = fltarr(sz_img[0], sz_img[1])
    for jj=0L,sz_img[0]-1 do begin
        for i=0L,nslit-1 do $
          if(slitstr[i].flg_anly) then $
          nrmimg[jj,nrm_yedg[jj,0,i]:nrm_yedg[jj,1,i]] = nrm[jj,i]
    endfor
    
; Elminiate unused regions!
    a = where(nrmimg EQ 0)
    nrmflat = dat
    if a[0] NE -1 then begin
        nrmflat[a] = 0.
        nrmimg[a] = 1.
    endif
    
    if not keyword_set( SILENT) then print, 'x_normflat: Normalizing the flat'
    nrmflat = nrmflat / nrmimg
; VAR image
    if flg_var EQ 1 then nrmvar = var / (nrmimg)^2
    
endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

if arg_present( OUTNRM ) then outnrm = nrm
; Release memory
delvarx, dat, nrmimg, nrm, tmp

return, nrmflat

end
