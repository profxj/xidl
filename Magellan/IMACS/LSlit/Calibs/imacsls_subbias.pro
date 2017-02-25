;+ 
; NAME:
; imacsls_subbias   
;     Version 1.1
;
; PURPOSE:
;   Subtract the overscan from an image (or set of images)
;
; CALLING SEQUENCE:
;  imacsls_subbias, imacsls, indx
;
; INPUTS:
;   imacsls   -  IMACS structure
;   indx  -  Index numbers of frame to subtract (default output is OV)
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  BIASFIL= - Name of bias file (default: Bias/BiasS[I].fits)
;  OVROOT=  - Root name of OV file (default: OV/ov_ )
;  /FORCE   - Overwrite existing OV files 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_subbias, imacsls, [47L,48L,49L]
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   17-Jul-2002 Written by JXP
;   01-Feb-2003 Polished (JXP)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_subbias, imacsls, indx, OVROOT=OVROOT, BIASFIL=biasfil, FORCE=force
;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'imacsls_subbias, imacsls, indx, OVROOT=, BIASFIL=, /FORCE [v1.1]'
      return
  endif 
  
;  Optional Keywords
  if not keyword_set( OVROOT ) then ovroot = 'OV/ov_'
  if not keyword_set( BIASFIL ) then biasfil = 'Bias/Bias.fits'

; Open Bias file
;  print, 'imacsls_subbias: Using BIAS file: ', biasfil
;  if x_chkfil(biasfil+'*') EQ 0 then begin
;      print, 'imacsls_subbias: Bias file not found.  Create first!', biasfil
;      stop
;  endif
;  bias = xmrdfits(biasfil, /silent)
;  sz = size(bias, /dimensions)

; Loop

  for q=0,n_elements(indx)-1 do begin
      ;; Check for output
      outfil = ovroot+imacsls[indx[q]].img_root
      a = findfile(outfil+'*', count=na)
      if na NE 0 and not keyword_set( FORCE ) then begin
          print, 'imacsls_subbias: File ', outfil, ' found. Not resubtracting'
          imacsls[indx[q]].img_ov = outfil
          imacsls[indx[q]].flg_ov = 1
          continue
      endif
      ;; Open Raw image
      raw = xmrdfits(imacsls[indx[q]].rootpth+imacsls[indx[q]].img_root, 0, head, $
                    /silent, /fscale)
      print, 'imacsls_subbias: Subtracting the bias for ', $
        imacsls[indx[q]].rootpth+imacsls[indx[q]].img_root

      ;; Image edges, ov region
      x1 = 0
      x2 = (2048L/imacsls[indx[q]].cbin)-1
      a1 = x2+2
      y1 = 0
      y2 = (4096L/imacsls[indx[q]].rbin)-1
      a2 = y2+2
      ovimg = fltarr(x2-x1+1, y2-y1+1)
      ;; OV
      ov = djs_median(raw[a1:*,*],1)
      bset = bspline_iterfit(findgen(n_elements(ov)),ov, yfit=yfit,$
                             everyn=35, lower=2.5, upper=2.5) 
      raw[x1:x2,*] = $
        raw[x1:x2,*] - replicate(1.,x2-x1+1)#yfit
      ;; Bias row
      ov2 = djs_median(raw[x1:x2,a2:*],2)
      bset = bspline_iterfit(findgen(n_elements(ov2)),ov2, yfit=yfit,$
                             everyn=35, lower=2.5, upper=2.5) 
      ovimg = raw[x1:x2,y1:y2] - yfit#replicate(1.,y2-y1+1)

      ;; Update flg and header
      imacsls[indx[q]].flg_anly = 3
      sxaddpar, head, 'BIAS', 'T', biasfil
      
      ;; Write
      imacsls[indx[q]].img_ov = outfil
      imacsls[indx[q]].flg_ov = 1
      mwrfits, ovimg, outfil, head, /create, /silent
  endfor

  print, 'imacsls_subbias: All done!'

  return
end
