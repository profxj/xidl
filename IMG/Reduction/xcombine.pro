;+ 
; NAME:
; xcombine   
;   Version 1.2  
;
; PURPOSE:
;    Combines a set of images with a variety of options.
;
; CALLING SEQUENCE:
;   
;   xcombine, img, comb, [header], FCOMB=, GAIN=, RN=, MASKS=, SCALE=,
;     STATSEC=, MMEM=, SIGHI=, SIGLO=, MXITER=, OUTNM=, /SILENT, 
;     WEIGHT=, IMGINDX=
;
; INPUTS:
;   img -- Array of image names
;
; RETURNS:
;
; OUTPUTS:
;   comb -- Combined image
;
; OPTIONAL KEYWORDS:
;   fcomb -- Binary keyword describing combine method 
;          0,1 = Median, Avg
;          2 = SIG Rej
;          4 = Masks
;          8 = Weight
;   mmem -- Maximum memory to use at any one time [default: Approx 200Mb]
;   outnm -- File to write combined fits image
;   gain -- Gain (required for SIGMEDCLIP)
;   rn -- ReadNoise (required for SIGMEDCLIP)
;   scale -- scale:  Could be array of floats, 'MED', 'AVG', or
;            Keyword
;   statsec -- Section to perform stats on (for scale)
;   imgindx -- Index number of image [default: 0L]
;   siglo  -- Lower sigma for rejection [default: 3.]
;   sighi  -- Upper sigma for rejection [default: 3.]
;   
;
; OPTIONAL OUTPUTS:
;   header -- Header of one of the files
;
; COMMENTS:
;
; EXAMPLES:
;   xcombine, all_img, comb_img
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   18-Jul-2001 Written by JXP
;   28-Jul-2001 Added masks, sigma clipping
;   29-Jul-2001 Added scaling, median routines
;   07-Jan-2002 Fixed memory leak
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xcombine, img, cmbimg, header, FCOMB=fcomb, SCALE=scale, MMEM=mmem, $
              MASKS=masks, $
              SIGLO=siglo, SIGHI=sighi, MXITER=mxiter, OUTNM=outnm, $
              SILENT=silent, $
              GAIN=gain, RN=rn, STATSEC=statsec, IMGINDX=imgindx, WEIGHT=weight

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'xcombine, img, [cmbimg, header], FCOMB=, MMEM=, MASKS=, SIGLO=' 
      print, '                            SIGHI=, MXITER=, OUTNM=, GAIN= '
      print, '                            RN=, STATSEC=, WEIGHT= (v1.1)'
      return
  endif 
  
  nimg = n_elements(img)

;  Optional Keywords

  if not keyword_set( IMGINDX ) then imgindx = 0L
  if not keyword_set( RN ) then rn = 5.
  if not keyword_set( GAIN ) then gain = 1.
  if not keyword_set( MMEM ) then mmem = 2000. ;; 2Gb
  if not keyword_set( FCOMB ) then fcomb = 0
  if keyword_set( MASKS ) then begin
      if fcomb mod 8 LE 3 then print, 'fcomb not set for Masks; Proceeding without'
      if n_elements(masks) NE n_elements(img) then begin 
          print, 'Wrong number of mask images'
          if fcomb mod 8 GT 3 then fcomb = fcomb - 4
      endif
  endif else begin
      if fcomb mod 8 GT 3 then begin
          print, 'xcombine: fcomb set for Masks but none input! Proceeding without masks..'
          fcomb = fcomb - 4
      endif
  endelse

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Scale

  if keyword_set( SCALE ) then begin
      if size(scale,/type) EQ 4 OR size(scale,/type) EQ 5 then sclval = scale $
      else begin
          ; Fill in the scale array
          sclval = fltarr(nimg)
          case scale of 
              'MED' : begin
                  print, 'xcombine: Scaling the images with the median'
                  for i=0L,nimg-1 do begin
                      d1 = xmrdfits(img[i], IMGINDX, /silent, /fscale)
                      if keyword_set( STATSEC ) then begin
                          ssec = xregtovec(statsec,size(d1))
                          sclval[i] = float(median(d1[ssec[0]:ssec[1], $
                                                     ssec[2]:ssec[3]]))
                          delvarx, d1
                      endif else begin ;; Take ratio of entire image
                          if i EQ 0 then begin
                              d0 = temporary(d1)
                              sclval[i] = 1.
                          endif $
                          else sclval[i] = float(median(temporary(d1)/d0))
                      endelse
                  endfor
                  if keyword_set( d0 ) then delvarx, d0
              end
              'AVG' : begin
                  print, 'xcombine: Scaling the images with the average'
                  for i=0L,nimg-1 do begin
                      d1 = xmrdfits(img[i], IMGINDX, /silent, /fscale)
                      if keyword_set( STATSEC ) then begin
                          ssec = xregtovec(statsec,size(d1))
                          sclval[i] = float(mean(d1[ssec[0]:ssec[1], $
                                                      ssec[2]:ssec[3]],/double))
                      endif else sclval[i] = float(mean(d1,/double))
                      delvarx, d1
                  endfor
              end
              else : begin   ; Search for keyword for elapsed time
                  for i=0L,nimg-1 do begin
                      h1 = xheadfits(img[i])
                      elap = sxpar(h1, scale)
                      if count EQ 0 then begin
                          print, 'Header card', scale, ' not found!'
                          return
                      endif
                      sclval[i] = elap
                  endfor
              end
          endcase
      endelse
;    Normalize scale values to the mean + invert
      nmscl = mean(sclval)
      sclval = nmscl / sclval   
  endif

          
;  Investigate the images


  h0 = xheadfits(img[0])
  n1_0 = sxpar(h0, 'NAXIS1')
  n2_0 = sxpar(h0, 'NAXIS2')
;  bit_0 = sxpar(h0, 'BITPIX')
  
  for i=1L,nimg-1 do begin
      h1 = xheadfits(img[i])
      n1 = sxpar(h1, 'NAXIS1')
      n2 = sxpar(h1, 'NAXIS2')
;      bitp = sxpar(h1, 'BITPIX')

;      if( n1 NE n1_0 OR n2 NE n2_0 OR bitp NE bit_0 ) then begin
      if( n1 NE n1_0 OR n2 NE n2_0) then begin
          print, 'xcombine: Images arent all the same dimension or type'
          return
      endif
  endfor

;  Set Parameters according to the combine method
  
  if( fcomb mod 4 GT 1) then begin
      if not keyword_set(SIGHI) then sighi = 3.0
      if not keyword_set(SIGLO) then siglo = 3.0
      if not keyword_set(MXITER) then mxiter = 5
  endif

;  Output image
  comb = fltarr(n1_0, n2_0)

  ;  Restrict on memory
  dsiz = 4

;  case bit_0 of 
;      -32 : dsiz = 4
;      16 : dsiz = 2
;      else : begin
;          print, 'xcombine not set to handle bit_0 = ', bit_0
;          stop
;      end
;  endcase

  nmeg = (n1_0 * n2_0)*dsiz / 1024.^2   ; Number of Meg per image

  ; Allow for Masks
  if fcomb mod 8 GT 3 then nmeg = nmeg + (n1_0 * n2_0) / 1024.^2

  maximg =  (mmem-30) / nmeg - 2   ; Subtract 2 for output
  if maximg LE 2 then begin
      print, 'xcombine: Get some more memory and come back later!', maximg
      return
  endif

  ; Max number of rows at a time
  nrow = round(float(n2_0)*float(maximg)/float(nimg)) < n2_0 
;  nrow = 2000L

; Create Mask array as necessary
  if fcomb mod 8 GT 3 then mskarr = bytarr(n1_0,nrow,nimg)


; Status
  print, 'Combining images: '
  for i=0,nimg-1 do print, img[i]

                                ; Loop!
  srow = 0L                     ; Starting row
  lrow = nrow-1L                ; Ending row
  while (lrow LE n2_0-1) do begin

; Status
      print, 'xcombine: Reading rows: ', $
        strtrim(srow,2)+'-'+strtrim(lrow,2), $
	' of ', strtrim(n2_0-1,2)

      ; Create the arrays
      imgarr = fltarr(n1_0,lrow-srow+1,nimg)
      if fcomb mod 8 GT 3 then mskarr = bytarr(n1_0,lrow-srow+1,nimg)

      ; Read in the arrays
      for i=0,nimg-1 do imgarr[*,*,i] = $
        xmrdfits(img[i],IMGINDX,range=[srow,lrow], /silent, /fscale)

      ; Scale
      if keyword_set( SCALE ) then begin
          print, 'xcombine: Scaling'
          for i=0L,nimg-1 do imgarr[*,*,i] = imgarr[*,*,i] * sclval[i]
      endif

      ; Read in the Mask files
      if fcomb mod 8 GT 3 then begin
          print, 'xcombine: Reading in the masks'
          for i=0L,nimg-1 do mskarr[*,*,i] = $
            xmrdfits(masks[i],IMGINDX,range=[srow,lrow], /silent)
      endif
    
      print, 'xcombine: Combining...'
      case fcomb of 
          0 : $                 ; Plain Median
            comb[*,srow:lrow] = djs_median(imgarr, 3)
          1 : $                 ; Plain Average
            comb[*,srow:lrow] = avg(imgarr, 2)
          2 : $                 ; Median with Sigma Rej
            comb[*, srow:lrow] = x_medsigclip(imgarr, $
                                             SIGLO=siglo, $
                                             SIGHI=sighi, $
                                             GAIN=gain,   $
                                             RN=rn,   $
                                             MAXITER=mxiter)
          3 : $                 ; Average with Sigma Rej
            comb[*, srow:lrow] = x_avsigclip(imgarr, $
                                             SIGLO=siglo, $
                                             SIGHI=sighi, $
                                             MAXITER=mxiter)
          6 : $                 ; Median with Sigma Rej and Masks
            comb[*, srow:lrow] = x_medsigclip(imgarr, $
                                             SIGLO=siglo, $
                                             SIGHI=sighi, $
                                             GAIN=gain,   $
                                             RN=rn,   $
                                             MAXITER=mxiter, $
                                             INMASK=mskarr)
          7 : $                 ; Average with Sigma Rej and Masks
            comb[*, srow:lrow] = x_avsigclip(imgarr, $
                                             SIGLO=siglo, $
                                             SIGHI=sighi, $
                                             MAXITER=mxiter, $
                                             INMASK=mskarr)
          else : begin
              print, 'xcombine: Not ready to handle fcomb = ', fcomb
              return
          end
      endcase
            
      ; Release memory
      delvarx, imgarr
      if fcomb mod 8 GT 3 then delvarx, mskarr

      ; Adjust the Rows
      if(lrow EQ n2_0-1) then break
      srow = lrow+1
      lrow = lrow+nrow < (n2_0-1)

  endwhile
          

; Output

  if keyword_set( OUTNM ) then mwrfits, comb, outnm, h0, /create
  if arg_present( header ) then header = h0
  if arg_present( cmbimg ) then cmbimg = comb else delvarx, comb

  print, 'xcombine: All done!'

end
