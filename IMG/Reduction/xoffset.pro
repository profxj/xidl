;+ 
; NAME:
; xoffset
;   Version 1.1
;
; PURPOSE:
;    Calculates integer offsets between a series of images.  The
;  program uses Sextractor to find objects in the image and then
;  uses a brute force algorithm to find the integer offset.
;
; CALLING SEQUENCE:
;   
;   xoffset, imglist, INSTR=
;
; INPUTS:
;   imglist    - List of images
;
; RETURNS:
;
; OUTPUTS:
;   offsets    -  ASCII file with offsets
;
; OPTIONAL KEYWORDS:
;   INSTR      - Instrument keyword
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xoffset, 'offR.list', INSTR='SITe3'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   26-Apr-2001 Written by JXP
;-
;------------------------------------------------------------------------------

pro xoffset, imglist, INSTR=instr, FWHM=fwhm, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'xoffset, imglist, INSTR=, FWHM=, OUTFIL= [v1.1]'
    return
  endif 

  if not keyword_set( INSTR ) then instr = ' '
;  Check for defaults
  chk = findfile('default.sex', COUNT=nfil)
  if nfil EQ 0 then begin
      print, 'xoffset: No default files found.  Copy them in!'
      return
  endif


;  Optional Keywords
  if not keyword_set( FWHM ) then fwhm = 1.0
  if not keyword_set( OUTFIL ) then outfil = 'offsets'

  close, 55
  openw, 55, outfil
  offx = 0L
  offy = 0L
  printf, 55, offx, offy, FORMAT='(2i6)'

; Instrument

  case INSTR of
      'SITe3': begin
         ax = -2000.
         ay = -1524.
         dx = 1.0
         dy = 1.0
         nx = 5000
         ny = 5000
         nstep = 1500
     end
      'MOSA': begin
         ax = -4500.
         ay = -4500.
         dx = 1.0
         dy = 1.0
         nx = 9000
         ny = 9000
         nstep = 4000
     end
     else : begin
         ax = -1024.
         ay = -1024.
         dx = 1.0
         dy = 1.0
         nx = 1000
         ny = 2048
         nstep = 1000
     end
 endcase
          

; Read images

  readcol, imglist, img_fil, FORMAT='A', /silent
  nimg = n_elements(img_fil)

; First image
  cmnd = 'sex '+img_fil[0]+' -SEEING_FWHM '+strtrim(fwhm,2)
  spawn, cmnd

; Parse the tmp file
  readcol, 'tmp.dat', xdum, ydum, area, SKIPLINE=5, /silent

  a = where(area GE 5, n1)
  x1 = xdum[a]
  y1 = ydum[a]

; Loop on remaining files

  for q=1L, nimg-1 do begin

      ; SEX first
      cmnd = 'sex '+img_fil[q]+' -SEEING_FWHM '+strtrim(fwhm,2)
      spawn, cmnd

      ; Parse the tmp file
      readcol, 'tmp.dat', xdum, ydum, area, SKIPLINE=5, /silent
      a = where(area GE 5, n2)
      x2 = xdum[a]
      y2 = ydum[a]
      ;
      zsv = fltarr(nx+1,ny+1)

      ; Match up
      for i=0L,n1-1 do begin
          zx = x1[i] - x2
          zy = y1[i] - y2

          k = round( (zx-ax)/dx)
          l = round( (zy-ay)/dy)

          gdpt = where(k GE 1 AND k LE nx AND $
                       l GE 1 AND l LE ny, ngd)
          ; Add em up
          for j=0L,ngd-1 do zsv[k[gdpt[j]],l[gdpt[j]]] = $
            zsv[k[gdpt[j]],l[gdpt[j]]] + 1
      endfor

      ; MAX
      zmx = max(zsv, imx)
      kmx = imx MOD (nx+1)
      lmx = imx / (nx+1)

      ; Offsets
      offx = long(ax + kmx*dx)
      offy = long(ay + lmx*dy)
      ; Output
      printf, 55, offx, offy, FORMAT='(2i6)'
  endfor
  close, 55

  ; Finished
  print, 'xoffset:  All done! Check ', outfil
end
      
      
      
