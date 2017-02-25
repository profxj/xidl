;+ 
; NAME:
; hires_pixflat
;     Version 1.1
;
; PURPOSE:
;  This routine creates a pixel flat using the 'special' flat images
;  made by closing the cross-disperser cover.   It may be best for the
;  standard user to use one of the archvied pixel flats.
;
; CALLING SEQUENCE:
;  hires_pixflat, hires, setup, [side]
;
; INPUTS:
;   hires     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite Output MilkyFlat
;   /OVCLOB  - Overwrite OV files if they exist for the flats
;   /SVOV    - Save the OV files created during this step
;   /USEBIAS - Use bias frame in OV subtraction
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_pixflat, hires, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;   hires_getfil
;   hires_subbias
;   hires_pixflat_work
;
; REVISION HISTORY:
;   16-May-2003 Adapted by JXP from existing programs by SB
;   24-Feb-2004 Switched to a series of median/linear interpolations (SB)
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function hires_pixflat_work, flat_fils, flatstddev=flatstddev, $
                             GAIN=gain, RN=rn 
   
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'flat = hires_pixflat_work(files, GAIN=, RN=) [v1.1] '
      return, -1
  endif 

     if NOT keyword_set(rn) then rn = 3.1
     if NOT keyword_set(gain) then gain = 2.
     if NOT keyword_set(fitsout) then fitsout=0

     nfiles = n_elements(flat_fils)

     if nfiles GT 1 then begin
         print, 'Doing   1st of ', nfiles, format='(a,i3,$)'
         first = hires_pixflat_work(flat_fils[0], GAIN=gain, rn=rn) 
         nx = (size(first))[1]
         ny = (size(first))[2]
         final = fltarr(nx,ny, nfiles)
         final[*,*,0] = temporary(first)
         for i=1,nfiles -1 do begin
           print, 'Doing ', i+1, (i GT 2) ? 'th of ': (i EQ 1) ? $
                   'nd of ' : 'rd of ',  nfiles, format='(a,i3,a,i3)'
           final[*,*,i] = hires_pixflat_work(flat_fils[i], GAIN=gain, rn=rn)
         endfor
         print, 'hires_pixflat: Final Combine'
         x_statarray, final, 3, mean=flat, stddev=flatstddev, sigrej=2., $
                      /OVRIDE

         ; factor of 2 kludge due to gross underestimation
         flatstddev = 2.0 * flatstddev / sqrt(nfiles)
         return, flat
     endif


     ; Open file
     flat = xmrdfits(flat_fils, /silent)
     flat = flat*gain
     sz = size(flat,/dimensions)

     if not keyword_set(width) then begin
         rbin = round(4096. / sz[1])
         width = round(50L * (2./rbin)) + 1
     endif

     
     ;; Bufarr
;     ny = sz[1] / width
;     buff = fltarr(sz[0], (ny+5)*width)
;     bsz = size(buff, /dimensions)
;     buff[*, 2*width:2*width+sz[1]-1] = flat
;     
;     buff = (transpose(buff))[*]
     
     ;; Median
;     mbuff = median(buff, width)
     
     ;;
;     norm = reform(mbuff, bsz[1], bsz[0])
;     norm = transpose(norm)

     ;; 
     norm = transpose( x_medianrow(transpose(flat),width) )
     
     ;; Final
;     nflat = flat / norm[*, 2*width:2*width+sz[1]-1]
     nflat = flat / norm
     
     return, nflat

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_pixflat, hires, setup, chip, SVOV=svov, $
                  OVCLOB=ovclob, CLOBBER=clobber, USEBIAS=usebias

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_pixflat, hires, setup, [chip],  /SVOV, ' + $
        ' /USEBIAS, /CLOBBER, /OVCLOB [v1.1]'
      return
  endif 
  
  ;; QA
  if setup LT 10 then fqa = 'QA/Flats0'+strtrim(setup,2) $
  else fqa = 'QA/Flats'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  if not keyword_set( CHIP ) then chip = [1L,2L,3L]
  
  for ii=0L,n_elements(chip)-1 do begin
      qq = chip[ii]
      case qq of
          1: begin
              print, 'hires_pixflat: Creating BLUE Pixel flat' 
              xd = -3.1
          end
          2: begin
              print, 'hires_pixflat: Creating GREEN Pixel flat'
              xd = -4.1
          end
          3: begin
              print, 'hires_pixflat: Creating RED Pixel flat' 
              xd = -5.5
          end
      endcase

      ;; Outfil
      outfil = hires_getfil('pixflt_fil', setup, CHIP=qq, /name, CHKFIL=chkf) 
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'hires_pixflat: Milky flat exists, moving on..'
          continue
      endif

      ;; Grab flats
      gdflt = where(hires.chip EQ qq AND hires.flg_anly NE 0 AND $
                    strtrim(hires.type,2) EQ 'PFLT' AND $
                    hires.setup EQ setup AND $
                    abs(hires.xdangl- xd) LT 0.2, nflt)
      if nflt EQ 0 then begin
          print, 'hires_pixflat: No Flats of type PFLT found with this chip!' 
          return
      endif

      ;; Bias Subtract
      hires_subbias, hires, gdflt, CLOBBER=ovclob, USEBIAS=usebias
      
      ;; Call Milky Flat
      flat = hires_pixflat_work(hires[gdflt].img_ov, $
                                GAIN=hires[gdflt[0]].gain, $
                                RN=hires[gdflt[0]].readno, $
                                flatstddev=sig_flat)
         
      var = sig_flat^2 
      ivar = 1./(var + (var EQ 0))
      nrow = (size(ivar))[2]

      ;; Zero out edge of blue chip
      if qq EQ 1 and not keyword_set( ALLBLUE ) then begin
          lrow = round(1970./hires[gdflt[0]].colbin) 
          ivar[lrow:*,*] = -1
      endif

;      ivar[*,0] =  0.
;      flat[*,0] =  0.
;      ivar[*,nrow-1] = 0.
;      ivar[*,nrow-1] = 0.0

      ;; Output
      mwrfits, flat, outfil, head, /create, /silent
      mwrfits, ivar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'hires_pixflat: Flat created ', outfil+'.gz'
      
      if not keyword_set(SVOV) then hires_delov, hires, gdflt, /silent

  endfor

  print, 'hires_pixflat: All done!'
  return
end
