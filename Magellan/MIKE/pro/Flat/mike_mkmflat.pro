;+ 
; NAME:
; mike_mkmflat   
;     Version 1.1
;
; PURPOSE:
;    This set of routines takes a series of flats observed through the
;    diffuser and creates a normalized Flat used to correct
;    pixel-to-pixel response variations.  A principal challenge with
;    MIKE in its current state (April 2004) is that it is difficult to
;    get sufficient counts on the blue side.  
;
;    mike_mkmflat :: The main routine simply does some basic
;    accounting, organizes bias subtraction and performs I/O.
;
;    mike_mkflat_work :: Accepts the name(s) of a series of OV
;    subtracted milky flats.  It then:
;      1.  Opens the file and creates an inverse variance array
;      2.  Takes out the low order variation in the image (lo_interp)
;      3.  Performs a series of medians along the columns rejecting bad
;          pixels (replace by local mean). 
;      4.  If multiple images were input, they are stacked with
;          rejection and the final image is returned
;
; CALLING SEQUENCE:
;   
;  mike_mkmflat, mike, setup, [side]
;
; INPUTS:
;   mike     -  MIKE structure
;   setup    -  Setup identifier 
;   [side]   -  Blue (1), Red (2), or both [1,2L]    (Default: [1,2L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per setup per side with names like
;  'Flats/Flat_B_01_M.fits.gz' 
;
; OPTIONAL KEYWORDS:
;   /CLOBBER - Overwrite Output MilkyFlat
;   /USEBIAS - Use bias frame in OV subtraction
;   /SMOOTH  - Smooth the Flat several times over.  This is only
;              recommended for milky flats of objects with many
;              absorption features (which are not recommended).
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_mkmflat, mike, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;   mike_getfil
;   mike_subbias
;   mike_mkflat_work
;   x_statarray
;
; REVISION HISTORY:
;   16-May-2003 Adapted by JXP from existing programs by SB
;   24-Feb-2004 Switched to a series of median/linear interpolations (SB)
;   22-Jun-2005 Returned to simple median smooth as the default
;-
;------------------------------------------------------------------------------

function lo_interp, image
     
     ncol = (size(image))[1]
     nrow = (size(image))[2]

     colboxsize = long(8*ncol/2048) > 4
     rowboxsize = long(4*nrow/1024) > 2
     newcol = ncol/colboxsize
     newrow = nrow/rowboxsize

     subcol = newcol * colboxsize
     subrow = newrow * rowboxsize

     repack = reform(image[0:subcol-1,0:subrow-1], colboxsize, newcol, $
                        rowboxsize, newrow)
     repack = transpose(repack, [1, 0, 2, 3])
     small_img = djs_avsigclip(reform(repack, newcol, $
                colboxsize*rowboxsize, newrow),2)

     xvalue = (findgen(ncol)-(colboxsize-1)/2)/colboxsize
     yvalue = (findgen(nrow)-(rowboxsize-1)/2)/rowboxsize
     result = bilinear(small_img, xvalue, yvalue)

     return, result
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

function mike_mkflat_work, milkyfile, flatstddev=flatstddev, $
  GAIN=gain, RN=rn, SMOOTH=smooth
   
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'flat = mike_mkflat_work(file, GAIN=, RN=) [v1.1] '
      return, -1
  endif 

     if NOT keyword_set(rn) then rn = 4.1
     if NOT keyword_set(gain) then gain = 1.06
     if NOT keyword_set(fitsout) then fitsout=0

     nfiles = n_elements(milkyfile)

     if nfiles GT 1 then begin
         print, 'Doing   1st of ', nfiles, format='(a,i3,$)'
         first = mike_mkflat_work(milkyfile[0], GAIN=gain, rn=rn) 
         nx = (size(first))[1]
         ny = (size(first))[2]
         final = fltarr(nx,ny, nfiles)
         final[*,*,0] = first
         for i=1,nfiles -1 do begin
           print, 'Doing ', i+1, (i GT 2) ? 'th of ': (i EQ 1) ? $
                   'nd of ' : 'rd of ',  nfiles, format='(a,i3,a,i3,$)'
           final[*,*,i] = mike_mkflat_work(milkyfile[i], GAIN=gain, rn=rn, $
                                          SMOOTH=smooth)
         endfor
         print, 'mike_mkmflat: Final Combine'
         x_statarray, final, 3, mean=flat, stddev=flatstddev, sigrej=2., $
           /OVRIDE

         ; factor of 2 kludge due to gross underestimation
         flatstddev = 2.0 * flatstddev / sqrt(nfiles)
         return, flat
     endif

     
     ;; Open file
     milky = xmrdfits(milkyfile, /silent)

     milky = milky*gain
     milkyivar =  1./(abs(milky)*gain + rn^2)

     ; Transpose for quicker execution
 
     image  = transpose(milky)
     ncol = (size(image,/dimen))[0]
     nrow = (size(image,/dimen))[1]

     invvar = transpose(milkyivar)

     if keyword_set(SMOOTH) then begin

         width = 2*(ncol/1024L) + 1L
         colw  = 2*long(2.5*width) + 1L
;
;     Step 1, take out broad order trend to put median on even footing
;
         lo_interp = lo_interp(median(image,width))
         image_norm = image/(lo_interp + (image*(lo_interp EQ 0)))
         ivar_norm = invvar*lo_interp^2 


         colprebuffer = image_norm[reverse(findgen(colw/2))+1,*]
         colpostbuffer = image_norm[ncol - 2 - reverse(findgen(colw/2)),*]
         colimage = [colprebuffer,image_norm,colpostbuffer]

;
;     Step 2, Take out residuals running down order with smooth along orders
;
         col_smooth = median_row(colimage, colw) > 0
         col2 = colimage / (col_smooth + (colimage > 1) * (col_smooth EQ 0))  

   
         rowwidth = width 
         rowprebuffer = col2[*,reverse(findgen(rowwidth/2))+1] 
         rowpostbuffer = col2[*,nrow - 2 - reverse(findgen(rowwidth/2))] 

         workimage = [[rowprebuffer],[col2],[rowpostbuffer]]
         workivar =[[fltarr((colw-1)+ncol,rowwidth/2)], $
                    [fltarr(colw/2,nrow), ivar_norm, fltarr(colw/2,nrow)], $
                    [fltarr((colw-1)+ncol,rowwidth/2)]]

;
;     Step 3, First iteration with median filter of box size (width x width)
;
         firstmed = median(workimage, width)
         ratio = workimage/(firstmed + (firstmed EQ 0)) * (firstmed GT 0)
         absresidual = abs(ratio - 1)
         maskratio = (absresidual*sqrt(workivar) GT 4) AND $
           (absresidual GT 0.03)

         masksmooth = (smooth(1.0*maskratio,width) GT 0.2) OR maskratio
         replace = where(masksmooth, nreplace, complement=good)
         print, '....Masking ', nreplace, ' Pixels', format='(a,i6,a, $)'
         if nreplace GT 0 then begin
;
;     Step 4, 2nd iteration with median filter of box size (width x width)
;             Replace masked pixels with local mean
;
             swidth = 2*width+1
             denom = smooth(1.0-masksmooth, swidth)
             smoothreplace = smooth((1.0-masksmooth) * workimage,swidth) / $
               (denom > 1./swidth^2) * (denom GT (2./swidth)^2) + $
               (denom LE (2./swidth)^2)
               
             workimage2 = workimage
             workimage2[replace] = smoothreplace[replace]
             firstmed = median(workimage2, width)
         endif
         ;;
         ;;    Construct final flat field and compute some stats
         ;;
         flat = workimage/(firstmed + (firstmed EQ 0)) * (firstmed GT 0)
         flat = transpose((flat * $
                           (flat GT 0.1))[colw/2:colw/2+ncol-1, $
                                          rowwidth/2:rowwidth/2+nrow-1]) > 0.1
         print, mean(abs(1.-flat[good])), mean(flat[good]), $
           format='(f9.5, f9.5)'
     endif else begin
         if not keyword_set(WIDTH) then begin
             sz = size(milky, /dimensions)
             cbin = round(2048. / sz[0])
             rbin = round(4096. / sz[1])
             width = round(62./rbin)
         endif
         norm = transpose( x_medianrow(image,width) )
         flat = milky/norm
         good = where(flat GT 0.)
         print, mean(abs(1.-flat[good])), mean(flat[good]), $
           format='(f9.5, f9.5)'
     endelse


     return, flat

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_mkmflat, mike, setup, side, CLOBBER=clobber, USEBIAS=usebias, $
                  SMOOTH=smooth, _EXTRA=extra

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'mike_mkmflat, mike, setup, [side], ' + $
        ' /USEBIAS, /CLOBBER, /SMOOTH [v2.0]'
      return
  endif 
  
  ;; QA
  if setup LT 10 then fqa = 'QA/Flats0'+strtrim(setup,2) $
  else fqa = 'QA/Flats'+strtrim(setup,2)
  a = findfile(fqa, count=count)
  if count EQ 0 then file_mkdir, fqa

;  Optional Keywords
  if not keyword_set( SIDE ) then side = [1L,2L]
  
  for ii=0L,n_elements(side)-1 do begin
      qq = side[ii]
      if qq EQ 1 then print, 'mike_mkmflat: Creating BLUE Milky flat' $
      else print, 'mike_mkmflat: Creating RED Milky flat'

      ;; Outfil
      outfil = mike_getfil('mflat_fil', setup, SIDE=qq, /name, CHKFIL=chkf) 
      if CHKF NE 0 AND not keyword_set( CLOBBER ) then begin
          print, 'mike_mkmflat: Milky flat exists, moving on..'
          continue
      endif

      ;; Grab flats
      gdflt = where(mike.side EQ qq AND mike.flg_anly NE 0 AND $
                    strtrim(mike.type,2) EQ 'MFLT' AND $
                    mike.setup EQ setup, nflt)
      if nflt EQ 0 then begin
          print, 'mike_mkmflat: No Flats of type MFLT found!' 
          return
      endif

      ;; Bias Subtract
      mike_subbias, mike, gdflt, CLOBBER=ovclob, USEBIAS=usebias, HEAD=head, $
                    _EXTRA=extra
                    
      
      ;; Call Milky Flat
      flat = mike_mkflat_work(mike[gdflt].img_ov, GAIN=mike[gdflt[0]].gain, $
                              RN=mike[gdflt[0]].readno, flatstddev=sig_flat, $
                              SMOOTH=smooth)
         
      if not keyword_set(sig_flat) then begin
          ;; Single flat!!
          print, 'mike_mkmflat:  You used only one milky flat!!'
          var = xmrdfits(mike[gdflt].img_ov,/sile)
          sig_flat = var
          sig_flat[*] = 1.
      endif else var = sig_flat^2 

      ivar = 1./(var + (var EQ 0)) * (sig_flat GT 0)
      nrow = (size(ivar))[2]

      ivar[*,0] =  0.
      flat[*,0] =  0.
      ivar[*,nrow-1] = 0.
      ivar[*,nrow-1] = 0.0

      ;; Output
      mike_taghead, head
      mwrfits, flat, outfil, head, /create, /silent
      mwrfits, ivar, outfil, /silent
      spawn, 'gzip -f '+outfil
      print, 'mike_mkmflat: Flat created ', outfil+'.gz'
      
      if not keyword_set(SVOV) then mike_delov, mike, gdflt

  endfor

  print, 'mike_mkmflat: All done!'
  return
end
