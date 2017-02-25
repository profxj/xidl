;+
; NAME:
;   automated_offset
;
; PURPOSE:
; Procedure that find offsets to aligne different frames
;
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
; The offset to apply
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
;           2011  Written by KK
;    05-Jun-2012  Revised by MF
;
;-
;------------------------------------------------------------------------------
;; 


pro automated_offset, files, outx, outy, offset_cut=offset_cut


  ;;some control parameters (defualt to slotis) 
  aperature = 15           ;;Number of pixels for aperature photometry
  skyannulus = [20,30]     ;;Pixels describing sky annulus
  n_brtstr = 50            ;;Set the number of brightest stars to match
  if not keyword_set(offset_cut) then offset_cut = 20 ;;Cutoff in pixels for offsets, anything greater is ignored
  binsize = 0.5            ;;Offset bin size
  
  
  ;;some space
  nfiles = N_ELEMENTS(files)
  x_offset = MAKE_ARRAY(nfiles, /FLOAT)
  y_offset = MAKE_ARRAY(nfiles, /FLOAT)
  

  ;;Find 50 brighest stars for the first image
  first_image = MRDFITS(files[0], 0, img_hdr)
  gain = 1.; FXPAR(img_hdr, 'GAIN')  ;; gain applied
  ;;Find all objects in first image
  find, first_image, first_x, first_y, flux, sharp, round, 0.5, 10.0, [-1.0,1.0], [0.2,1.0], /SILENT 
  ;;Perform aperature phot. on all objects first img
  aper, first_image, first_x, first_y, first_mag, first_stddev, sky, skyerr, gain, aperature, $
        skyannulus, [-5000,5000], /SILENT 
  ;;Find brightest stars in template
  first_brtstr = sort(first_mag) 
  first_brtstr = first_brtstr[0:n_brtstr-1]
  
  ;;Loop finds offset for each of the extra images compared to the first image
  for i=1, nfiles-1 do begin
	image = mrdfits(files[i], 0, img_hdr)
	gain = 1. ;fxpar(img_hdr, 'GAIN') '' gain applied
	FIND, image, x, y, flux, sharp, round, 0.5, 10.0, [-1.0,1.0], [0.2,1.0], /SILENT 
	APER, image, x, y, mag, stddev, sky, skyerr, gain, aperature, skyannulus, [-5000,5000], /SILENT 
        ;;Find brightest stars in template
        brtstr = SORT(mag) 
        brtstr = brtstr[0:n_brtstr-1]
        ;;Make arrays that stores the the differences in x,y between current image stars and first image stars
        delta_x = MAKE_ARRAY(n_brtstr, n_brtstr, /FLOAT) 
	delta_y = MAKE_ARRAY(n_brtstr, n_brtstr, /FLOAT)
        ;;Find delta x,y for the brightest image stars for all the brightest template stars
        for j=0, n_brtstr-1 do begin  
		delta_x[j,*] = x[brtstr[j]] - first_x[first_brtstr]
		delta_y[j,*] = y[brtstr[j]] - first_y[first_brtstr]
	endfor	
        ;;Cut all offsets greater than the set cut (in pixels)
	delta_x = delta_x[WHERE(ABS(delta_x) LE offset_cut)] 
	delta_y = delta_y[WHERE(ABS(delta_y) LE offset_cut)]
        ;;Bin up delta x,y because the largest bin contains the true x,y offsets of the image from the template
	delta_xbins = HISTOGRAM(delta_x, BINSIZE=binsize, LOCATION=delta_xbin_values)
	delta_ybins = HISTOGRAM(delta_y, BINSIZE=binsize, LOCATION=delta_ybin_values)
	max_xbin_value = delta_xbin_values[WHERE(delta_xbins EQ MAX(delta_xbins))] 
        ;;Find the location of the largest bins
	max_ybin_value = delta_ybin_values[WHERE(delta_ybins EQ MAX(delta_ybins))]
        ;;Store indicies in delta x,y found only in the largets bins
	xbin_points = WHERE(delta_x GE max_xbin_value[0] AND delta_x LE max_xbin_value[0]+binsize) 
        ybin_points = WHERE(delta_y GE max_ybin_value[0] AND delta_y LE max_ybin_value[0]+binsize)
        ;;Take median of largest bins to be the x,y offsets, should be VERY robust!
	x_offset[i] = MEDIAN(delta_x[xbin_points]) 
	y_offset[i] = MEDIAN(delta_y[ybin_points])
ENDFOR

outx = -x_offset
outy = -y_offset

END
