;+
; NAME:
;    slit_merge
;
; PURPOSE:
;    combines slit*.fits.gz files into one large fits file, for viewing
;    all of a DEIMOS mask at once. B images on left, R images on right
;    last column of returned array is float(slitn)
;
; CALLING SEQUENCE:
;    image = slit_merge(filelist)
;
; INPUTS:
;    filelist -- list of slit.*.fits files
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
; OUTPUTS:
;   image-- a full frame of the DEIMOS data
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   this is intended to be used to merge into one fits file all the 2d
;   sky subtracted slitlets for a given DEIMOS mask
;
; REVISION HISTORY:
;   09jun02 MD
;----------------------------------------------------------------------

function  slit_merge, filelist


  filelist = filelist[sort(filelist)] ;make certain list is sorted
  nrows = intarr(n_elements(filelist)) ;rows per slitlet
  imageB_cnt = 0
  imageR_cnt = 0
  slitB = 0
  slitR = 0
  imageB = fltarr(4096)
  imageR = fltarr(4096)
  bar = fltarr(4096) -1000. ;separating bar
  last_column = -1.

  for ifile=0, n_elements(filelist)-1 do begin
     data_in = mrdfits(filelist[ifile], 1, /silent )
     nrows[ifile] = (size(data_in.flux, /dimen))[1]
     print, 'reading file, rows: ', filelist[ifile], nrows[ifile]
     filename = strsplit(filelist[ifile], '.', /extract) 
     slitn = strmid(filename[2], 0, 3) ;first 3 characters, slit number
     color = strmid(filename[2], 3, 1) ; R or B?

     if color eq 'B' then begin ;B side spectrum
        if imageB_cnt eq 0 then imageB = data_in.flux else $
             imageB = [[imageB], [data_in.flux]]
        imageB = [[imageB], [bar]] ;put in separating bar
        slitb = slitn

        if imageB_cnt gt imageR_cnt then begin ;R was missing on previous
           imageR = [[imageR], [fltarr(4096, nrows[ifile-1])]] ;blank R
           imageR_cnt = imageR_cnt + 1 ;catch up the count
        endif
        imageB_cnt = imageB_cnt + 1
     endif

     if color eq 'R' then begin ;R side spectrum
        if imageR_cnt eq 0 then imageR = data_in.flux else $
             imageR = [[imageR], [data_in.flux]]
        if slitb ne slitn then begin $ ;B was missing on this frame
           imageB = [[imageB], [fltarr(4096, nrows[ifile])]] ;blank B
           imageB_cnt = imageB_cnt + 1 ;catch up the count
        endif
        imageR = [[imageR], [bar]] ;put in separating bar
        slitr = slitn
        imageR_cnt = imageR_cnt + 1

;check to ensure that number of R and B rows match
         sizer = (size(imager, /dimen))[1]
         sizeb = (size(imageb, /dimen))[1]
         if sizer gt sizeb then $   ;more R than B
           for j=0, sizer-sizeb -1 do $
              imageB = [[imageB], [bar]] ;make more separation bars
         if sizeb gt sizer then $ ;more B than R
            for j=0, sizeb-sizer -1 do $
              imageR = [[imageR], [bar]] ;make more separation bars

;add last column, which specifies slitn
        append = (nrows[ifile] >  nrows[ifile-1])+1 
        if last_column[0] eq -1. then last_column = $
                   fltarr(append) + fix(slitb) else $
           last_column = [last_column, fltarr(append) + fix(slitb) ]
;extend last column with appropriate slit number
     endif

  endfor

  image = [imageB, imageR] ;merge side by side

;TBD fix this later
;  image = [image, transpose(last_column)] ;merge column of slit id

return,  image
end


