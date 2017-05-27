;+
; NAME:
;   spslit_combine
;
; PURPOSE:
;   combines the separate exposures of a slitlet into individual sky
;   subtracted 2-d spectra, writes the combined file back as net FITS
;   table.
;
; CALLING SEQUENCE:
;    spslit_combine, filelist, [result]
; 
; INPUTS:
;    filelist -- array of filenames containing spSlit.*.fits* data
;
; OPTIONAL INPUTS:
;	
; KEYWORDS:
;
;    NLSKY -- if set, nonlocal sky is put in HDU 3
; OUTPUTS:
;    a slit.*.fits file is written for each input file
;    [result] -- last structure computed is returned (could be only
;                one called)
; EXAMPLES:
;
; COMMENTS:
;    calls slitexam
;    output file is sky-subtracted, weighted average file based on
;    separate frames contained within the spSlit file
;
; REVISION HISTORY:
;   09jun02 MD  
;   05jun02 MD added WCS keywords for Lambda
;----------------------------------------------------------------------

; At this stage we should probably be extracting the OBJECTS
; themselves and writing these out as spObj files...   - DPF


pro spslit_combine, filelist,  result, nlsky = nlsky

if n_elements(nlsky) eq 0 then nlsky = 0


  for ifile=0, n_elements(filelist)-1 do begin

     filename = strsplit(filelist[ifile], '.', /extract) 
     chars = strlen(filename[0]) ;get length of string up to first '.'
;split file into pieces to reconstruct new file name
     type = strmid(filename[0], n_elements(chars)-6, 6)
     if type ne 'spSlit' then message, 'wrong type of file: '+filelist[ifile] $
     else begin
        prefix =  (chars le 6) ? '' : strmid(filename[0], 0, chars-7)  
     ;   match = strmatch(filename[0], 'spSlit')
       header = headfits(filelist[ifile], ext=1)

;       errmsg = '' ; to avoid crashes
;       fxhclean, header, err=errmsg ; remove BINTABLE info from structure
; copy thr header from the spSlit file...change the extname to 'slit'
; from 'spSlit'.
       header = copy_header(header, 'slit')

; this is the workhorse of combing frames 
       result = 0 ;reset  before calling
       slitexam, filelist[ifile], result, lineparams=lineparams, /noplot, $
         totalexp=totalexp
       fileout = prefix + 'slit.' 
 
       i = 0
       while i lt 10 do begin 
;construct remainder of filename
          i = i+1
          if filename[i] ne 'fits' then fileout = fileout + filename[i] +'.' $
          else begin
             fileout = fileout+ 'fits'
             i = 100
          endelse
       endwhile
     endelse

     sxaddpar, header, 'EXPTIME', totalexp,'Total exposure time'
     fxaddpar, header, 'COMMENT', 'spslit frames combined, sky subtracted'
     fxaddpar, header, 'COMMENT', 'CR removal by means of chisqr analysis'

;add WCS information for lambda
; print, 'Making the header ...'
;  mkhdr, hdr, image, /extend

  tags = tag_names(result)
  if total(tags eq 'TILTX') then begin 
     ll = polyleg([-1., 1.], result.lambdax)
     delt=[(ll[1]-ll[0])/4096,0]
     CRPIX=[0,0]
     CRVAL = [ll[0], 0]
  endif 
  if total(tags eq 'COEFF') then begin 
     delt = result.coeff[1,2]/2048.  ; Ang per pixel
     CRPIX=[2048,0]
     CRVAL = [result.coeff[0,2], 0] ; row 1
  endif 


  CTYPE=['LAMBDA','LAMBDA']
  make_astr,astr,delt=delt,crpix=CRPIX,crval=CRVAL,ctype=CTYPE
  putast,header,astr

     dimens = (size(lineparams, /dim))
     ndim = n_elements(dimens)
     nrow = dimens[0]
     if ndim gt 1 then nline = dimens[1] else nline = 1

       if nline lt 2 then begin
              lineparams = {amp: 0., cen: 0., sig: 0., base: 0., mask: 0B}
              print, 'unable to perform skytweak in file: ', fileout
              message, 'unable to perform skytweak in  file: '+fileout, /INFO
              ; 2^4 bit of infomask: no skytweak
              result.infomask=result.infomask + 16b
        endif


; formerly:if  this slit is really bad, don't output the file

        mwrfits, result, fileout, header, /create, /no_comment
        print, 'writing output file: ', fileout

 

;     if nline ge 2 then begin

; append lineparams structure containing params of gauss fit to sky lines

        h2 = header
        fxhclean, h2
        sxaddpar, h2, 'NROW', nrow, ' Number of rows'
        sxaddpar, h2, 'NLINE', nline, ' Number of lines'
        mwrfits, lineparams, fileout, h2

;     endif 
                               ; Now, do the nonlocal sky fits and
                                ; put the resulting combined frame in
                                ; HDU 3

        if nlsky then begin
          slitexam, filelist[ifile], resultnl, /NLSKY, /noplot
          if n_elements(resultnl.slitfn) gt 1 then $
            mwrfits, resultnl, fileout, header, /no_comment
        endif

; save gzipping for end of entire process
;     spawn, 'gzip ' + fileout +' &'
;     endif else begin
     
;     endelse

  endfor

return
end





