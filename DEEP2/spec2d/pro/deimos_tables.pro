function  deimos_tables,  filename,  bluslits=bluslits, slitobjmap=slitobjmap
;+
; NAME:
;   deimos_tables
;
; PURPOSE:
;   reads in tables of object catalogs and slitmask coordinates for a
;   given datafile. DEIMOS fits headers are very extensive, and this
;   program reads the relevant FITS tables needed to associate objects
;   and slitlets.
;
; CATEGORY:
;   spec2d
;
; CALLING SEQUENCE:
;   catalog = deimos_tables( filename, bluslits=bluslits,slitobjmap=slitobjmap)
;
; INPUTS:
;   filename  -- string name of DEIMOS fits data file
;
; KEYWORD PARAMETERS:

;
; OUTPUTS:
;   catalog  -- structure of object list in this mask
;
;
; OPTIONAL OUTPUTS:
;   bluslits  -- structure of slitmask design (x,y coords of center of
;                each slitlet, plus length of slitlet)
;   slitobjmap -- slitobjmap structure
;
; MODIFICATION HISTORY:
;  md 11feb02
;-

; -------- check inputs
  if NOT keyword_set(filename) then message, 'Please pass a filename'

  foo = findfile(filename, count=filect)
  if filect EQ 0 then begin 
     print, 'FILE NOT FOUND:'
     message, filename
  endif 
  fits_info,  filename, /silent, n_ext=n_ext ;get number of table extensions
  if n_ext eq 0 then begin 
     print,  'ERROR: no extensions in this file!'
     return, 0
  endif 

; -------- determine if topdist and botdist are correct. that is,
;          catch possible error in mask design. we must determine the
;          creation date of the mask so as to adjust for a bug in
;          Drew's mask code. Only do this if a science frame is passed
;          to the routine...this problem has already been corrected in
;          the bin tables.
  isabintab = total(strpos(findfile(filename), 'bintab')) gt 0
; just to be safe, initialize the ddumb variable.
      ddumb = 0

  if not(isabintab) then begin
      fits_open, filename, fcb


; determine which extensions are bin tables.
      dex = where(fcb.xtension eq 'BINTABLE' or $
                  fcb.xtension eq 'TABLE', tabcnt)
; find mask design bin table.  
      maskdex = where(fcb.extname eq 'MaskDesign', cnt)
      if cnt eq 0 then begin 
          fits_close, fcb
          message, 'MaskDesign bin table not found!' 
      endif else maskdex = maskdex[0]
; read in mask design table.
      masktable = mrdfits(filename, maskdex, /silent)
; get mask design data and then split it into year, month, and day.
      date = (strsplit(masktable.DESDATE, 'T', /extract))[0] ;date has form 'yyyy-mm-dd'
      subdate = strsplit(date, '-', /extract)
      print, '(deimos_tables.pro) Printing mask design date: ' + $
        strcompress(date, /remove_all)
; catch possible error.
      if n_elements(subdate) ne 3 then begin
          print, '(deimos_tables.pro) ERROR: error in extraction of mask design date!!' 
          ddumb = 0
      endif else begin
; all masks made between the dates July 27, 2002 and May 15,
; 2003 have their botdist and topdist mixed up!
          if ( long(subdate[0]) eq 2002 and long(subdate[1]) ge 7 ) or $
            ( long(subdate[0]) ge 2003 ) then begin
              if long(subdate[0]) eq 2002 and $
                long(subdate[1]) eq 7 and $
                long(subdate[2]) lt 25 then ddumb = 0 $
              else ddumb = 1
          endif else ddumb = 0
      endelse
; close the fcb.
      fits_close, fcb
  endif

  catflag = 0
  bluflag = 0
  desiflag = 0

  minext = 9
  isabintab = total(strpos(findfile(filename), 'bintab')) gt 0
  if isabintab then minext = 1

  for i=minext, n_ext  do begin ;read only table extensions
     header = headfits(filename,  exten=i)
     type = fxpar(header, 'XTENSION')
     extname = fxpar(header, 'EXTNAME')
     if (type eq 'IMAGE   ') then extname = 'Image'

     table = (type eq 'TABLE   ') OR (type eq 'BINTABLE')
     if table and extname eq 'ObjectCat' then begin
        catalog = mrdfits(filename, i,  /silent)
        catflag = 1
     endif
     if table and extname eq 'BluSlits' then begin
       blu_in = mrdfits(filename, i,  /silent)
       temp = {slitno:0L, dslitid:0L, xmm:0.0,  ymm:0.0, $
               pa:0.0, slitlength:0.0,  slitwidth:0.0, $
              xmmb:0.0,ymmb:0.0,xmmt:0.0,ymmt:0.0}
       bluslits = replicate(temp, n_elements(blu_in) )

       bluslits.dslitid = blu_in.dslitid
       bluslits.xmm = (blu_in.slitx1 + blu_in.slitx2 + $
                       blu_in.slitx3 + blu_in.slitx4)/4.
       bluslits.ymm = (blu_in.slity1 + blu_in.slity2 + $
                       blu_in.slity3 + blu_in.slity4 )/4.
       bluslits.xmmt=(blu_in.slitx1 + blu_in.slitx4)/2.
       bluslits.ymmt=(blu_in.slity1 + blu_in.slity4)/2.
       bluslits.xmmb=(blu_in.slitx2 + blu_in.slitx3)/2.
       bluslits.ymmb=(blu_in.slity2 + blu_in.slity3)/2.


       bluslits.pa  = atan(blu_in.slity2 - blu_in.slity1, $
                           blu_in.slitx2 - blu_in.slitx1)*!radeg


       bluslits.pa =  bluslits.pa*(abs(bluslits.pa) lt 90.) + $
                      (bluslits.pa+180.)*(bluslits.pa lt -90.) + $
                      (bluslits.pa-180.)*(bluslits.pa gt  90.)
       bluslits.slitlength = abs(blu_in.slitx1 - blu_in.slitx2)
       bluslits.slitwidth  = float(abs(blu_in.slity2 - blu_in.slity3))
       bluslits = bluslits[sort(bluslits.dslitid)] ;sort on dslitid
       bluflag = 1
     endif
     if table and extname eq 'SlitObjMap' then begin
        if arg_present(slitobjmap) then begin
          slitobjmap = mrdfits(filename, i, /silent)
          if ddumb then begin
            print, '(deimos_tables.pro) Correcting error in mask design...'
            print, '(deimos_tables.pro) Switching topdist and botdist.'
            tmpdist = slitobjmap.topdist
            slitobjmap.topdist = slitobjmap.botdist
            slitobjmap.botdist = tmpdist
          endif
        endif
     endif
     if table and extname eq 'DesiSlits' then begin
        desislits = mrdfits(filename, i, /silent)
        desislits = desislits[sort(desislits.dslitid)] ;sort 
        desiflag = 1
     endif
     
     
  endfor
  if catflag eq 0 then print, 'ERROR: no object catalog in this file'
  if bluflag eq 0 then print, 'ERROR: no blueprint table in this file'
  if desiflag eq 0 then print, 'ERROR: no DesiSlits table in this file'
; match slitname to UCB slitname convention, not SYBASE tag
  bluslits.slitno =  long(desislits.slitname) ;copy into bluslits
  
  return,  catalog 
end



