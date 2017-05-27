;+
;
; NAME
;      find_objpos.pro
;
; PURPOSE
;      Determine the position (pixel column) of an object in the slit
;      from the design specifications of the slit mask. The function
;      find_objpos takes two arguments, a slit number and the number
;      of pixel rows (spatially) in the slit. the function returns the
;      location of the object in the slit in terms of the pixel row
;      (as a float) of the object relative to the slit edge. The
;      information needed to determine the object position for a given
;      slit is drawn from one of three possible areas: (1) the
;      xxxx.bintab.fits file created in the reduction process by the
;      procedure make_bintab_file.pro, (2) the files 1HSmask.xxxx.fits
;      that specify the design specs for the mask and all the slits on
;      the mask (slitinfo and maskinfo), and lastly (3) a raw DEIMOS
;      frame and the appropriate pcat structure. Option (1) is the
;      default setting, but option (2) can be activated by passing the
;      maskinfo and slitinfo tables to the routine and option (3) can
;      be selected by passing a raw DEIMOS FITS file and the pcat to
;      the routine.
;
;
; SYNTAX
;      objpos = find_objpos(slitn, nrows, [fwhm=fwhm, bfile=bfile, $
;                                   objnum=objnum, objtype=objtype]
; INPUTS
;      slitn = the number of a slit (not an array of slit numbers) from
;              a single slit mask. This slit number can be easily read
;              from the header of any slit file.
;      nrows = the number of pixel rows in the spatial direction of a
;              slit. This value is used to determine how well the slit
;              was traced out by the pipeline (this defines nrows)
;              relative to the length of the slit as expected from the
;              design specs. In determining the object position this
;              is an important factor.
;      fwhm = an optional variable input that will be set equal to
;             the fwhm of the spatial profile of the object. This
;             width is predicted using the design specifications of
;             the slit mask and source information from the DEEP
;             photometry catalog. If multiple slit numbers are passed,
;             then an array of fwhm values will be returned.
;      bfile  = path/filename of a bintab file.  File is found in
;      current directory if not set.
;      objnum = an optional variable argument that will be set equal
;               to the object number(s) for the objects found in slit
;               number slitn. Recall that the object number is of the
;               form ppxxxxxx where pp designates the photometry
;               field. Also note that if the first digit of the xxxxxx
;               number is a 5 (or not zero), then the slit is a sky
;               slit. objnum will be returned as a STRING.
;      objtype = DEIMOS type of the object (optional)
;
; KEYWORDS
;      None.
;
; OUTPUTS
;      objpos = the position of the object in the slit. 
;      fwhm (see INPUTS)
;
; PROCEDURES CALLED 
;      FITS_INFO
;      HEADFITS
;      MRDFITS
;
; EXAMPLES
;      
;
; COMMENTS
;
; HISTORY
;      Created July 19, 2002 by mcc. 
;      Revised July 27. 2002 by mcc - employed information from the
;         pcat structure to determine the object position. Function
;         now takes into account the position angle and shape of the
;         object as well as the orientation of the mask.
;      Revised August 26, 2002 by mcc - revised routine to handle
;         cases when there are two objects in a single slit. In these
;         cases, the returned fwhm and pos variables are vectors of
;         length equal to the number of objects in the slit. The
;         objects are sorted according to magnitude (bright-to-faint).
;      Revised March 2, 2002 by mcc - major revision of
;         routine. Removed the capability of using the 1HSmask file or
;         DEIMOS rawframe to get the design positions...now relies on
;         existence of bintab file. Also, added code so that the
;         seeing from the pcat is used to put a lower limit on the
;         extent (fwhm) of the object.
;-

; make a simple procedure to set the return variables to null
; results in case of errors.
pro set_nulls, fwhm=fwhm, objnum=objnum, objtype=objtype
  fwhm = 0.0
  objnum = 0
  objtype = ''
end
; ---------------------------------

function find_objpos, slitn, nrows, fwhm=fwhm, bfile=bfile,$
                      objnum=objnum, objtype=objtype

; check that the user passed the needed parameters.
  if n_elements(slitn) eq 0 then begin
      print, 'You must supply a slit number!'
      message, 'See calling sequence in DOC_LIBRARY info.'
  endif
  if n_elements(nrows) eq 0 then begin
      print, 'You must supply the number of pixel rows in slit!'
      message, 'See calling sequence in DOC_LIBRARY info.'
  endif

; convert slitn from any input type to type long.
  slitn = long(slitn)

; find the bintab file.
  if n_elements(bfile) ne 0 then nfile =1 $
    else bfile = findfile('*.bintabs.fits*', count=nfile)
  if nfile eq 0 then $
    message, 'No bintable file found in current directory!!' $
  else begin
; gunzip the bintab file if needed...this is not needed if the user is
; running the most recent version of the GODDARD FITS I/O IDL library.
      bfile = bfile[0]
      len = strlen(bfile)
      if strmid(bfile, len-3, 3) eq '.gz' then begin
          spawn, 'gunzip -v1 ' + bfile
          bfile = strmid(bfile, 0, len-3)
      endif
  endelse

; open the bintab file and extract the needed tables.
  fits_open, bfile, fcb
; extract the names of the tables.
  extnames = strcompress(fcb.extname, /rem)
; remember to close the open file!
  fits_close, fcb

; ---------------------------------
; read in the slit design (DesiSlits) table and sort it so that the
; entries will match-up correctly with the entrees in the ObjectCat
; table.
  desidex = where(extnames eq 'DesiSlits', desicnt)
  if desicnt gt 0 then design = mrdfits(bfile, desidex[0], /silent) $
; if the design table is NOT found, then return null values.
  else begin
      print, 'No DesiSlits table found in ' + bfile + ' file!!!'
      set_nulls, fwhm=fwhm, objnum=objnum, objtype=objtype
      return, 0.
  endelse
; read in the slit object map (SlitObjMap) table.
  slitdex = where(extnames eq 'SlitObjMap', slitcnt)
  if slitcnt gt 0 then slitobjmap = mrdfits(bfile, slitdex[0], /silent)
; read in the mask design (MaskDesign) table.
  maskdex = where(extnames eq 'MaskDesign', maskcnt)
  if maskcnt gt 0 then mbasics = mrdfits(bfile, maskdex[0], /silent)
; check if this is a DEEP2 mask...
  deimos_isdeep, isdeep
; read in the object catalog (ObjectCat) table and sort it so that the
; entries will match-up correctly with the entrees in the DesiSlits
; table.
  objdex = where(extnames eq 'ObjectCat', objcnt)
  if objcnt gt 0 then objectcat = mrdfits(bfile, objdex[0], /silent) 
; read in the pcat (PhotCat) table and its header.
  pdex = where(extnames eq 'PhotCat', pcnt)
  if pcnt gt 0 then pcat = mrdfits(bfile, pdex[0], phdr, /silent)


; ---------------------------------
; match the given slit (slitn) to the DesiSlits table by matching the
; slitn value to the slit number entry in the DesiSlits table. Note
; that design.slitname is an array of strings - ['001', '002', '003',
; etc.] - let's convert them to longs.
  slitnumbers = long(design.slitname)
  desidex = where(slitnumbers eq slitn, cnt)
; if slit is not found in the DesiSlits table then return an error
; message and null values for the output parameters.
  if cnt eq 0 then begin
      print, 'Slit(s) ' + string(slitn, format='(I3.3)') + $
        ' not found in DesiSlits table!'
      set_nulls, fwhm=fwhm, objnum=objnum, objtype=objtype
      return, 0.
  endif
  objtype = design[desidex].slittyp
; ---------------------------------
; find the index in the SlitObjMap table for the given
; slit. slitobjmap.dslitid and design.dslitid are arrays of longs such
; as [7391, 7392, 7393, etc.]. I do NOT know the significance of these
; numbers.
  mapdex = where(slitobjmap.dslitid eq design[desidex].dslitid, cnt)
; if slitlet is NOT found in the SlitObjMap table, then return and
; error message and the requisite null values. Note that if this is a
; DEEP2 mask, do NOT return the objtype parameter. This value is worng
; in the cases of sky-only slits and the code is not able to check if
; this slit is a sky-only slit. 
  if cnt eq 0 then begin
      print, 'Slit(s) ' + string(slitn, format='(I3.3)') + $ 
        ' not found in SlitObjMap table!'
      if isdeep then $
        set_nulls, fwhm=fwhm, objnum=objnum, objtype=objtype $
      else begin
          set_nulls, fwhm=fwhm, objnum=objnum
      endelse
      return, 0.
  endif

; calculate the position of the object from the design specs of the
; mask. define the pixel scale for DEIMOS (arcseconds per pixel).
  pixscl = 0.117371
; get the length of the slit (in pixels) as given in the design files.
  desilen = (slitobjmap[mapdex].topdist + $
             slitobjmap[mapdex].botdist) / pixscl
; check how this design length compares to the actual measured length
; as traced out in the pipeline.
  edgeloss = (desilen - nrows) / 2.
; finally, calculate the object position.
  objpos = slitobjmap[mapdex].botdist / pixscl - edgeloss

; ---------------------------------
; now find the correct entry in the ObjectCat for this slit. Note that
; slitobjmap.objectid and objectcat.objectid are arrays of longs such
; as [11119, 11120, 11121,...]. I don't know the significance of these
; numbers.
  if cnt eq 1 then $
    objdex = where(objectcat.objectid eq $
                   slitobjmap[mapdex].objectid, cnt) $
  else begin
; if there is >1 object in the slit, then find (and return) the
; objpos, fwhm, objnum, and objtype for all objects.
      print, '(find_objpos.pro) ' + strcompress(string(cnt), /rem) + $
        ' objects in slit ' + string(slitn, format='(I3.3)') + '!!!'
      num = cnt
      cnt = 0
      objdex = intarr(num)
      for i=0,num-1 do begin
          objdex[i] = $
            where(objectcat.objectid eq $
                  slitobjmap[mapdex[i]].objectid, any)
          cnt = cnt + any
      endfor
; sort objects according to R magnitude.
;      if cnt gt 0 then begin 
;        gdex = where(objdex ne -1)
;        objdex = objdex[gdex]
;        sub = sort(objectcat[objdex].mag)
;        objdex = objdex[sub]
;        ind = where(slitobjmap[mapdex].objectid eq $
;                    objectcat[objdex].objectid)
;        mapdex = mapdex[ind]
;      endif
  endelse
; if object's slitlet is not found in the the ObjectCat, then return
; an error message and null values.
  if cnt eq 0 then begin
      print, '(find_objpos.pro) Object(s) ' + $
        strcompress(string(slitobjmap[mapdex].objectid), /rem) + $
        ' in slit(s) ' + string(slitn, format='(I3.3)') + $
        ' not found in ObjectCat table!'
      if isdeep then $
        set_nulls, fwhm=fwhm, objnum=objnum, objtype=objtype $
      else begin
          objtype = design[desidex].slittyp
          set_nulls, fwhm=fwhm, objnum=objnum
      endelse
      return, objpos
  endif

; make sure that the object numbers are strings.
  objstr = objectcat[objdex].object
  objnum = strcompress(string(objstr), /rem)

; check if this is a DEEP2 mask. if so, then check if the object number
; denotes a sky-only slitlet. for some reason, the sky-only slits are
; given a 'P' for "program object" in the SLITTYP field in the
; DesiSlit bin table.
  if isdeep then begin
      if n_elements(objnum) eq 1 then begin
          if strlen(objnum) ge 3 then thrdig = strmid(objnum, 2, 1) $
          else thrdig = ''
          if thrdig[0] eq '5' then begin 
              objtype = 'S'
; if it is a sky slit then give it a fwhm value of 0.
              fwhm = 0.
              return, objpos
          endif
      endif
  endif
; get the object type string.
  objtype = design[desidex].slittyp
  
; ---------------------------------
; finally find the index in the pcat for each object. Note that
; objectcat.object and pcat.objno are arrays of strings and longs
; respectively that give the object number for each DEEP2 object.
  if keyword_set(pcat) then begin
      if cnt eq 1 then $
        phdex = where(pcat.objno eq $
                      long(objectcat[objdex].object), cnt) $
      else begin
          num = cnt
          cnt = 0
          phdex = intarr(num)
          for i=0,num-1 do begin
              phdex[i] = where(pcat.objno eq $
                               long(objectcat[objdex[i]].object), any)
              cnt = cnt + any
          endfor
      endelse
; if object's slitlet is not found in pcat, then return an
; error message.
      if cnt eq 0 then begin
          print, '(find_objpos.pro) Object(s) ' + $
            strcompress(objectcat[objdex].object, /rem) + $
            ' in slit(s) ' + string(slitn, format='(I3.3)') + $
            ' not found in PhotCat table!'
; return the values, but return a fwhm = 0 and objpos = 0.
          set_nulls, fwhm=fwhm
          return, objpos
      endif
; define the cfht pixel scale (arcseconds per pixel).
      cfht_pixscl = 0.207
; extract the CFHT seeing value from the header of the pcat and
; convert the value from CFHT pixels to arcseconds so that it can be
; compared to the gaussian radii of the objects.
      seeing = sxpar(phdr, 'SEEING')
      seeing = float(seeing) * cfht_pixscl * 2.35
; translate the gaussian radius (in arcseconds) to a fwhm in DEIMOS
; pixels. major axis = 1.2 * rg in units of arcseconds.
      maskpa = mbasics.pa_pnt
      ellip = pcat[phdex].e2
      pa = pcat[phdex].pa
      rg = pcat[phdex].rg       ;or rg = objectcat[desidex].majaxis / 1.2
; calculate the fwhm value.
      fwhm = 2.35 / pixscl * rg * cfht_pixscl * $
        sqrt( cos((pa - maskpa)*!dtor)^2 + $
              ( (1.- ellip) * $ 
                sin((pa - maskpa)*!dtor) )^2 )
; convert the fwhm value from a double to a float.
      fwhm = float(fwhm)
; make sure that fwhm is greater than or equal to the seeing.
      fwhm = fwhm > seeing
; else return an array of null values.
  endif else fwhm = fltarr(n_elements(objpos))
; return the position of the object (in units of DEIMOS pixels).
  return, objpos

  
end





