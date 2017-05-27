
;+
;
; NAME
;      find_objpos_mod.pro
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

function find_objpos_mod, slitn, nrows, fwhm=fwhm, bfile=bfile,$
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
;  print,'bintab file : ',bfile
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
;nzs  deimos_isdeep, isdeep
  isdeep=0                      ; nzs
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

;  print,objpos
;  print,mapdex
;  print,slitobjmap[mapdex].topdist 
;  print,slitobjmap[mapdex].botdist

 for mapidex=0,cnt-1 do begin
  if slitobjmap[mapdex[mapidex]].topdist eq slitobjmap[mapdex[mapidex]].botdist then begin
     ;;objpos=0.
     print,'*** error in position, OK if Alignment Star *** : slit : ',slitn
  endif 
 endfor

  
;print,objpos,mapdex,slitobjmap[mapdex].topdist,slitobjmap[mapdex].botdist,pixscl,desilen,nrows,edgeloss
;if slitn le 94 then begin
; help,/struct,slitobjmap[mapdex]
;endif
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
     rg = pcat[phdex].rg        ;or rg = objectcat[desidex].majaxis / 1.2
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
  return, objpos(0)
  
  
end



pro get_object_positions,maskname,outputfile

  bintab = maskname + '.bintabs.fits'

  print,'processing mask :',maskname
  color = ['B','R']
  
  srctable = mrdfits(bintab,1,/silent)

  jslit=where(strtrim(srctable.objclass,2) eq 'Program_Target' or strtrim(srctable.objclass,2) eq 'Alignment_Star',ct)
  nslit=ct
  names=make_array(nslit,/string) 
  slit=make_array(2,nslit,/int)
  num=make_array(2,nslit,/int)
  predicted=make_array(2,nslit,/float)
  offsets=make_array(2,nslit,/float,value=1.e4)
  obj=make_array(2,nslit,/int)
  slit_widths=make_array(2,nslit,/float,value=0.)
  
  
  nobj_max=500
  meas_center=make_array(nobj_max,/float)
  meas_disp=make_array(nobj_max,/float) 
  meas_chisq=make_array(nobj_max,/float) 
  meas_slit=make_array(nobj_max,/int)   
  meas_color=make_array(nobj_max,/int)   
  


  tol=20.
  iobj=0
;nslit=15
  FOR islit=0,nslit-1 DO BEGIN
     FOR ic=0, 1 DO BEGIN
        j_best=-1
        slitno = jslit[islit]   ;
        good=1                  ; assume good, flag if bad
        
        if slitno LT 1000 THEN slitname = '.' + strtrim(slitno,2)
        if slitno LT 100 THEN slitname = '.0' + strtrim(slitno,2)
        if slitno LT 10 THEN slitname = '.00' + strtrim(slitno,2)
        
        filename = 'slit.' + maskname + slitname + color(ic) + '.fits'
        names[islit]='slit.' + maskname + slitname
print,filename 
        goodfile= file_test(filename)
        if goodfile eq 0 then begin
           filename = 'spSlit.' + maskname + slitname + color(ic) + '.fits'
           goodfile= file_test(filename)
        endif

        if goodfile eq 0 then print,' **** could not find file : ',filename
        if goodfile GT 0 THEN BEGIN
           profivar=1
           sprof=find_object(filename,profivar=profivar, $
                 /cr,/bpm)
           plot,sprof
;
; in searching for a source, exclude end 2 pix in wavelength
;
           exclude=2
           ns=n_elements(sprof)-2*exclude
           chan=exclude+indgen(ns)
           sprof_short=sprof[chan]
;
; get slit dimensions and predicted position
;
           filename = 'slit.COS' + slitname + color(ic) + '.fits'
           
           goodfile= file_test(filename)
           if goodfile eq 0 then begin
              filename = 'slit.' + maskname + slitname + color(ic) + '.fits'
              goodfile= file_test(filename)
           endif
           
           
           if goodfile eq 0 then print,' **** could not find file : ',filename
            if goodfile GT 0 THEN BEGIN

               slit = mrdfits(filename,1,/silent)
               slitsize=slit.rawsize
               dithersize=slit.dithersize/2.0
               slit_widths[ic,islit]=1.*n_elements(sprof)
               predicted(ic,islit)=find_objpos_mod(slitno,slitsize)+dithersize
               goodfit=1
               while goodfit do begin

                  fitprofile=gaussfit(indgen(ns),sprof_short, fitparam,sigma=sigma,nterms=6,chisq=chisq)
                  
                  resid = sprof_short-fitprofile
                  oplot,chan,fitprofile    
                  
                  if fitprofile[0] eq -1 or fitparam[0] lt 0. or $
                     fitparam[1] lt 3 or fitparam[1] gt ns-3 or $
                     fitparam[2] lt 4. or fitparam[2] gt (ns/2) or $
                     fitparam[0]/sigma[0] lt 4. then begin
                     goodfit=0
                  endif else begin
                     print,' found an object -- offset,dithersize & fitted parameters : ', $
                           fitparam[1]-predicted(ic,islit),dithersize
                     print,fitparam
                     print,sigma
                     sprof_short=sprof_short-fitprofile
                     meas_center[iobj]=fitparam[1]
                     meas_disp[iobj]=fitparam[2]
                     meas_chisq[iobj]=chisq
                     meas_slit[iobj]=islit
                     meas_color[iobj]=ic
                     num[ic,islit]=num[ic,islit]+1
                     offset=(meas_center[iobj]-predicted[ic,islit])
                     if abs(offset) le abs(offsets[ic,islit]) then begin
                        obj[ic,islit]=iobj
                        offsets[ic,islit]=offset
                     endif
                     iobj=iobj+1
                  endelse
                  
;  read,ig
               endwhile
            endif
         endif
        
        
     endfor
  endfor
  meas_center=meas_center[0:iobj-1]
  meas_disp=meas_disp[0:iobj-1]
  meas_chisq=meas_chisq[0:iobj-1]
  meas_slit=meas_slit[0:iobj-1]
  meas_color=meas_color[0:iobj-1]
;
; offsets contains the min offset to each predicted source 
;for each slit and color
;
  j=where(abs(offsets) lt 100.,ct)
  print,'offsets :',offsets[j]
  med_blue=median(offsets[j])
  med_red=med_blue
  fwhm=2.35*avg(meas_disp[j])
  print,'mean FWHM (pixels) : ',fwhm
  j=where(abs(offsets[0,*]-med_blue) lt tol,ct)
  if ct gt 2 then begin
     print,'BLUE low offsets : ',reform(offsets[0,j])
     med=median(offsets[0,j])
     moments=moment(offsets[0,j],sdev=sdev)
     print,'BLUE : ',ct,' objects -- mean, median, sdev : ',moments[0],med,sdev
     print,'removing outliers'
     j=where(abs(offsets[0,*]-med) lt tol/4.,ct)
     print,'BLUE low offsets : ',reform(offsets[0,j])
     med=median(offsets[0,j])
     moments=moment(offsets[0,j],sdev=sdev)
     print,'BLUE : ',ct,' objects -- mean, median, sdev : ',moments[0],med,sdev
     med_blue=med
  endif
  
  j=where(abs(offsets[1,*]-med_red) lt tol,ct)
  if ct gt 2 then begin
     print,'RED low offsets : ',reform(offsets[1,j])
     med=median(offsets[1,j])
     moments=moment(offsets[1,j],sdev=sdev)
     print,'RED : ',ct,' objects -- mean, median, sdev : ',moments[0],med,sdev
     print,'removing outliers'
     j=where(abs(offsets[1,*]-med) lt tol/4.,ct)
     print,'RED low offsets : ',reform(offsets[1,j])
     med=median(offsets[1,j])
     moments=moment(offsets[1,j],sdev=sdev)
     print,'RED : ',ct,' objects -- mean, median, sdev : ',moments[0],med,sdev
     med_red=med
  endif
  med_ave=(med_blue+med_red)/2.
  print,'blue , red and average offsets : ',med_blue,med_red,med_ave
  use_ave=1
  if abs(med_blue-med_red) gt 5. then use_ave=0
  blue_offset=med_blue
  red_offset=med_red
  if use_ave then begin
     blue_offset=med_ave
     red_offset=med_ave
  endif
;
; generate output list including program sources + any 2ndary sources
;
  openw,unit,outputfile,/get_lun
  centers=predicted
  for i=0,nslit-1 do begin
     centers[0,i]=centers[0,i]+blue_offset
     centers[1,i]=centers[1,i]+red_offset
;
; print out a line for each program source if within range
;
     if (centers[0,i]) gt 2. and (centers[1,i]) gt 2. and $
        (centers[0,i]) lt (slit_widths[0,i]-2.) and $
        (centers[1,i]) lt (slit_widths[1,i]-2.) then begin
        
        line=names[i]+' '+string(centers[0,i],format='(f7.1)')+ $
             ' '+string(centers[1,i],format='(f7.1)')+ $
             ' '+string(fwhm,format='(f7.1)')+ $
             ' '+string(fwhm,format='(f7.1)')
        printf,unit,line
     endif
;
; now treat cases w/ a possible secondary
; require the secondary be displaced by at least sec_tol pixels
; from program source
;
     sec_tol=7.
     j=where(meas_slit eq i,ct)
     if ct gt 0 then begin
        jb=where(meas_color[j] eq 0 and abs(meas_center[j]-centers[0,i]) gt sec_tol,ctb)
        if ctb gt 0 then jsb=j[jb]
        jr=where(meas_color[j] eq 1 and abs(meas_center[j]-centers[1,i]) gt sec_tol,ctr)
        if ctr gt 0 then jsr=j[jr]
;
;   now make a list of matched secondaries
;   a lot of awful code since we want to make the red and blue sec
;   correspond but not have duplicates
;
        if ctb ne 0 and ctr ne 0 then begin
           meas_2b=0.0 & meas_2r=0.0 & fwhm_2b=0. & fwhm_2r=0.
           if ctb gt 0 then begin
              meas_2b=[meas_2b,meas_center[jsb]]
              meas_2r=[meas_2r,meas_center[jsb]+red_offset-blue_offset]
              fwhm_2b=[fwhm_2b,2.35*meas_disp[jsb]]
              fwhm_2r=[fwhm_2r,2.35*meas_disp[jsb]]
           endif
           if ctr gt 0 then begin
              meas_2r=[meas_2r,meas_center[jsr]]
              meas_2b=[meas_2b,meas_center[jsr]+blue_offset-red_offset]
              fwhm_2b=[fwhm_2b,2.35*meas_disp[jsr]]
              fwhm_2r=[fwhm_2r,2.35*meas_disp[jsr]]
           endif
           meas_2b=meas_2b[1:*]
           meas_2r=meas_2r[1:*]
           fwhm_2b=fwhm_2b[1:*]
           fwhm_2r=fwhm_2r[1:*]
           if ctb gt 0 and ctr gt 0 then begin
              ind=indgen(n_elements(meas_2b))
              for ii=0,n_elements(meas_2b)-1 do begin
                 jsame=where(abs(meas_2b[ii]-meas_2b) lt sec_tol $
                             and ind ne ii and meas_2b ne -1000.,ct)
                 if ct gt 0 then begin
                    meas_2b[jsame]=-1000.
                    meas_2r[jsame]=-1000.
                 endif
              endfor
              ind=indgen(n_elements(meas_2r))
              for ii=0,n_elements(meas_2r)-1 do begin
                 jsame=where(abs(meas_2r[ii]-meas_2r) lt sec_tol $
                             and ind ne ii and meas_2r ne -1000.,ct)
                 if ct gt 0 then begin
                    meas_2b[jsame]=-1000.
                    meas_2r[jsame]=-1000.
                 endif
              endfor
              j=where(meas_2b ne -1000. and meas_2r ne -1000.,ct)
              if ct gt 0 then begin
                 meas_2b=meas_2b[j]
                 meas_2r=meas_2r[j]
                 fwhm_2b=fwhm_2b[j]
                 fwhm_2r=fwhm_2r[j]
              endif
           endif
           nsec=n_elements(meas_2r)
           for ii=0,nsec-1 do begin
;
; print out a line for each secondary source if it's w/i range
;
              if (meas_2b[ii]) gt 2. and (meas_2r[ii]) gt 2. and $
                 (meas_2b[ii]) lt (slit_widths[0,i]-2.) and $
                 (meas_2r[ii]) lt (slit_widths[1,i]-2.) then begin
                 
                 line=names[i]+' '+string(meas_2b[ii], $
                                          format='(f7.1)')+' '+string(meas_2r[ii],format='(f7.1)')+ $
                      ' '+string(fwhm_2b[ii],format='(f7.1)')+ $
                      ' '+string(fwhm_2r[ii],format='(f7.1)')
                 printf,unit,line
                 
              endif
           endfor
        endif
     endif
  endfor
  close,unit
  free_lun,unit
  print,'***** outputs in file : ',outputfile
  
end
