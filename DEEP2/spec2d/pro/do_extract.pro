;+
;
; NAME
;      do_extract.pro
;
; PURPOSE
;      Makes iterative calls to extract1d.pro to extract 1d spectra
;      from slit files. The extracted spectra are packaged according
;      to slit number (blue and red portions together) and saved to a
;      .fits file.
;
;
; SYNTAX
;      do_extract, [files=files, _EXTRA=EXTRA]
;
; INPUTS
;      files = an optional parameter giving a vector of slit file
;              names. The procedure calls the extract1d.pro and
;              extracts a 1d spectrum from each of the slit files. 
;
; KEYWORDS
;      
;
; OUTPUTS
;      A .fits file which contains the extracted 1d spectra for a
;      given slit on a given mask. The files are named
;      "spec1d.xxxx.nnn.fits" where xxxx denotes the mask number and
;      nnn gives the slit number. Within the fits file, the blue
;      portion of the spectra is saved as the first extension and the
;      red portion of the spectrum is saved as the second
;      extension. Each is a structure containing flux, wavelength, and
;      inverse variance information. For example:
;           bfoo = mrdfits('spec1d.xxxx.nnn.fits', 1)
;           rfoo = mrdfits('spec1d.xxxx.nnn.fits', 2)
;           bfoo.spec = the flux as a function of wavelength on the
;                       blue end.
;           bfoo.lambda = the wavelength solution for the blue end (in
;                         linear lambda form). 
;           bfoo.ivar = the inverse variance as a function of
;                       wavelength on the blue end.
;
; PROCEDURES CALLED 
;      extract1d.pro
;      mcc_gauss1dfit.pro
;      find_object.pro
;      peakinfo.pro
;      find_objpos.pro
;      fxaddpar
;      headfits
;      mrdfits
;      mwrfits
;
; EXAMPLES
;      None.
;
; COMMENTS
;      The do_extract.pro routine assumes that you are already in the
;      directory which contains the slit files. This is a reasonable
;      assumption since this routine is intended to be called within
;      the domask.pro wrapper which makes the same assumption. Note
;      that the PBS shell script cd's to the proper directory. 
;
;      If the user passes do_extract 100 or more slitfiles, then the
;      routine will determine the seeing difference between the
;      photmetric observations (design fwhm) and the spectroscopy
;      observations (fwhm from spatial profile). This difference will
;      be applied to the design fwhm's and used in the extraction. The
;      position of the object will be chosen according to the significance
;      of the peak in the spatial profile. If a strong peak, then
;      pkcol (from peakinfo) will be used. If not so strong, the
;      design position will be employed.
;      If you feed do_extract less than 100 slitfiles, then it will
;      use the fwhm according to the design specifications (without
;      adjusting for seeing changes between the photometric
;      observations and the present spectroscopy data). The position
;      will be chosen in the same manner according to the significance
;      of the peak in the spatial profile.
;
; HISTORY
;      Created July 22, 2002 by mcc.
;      Revised July 28, 2002 by mcc - routine adjusted to accomodate     
;         changes in the routines find_object, peakinfo, and
;         extract1d. 
;      Revised July 29, 2002 by mcc - routine revised to handle the
;         varying lengths of structures. The structures read-in from
;         the slitfiles are now stored using pointers.
;      Revised August 4, 2002 by mcc - reversed revision from July
;        29th. Saving structures using pointers was too inefficient
;        and so now the routine simply reads the files multiple times.
;      Revised August 26, 2002 by mcc - added /CREATE keyword to call
;         to MWRFITS so that the spec1d.xxxx.nnn.fits files are
;         re-written each time that the pipeline is run. Also added
;         code to account for multiple objects in a single slit..
;      Hacked 21sep02 by MD to simplify object finding.
;      Revised October 21, 2002 by mcc - complete revision of routine!
;      Revised November 4, 2002 by mcc - another complete revision!
;      Revised March 1, 2003 by mcc - routine now extracts spectra
;         according to the /boxsprof and /horne extractions. Tests
;         with a fake data showed that these extractions gave better
;         results in the cases of bad pixel columns. The /boxsprof was
;         far superior to the /boxcar and the /horne was very
;         minimally better than the /optimal.
;-

; a simple procedure to find only the unique solutions in a array. 
pro select_uniq, bstr, rstr, bfin=bfin, rfin=rfin, isdeep=isdeep
  flagval = 0
  sz = n_elements(bstr)
  for i=0,sz-1 do begin
      if bstr[i].objpos eq flagval then flag = 1 else flag = 0 
      dex = where(bstr.objpos eq bstr[i].objpos, dcnt)
      if i eq 0 then subscr = dex[0] $
      else begin
          junk = where(subscr eq dex[0], jcnt)
          if jcnt eq 0 or flag then subscr = [subscr, dex[0]]
      endelse
  endfor
  bfin = bstr[subscr]
  rfin = rstr[subscr]
  num = n_elements(subscr)
  if keyword_set(isdeep) then begin
      alpha = ['b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k']
      lchar = alpha[indgen(num)]
      objno = strmid(bfin[0].objno, 0, 9) + lchar
      bfin.objno = objno
      rfin.objno = objno
  endif else begin
      objnames = 'serendip' + string((indgen(num) + 1), format='(I0.3)')
      bfin.objno = objnames
      rfin.objno = objnames
  endelse
end
;--------------------------------
; a function to determine the name of a serendip object.
function get_objname, inum, objpos, ss, pcat, ra=ra, dec=dec
; initialize the variables ra and dec to be null-strings.
  ra = ''
  dec = ''
; construct a way to convert numbers to letters...this is bad!
  alpha = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', $
           'i', 'j', 'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r']
; determine the number of catalog objects in the slit.
  catnum = n_elements(ss)
  for i=0,catnum-1 do begin
; determine the distance between the serendip and the known object.
      if ss[i].objpos gt 0. then pixdist = abs(objpos - ss[i].objpos) $
      else pixdist = abs(objpos - ss[i].cat_objpos)
; get the RA and DEC for the catalog object.
      dex = where(pcat.objno eq long(ss[i].objno), cnt)
      if cnt gt 0 then begin
          dec0 = pcat[dex[0]].dec
          ra0 = pcat[dex[0]].ra
; determine the RA and DEC offsets to the serendip.
          slit_xy2radec, pixdist, dec0, pixel_scale=pixscl, $
            slitpa=slitpa, maskpa=maskpa, $
            delta_ra=delta_ra, delta_dec=delta_dec, /degrees
; calculate the RA and DEC of the serendip.
          if keyword_set(dec) then dec = [dec, dec0 + delta_dec] $
            else dec = [dec0 + delta_dec]
          if keyword_set(ra) then ra = [ra, ra0 + delta_ra]
          ra = [ra0 + delta_ra]
      endif
  endfor
; if there are multiple catalog objects then take the average RA and
; DEC values.
  if keyword_set(dec) then dec = mean(dec) 
  if keyword_set(ra) then ra = mean(ra)
; now convert the RA and DEC from decimal degrees to HR:MIN:SEC
; format.
  if ra ne '' and dec ne '' then begin
      radec, ra, dec, hr, mn, sc, deg, min, sec
      ra = strn(hr,  padtype=1, padchar='0', format='(I2.2)') + ':' + $
        strn(mn,  padtype=1, padchar='0', format='(I2.2)') + ':' + $
        strn(sc, padtype=1, padchar='0', format='(F5.2)', length=5)
      if deg lt 0 then sign = '-' else sign = '+'
      dec = sign + $
        strn(deg,  padtype=1, padchar='0', format='(I2.2)') + ':' + $
        strn(min, padtype=1, padchar='0', format='(I2.2)') + ':' + $
        strn(sec, padtype=1, padchar='0', format='(F4.1)', length=4)
  endif else begin
      ra = ''
      dec = ''
  endelse
; lastly construct the name for the serendip object.
  name = 's' + strcompress(string(ss[0].objno), /rem) + alpha[inum]
  return, name
end
;--------------------------------
; a procedure to extract the RA and DEC of an object from the objectcat.
function get_radec, objno, objcat=objcat
  if keyword_set(objcat) then begin
      dex = where(strcompress(objcat.object, /rem) eq objno, cnt)
      if cnt gt 0 then begin
          dec = float(objcat[dex[0]].dec_obj)
          ra = float(objcat[dex[0]].ra_obj)
          radec, ra, dec, hr, mn, sc, deg, min, sec
          ra = strn(hr,  padtype=1, padchar='0', format='(I2.2)') + ':' + $
            strn(mn,  padtype=1, padchar='0', format='(I2.2)') + ':' + $
            strn(sc, padtype=1, padchar='0', format='(F5.2)')
          if deg lt 0 then sign = '-' else sign = '+'
          dec = sign + $
            strn(deg,  padtype=1, padchar='0', format='(I2.2)') + ':' + $
            strn(min, padtype=1, padchar='0', format='(I2.2)') + ':' + $
            strn(sec, padtype=1, padchar='0', format='(F4.1)')
      endif else begin
          dec = ''
          ra = ''
      endelse
  endif
;  if keyword_set(pcat) then begin
;      dex = where(pcat.objno eq long(objno), cnt)
;      if cnt gt 0 then begin
;          dec = pcat[dex[0]].dec
;          ra = pcat[dex[0]].ra
;      endif else begin
;          dec = ''
;          ra = ''
;      endelse
;  endif
  return, [ra, dec]
end

;--------------------------------
; a procedure to extract the magnitudes (B,R,I) of an object from the
; pcat.
function get_mags, objno, pcat=pcat, objpa=objpa
  dd = where(strcompress(string(pcat.objno), /rem) eq $
             strcompress(objno, /rem), cnt)
  if cnt eq 0 then begin
      objpa = 999.0
      output = [0.0, 0.0, 0.0]
  endif else begin
      output = [pcat[dd[0]].magb, $
                pcat[dd[0]].magr, $
                pcat[dd[0]].magi]
      objpa = pcat[dd[0]].pa
  endelse
  return, output
end


pro do_extract, files=files, nonlocal=nonlocal, $
                nsigma_optimal=nsigma_optimal, $
                nsigma_boxcar=nsigma_boxcar, _extra=extra

; if no slit files are passed to the routine, then simply find all
; slit files in the current directory.
  if n_elements(files) eq 0 then files = findfile('slit*.*')
; if no slit files are passed or found, then return error message.
  nfiles = n_elements(files)
  if nfiles eq 0 then message, 'ERROR: no slitfiles supplied by ' + $
    'user or found in current directory!'

; if the nonlocal keyword is set, then extract from the extension
; in the slitfile that contains the non-local-sky-subtracted data.
  if keyword_set(nonlocal) then goto, jump_nonlocal

;--------------------------------
; step 1: match the pairs of slit files (red and blue pairs) and
; determine objects positions, fwhms, and serendips. tabulate all the
; information in an array of structures.

; define the scale in pixels at which two peaks in the spatial
; profile will be considered to be denoting the same object.
  resolu = 5.
; define the scale over which to smooth the spatial profile when
; searching for serendips.
  smthscl = 3.
; define the threshold level for object detection.
  thresh = 15.
  thresh_serendip = 10.
; define the window and nbuf values passed to the peakinfo routine.
  window = 5.
  nbuf = 2.
; define the template structure to catch the plethora of info
; regarding each object.
  template = {objno:'', slitno:long(0), slitfile:'', $
              objtype:'', color:'', cat_objpos:float(0), $
              cat_fwhm:float(0), objpos:float(0), fwhm:float(0), $
              corr_fwhm:float(0), nrows:float(0), $
              s2n_fwhm:float(0), s2n_window:float(0), $
              ra:'', dec:'', Xmm:float(-1), Ymm:float(-1), $
              magb:999.0, magr:999.0, magi:999.0, objpa:999.0}
; make the inital done array.
  done_arr = [-1]

; loop through the slitfiles...
  for i=0,nfiles-1 do begin
; check if this particular slitfile has already been processed with
; its counterpart.
      ddex = where(done_arr eq i, dcnt)
      if dcnt eq 0 then begin
; now for the ith file in the list of slit files, find all files
; (blue/red) that correspond to the same slit number. that is, find
; the pair (blue and red) files that contain all of the data for the
; slit.
          loc = strpos(files[i], '.fits')
          slitdex = where( strmid(files, 0, loc-1) eq $
                           strmid(files[i], 0, loc-1) and $
                           files ne files[i], slitcnt )
; for each of the slitfiles (at most 2 - red & blue), get the object
; position, fwhm, etc.
          for j=0,slitcnt do begin
; get the right file.
              if j eq 0 then this_file = files[i] $
              else this_file = files[slitdex[j-1]]
; extract the data from the file.
              slit = mrdfits(this_file, 1, hdr, /silent) 
; get the mask number from the file header. only do this once for the
; set of slitfiles (we can assume that all files are from the same
; mask).
              if i eq 0 then begin
                  mask = strcompress( sxpar(hdr, 'SLMSKNAM'), /rem)
; and if the mask is a DEEP2 mask, then remove the ".E" or ".W" at the
; end of the mask name (only keep the xxxx numbers). to do this we
; must check the bin table file to see if the mask is a DEEP2 mask!
                  bintab = findfile('*.bintabs.fits*', count=bincnt)
                  if bincnt eq 0 then $
                    message, '(do_extract.pro) bintab file not found!' $
                  else begin
                      bintab = bintab[0]
; open bin table file and get the extension names...remember to close
; the opened file!!!
                      fits_open, bintab, fcb
                      extnames = fcb.extname
                      fits_close, fcb
; determine if this is a DEEP2 mask...find the mask design table in
; the bintab file and from it extract the PROJECT name.
                      deimos_isdeep, isdeep, mask
; if this is a DEEP2 mask, then get the pcat info.
                      if isdeep then begin
                          tabdex = where(extnames eq 'PhotCat', cnt)
                          if cnt eq 0 then $
                            print, '(do_extract.pro) PhotCat ' + $
                            'table not found in bintab file!' $
                          else begin
; from the pcat table (idl structure), read in the pcat and get the
; seeing value from the header....seeing is in CFHT pixels.
                              pcat = mrdfits(bintab, tabdex[0], phdr, /silent)
                              psee = sxpar(phdr, 'SEEING')
                          endelse
                      endif
; also extract the ObjectCat table from the bintab file.
                      tabdex = where(extnames eq 'ObjectCat', cnt)
                      if cnt eq 0 then $
                        print, '(do_extract.pro) ObjectCat table ' + $
                        'not found in bintab file!' $
                      else objcat = mrdfits(bintab, tabdex[0], /silent)
; extract the Xmm and Ymm positions from the bintab file. The Xmm and
; Ymm values give the position of the slit (at its center) on the
; DEIMOS detector in units of milimeters.                      
                      bludex = where(extnames eq 'BluSlits', blucnt)
                      desdex = where(extnames eq 'DesiSlits', descnt)
                      if blucnt eq 0 or descnt eq 0 then begin
                          if blucnt eq 0 then $
                            print, '(do_extract.pro) BluSlits table ' + $
                            'not found in bintab file!'
                          if descnt eq 0 then $
                            print, '(do_extract.pro) DesiSlits table ' + $
                            'not found in bintab file!'
                      endif else begin
                          blutab = mrdfits(bintab, bludex[0], /silent)
                          desitab = mrdfits(bintab, desdex[0], /silent)
; sort the BluSlits and DesiSlits tables by dslitid number.
                          blutab = blutab[sort(blutab.dslitid)]
                          desitab = desitab[sort(desitab.dslitid)]
                          nslit = n_elements(blutab)
                          slitinfo = {slitn:lonarr(nslit), $
                                      Xmm:fltarr(nslit), $
                                      Ymm:fltarr(nslit)}
                          slitinfo.Xmm = (blutab.slitx1 + blutab.slitx2 + $
                                          blutab.slitx3 + blutab.slitx4) / 4.0
                          slitinfo.Ymm = (blutab.slity1 + blutab.slity2 + $
                                          blutab.slity3 + blutab.slity4) / 4.0
                          slitinfo.slitn = long(desitab.slitname)
                      endelse
                  endelse
              endif ;if i eq 0 (only for first slitfile in list)
; determine the number of pixels in the slitfile.
              npix_rows = n_elements(slit.flux[0,*])
; get the slit number from the file header and extract the object
; position, fwhm, object number, object type from the mask design
; tables.
              slitnum = long( sxpar(hdr, 'SLITNO') )
              desipos = find_objpos(slitnum, npix_rows, fwhm=cat_width, $
                                    objnum=objnum, objtype=objkind)
              objnum = strcompress( string(objnum), /rem)
; call peakinfo to determine position of object based on the spatial
; profile and to measure the fwhm of the profile.
              sprof = find_object(slit, profivar=profivar, npix=npix)
              peakinfo, sprof, pk_pixel, width, pk_quad=pkq, pk_cent=pkc, $
                profivar=profivar, npix=npix, s2n_window=s2n_win, $
                nbuf=nbuf, window=window, s2n_fwhm=s2n_width, signif=thresh
; if the peak position determined via centroiding is good, then use
; it. otherwise use the peak pixel position. recall that the peakinfo
; routine returns the value -1 for pk_cent if the centroid fails.
              pkcol = float(pk_pixel)
              gdex = where(pkc ge 0., ngood)
              if ngood gt 0 then pkcol[gdex] = pkc[gdex]
; determine if this a blue or red portion of the slit.
              loc = strpos(this_file, '.fits')
              hue = strmid(this_file, loc-1, 1)
; determine how many objects are expected to be in this slit.
              catnum = n_elements(desipos)
; create an array of structures to hold the info about these objects.
              str = replicate(template, catnum)
              for k=0,catnum-1 do begin
                  str[k].slitfile = this_file
                  str[k].slitno = slitnum
                  str[k].objno = objnum[k]
                  str[k].objtype = objkind
                  str[k].cat_objpos = desipos[k]
                  str[k].cat_fwhm = cat_width[k]
                  str[k].nrows = npix_rows
                  str[k].color = hue
                  dd = where(slitinfo.slitn eq slitnum, ddcnt)
                  if ddcnt gt 0 then begin
                      str[k].Xmm = slitinfo.Xmm[dd[0]]
                      str[k].Ymm = slitinfo.Ymm[dd[0]]
                  endif
                  if keyword_set(objcat) then begin
                      radec = get_radec(objnum[k], objcat=objcat)
                      str[k].ra = radec[0]
                      str[k].dec = radec[1]
                  endif else begin
                      str[k].ra = ''
                      str[k].dec = ''
                  endelse
                  if isdeep and keyword_set(pcat) then begin
                      deep_mags = get_mags(objnum[k], pcat=pcat, objpa=objpa)
                      str[k].magb = deep_mags[0]
                      str[k].magr = deep_mags[1]
                      str[k].magi = deep_mags[2]
                      str[k].objpa = objpa
                  endif
; now compare the locations of the peaks found in the spatial profile
; to the object position given by the design specifications.
                  mindiff = min(abs(desipos[k] - pkcol), minpos)
                  print,desipos[k]
                  if mindiff le resolu then begin
                      str[k].objpos = pkcol[minpos]
                      str[k].fwhm = width[minpos]
                      str[k].s2n_fwhm = s2n_width[minpos]
                      str[k].s2n_window = s2n_win[minpos]
                  endif
; if this is DEEP2 data, then check if the slit is a sky-only slit.
                  if isdeep then begin
                      if str[k].objtype eq 'S' then begin
; if it is a sky-only slit, then extract at the middle of the slit.
                          str[k].objpos = npix_rows / 2.
                          str[k].fwhm = width[0]
                          str[k].s2n_fwhm = s2n_width[0]
                          str[k].s2n_window = s2n_win[0]
                      endif
                  endif
              endfor ;k-index
; now go back and look for serendips. to do this, let's smooth the
; spatial profile a bit.
              sprof = ivarsmooth(sprof, profivar, smthscl, profivar)
; and call peakinfo again to find possible serendip detections.
              peakinfo, sprof, pkpix, fwhm, npix=npix, pk_cent=pkc, $
                pk_quad=pkq, profivar=profivar, s2n_window=s2n_win, $
                s2n_fwhm=s2n_fwhm, nbuf=nbuf, window=window, $
                signif=thresh_serendip, /serendip

; exclude all peaks found near the original objects. also exclude all
; peaks for which the fwhm is bogus.
              pkcol = float(pkpix)
              gdex = where(pkc gt 0., gdcnt)
              if gdcnt gt 0 then pkcol[gdex] = pkc[gdex]
              pknum = n_elements(pkcol)
              str_ser = replicate(template, pknum)
              qnum = 1
              for q=0,pknum-1 do begin
                  mindiff = min( abs(str.cat_objpos - pkcol[q]), mindex )
; only take the serendip if it is at least 5 pixels away from any
; program object position and has a fwhm with a reasonable value.
                  if mindiff gt resolu and $
                    fwhm[q] gt 2 and fwhm[q] lt 20 then begin
                      if isdeep then $
                        str_ser[q].objno = $
                        get_objname(qnum, pkcol[q], str, pcat, $
                                    ra=ra, dec=dec) $
                      else begin
                          str_ser[q].objno = 'serendip' + $
                            string(qnum, format='(I0.3)')
                          ra = '0.0'
                          dec = '0.0'
                      endelse
                      str_ser[q].slitfile = this_file
                      str_ser[q].objtype = 'Q'
                      str_ser[q].slitno = slitnum
                      str_ser[q].color = hue
                      dd = where(slitinfo.slitn eq slitnum, ddcnt)
                      if ddcnt gt 0 then begin
                          str_ser[q].Xmm = slitinfo.Xmm[dd[0]]
                          str_ser[q].Ymm = slitinfo.Ymm[dd[0]]
                      endif
                      str_ser[q].objpos = pkcol[q]
                      str_ser[q].fwhm = fwhm[q]
                      str_ser[q].s2n_fwhm = s2n_fwhm[q]
                      str_ser[q].s2n_window = s2n_win[q]
                      str_ser[q].nrows = npix_rows
                      str_ser[q].ra = ra
                      str_ser[q].dec = dec
                      qnum = qnum + 1
                  endif
              endfor ;q-index
; trim off unfilled structures. or rather only keep the filled
; structures. 
              fill = where(str_ser.objno ne '', fillcnt)
              if fillcnt gt 0 then $
                str = [str, str_ser[fill]]

; add the information for this slit to the composite structure.
              if i eq 0 and j eq 0 then finstr = str $
              else finstr = [finstr, str]
; lastly before iterating, add these entries in the slit files array
; to the array containing all the entries which have been analyzed.
              if i eq 0 then begin
                  if slitcnt gt 0 then done_arr = [i, slitdex] $
                  else done_arr = [i] 
              endif else begin
                  if slitcnt gt 0 then done_arr = [done_arr, i, slitdex] $
                  else done_arr = [done_arr, i]
              endelse
          endfor ;j-index
; compare the serendips detections on the red and blue side of the
; slit. first, step back and find where in the output structures the
; entries for the serendips are for this slit.
          bdex = where(finstr.slitno eq slitnum and $
                       finstr.color eq 'B' and $
                       finstr.objtype eq 'Q', bcnt)
          rdex = where(finstr.slitno eq slitnum and $
                       finstr.color eq 'R' and $
                       finstr.objtype eq 'Q', rcnt)
          totcnt = bcnt + rcnt
          if bcnt gt 0 and rcnt gt 0 then begin
              for p=0,bcnt-1 do begin
                  mindiff = min(abs(finstr[rdex].objpos - $
                                    finstr[bdex[p]].objpos), mindex)
                  if mindiff lt resolu then begin
                      if p eq 0 then rtmp = (finstr[rdex])[mindex] $
                      else rtmp = [rtmp, (finstr[rdex])[mindex]]
                  endif else begin 
                      if p eq 0 then rtmp = finstr[bdex[p]] $
                        else rtmp = [rtmp, finstr[bdex[p]]]
                  endelse
              endfor
              for p=0,rcnt-1 do begin
                  mindiff = min(abs(finstr[bdex].objpos - $
                                    finstr[rdex[p]].objpos), mindex)
                  if mindiff lt resolu then begin
                      if p eq 0 then btmp = (finstr[bdex])[mindex] $
                      else btmp = [btmp, (finstr[bdex])[mindex]]
                  endif else begin
                      if p eq 0 then btmp = finstr[rdex[p]] $
                      else btmp = [btmp, finstr[rdex[p]]]
                  endelse
              endfor
; remember to reset the color tags!
              bnum = n_elements(btmp)
              btmp.color = strarr(bnum) + 'B'
              btmp.slitfile = strarr(bnum) + finstr[bdex[0]].slitfile
              rnum = n_elements(rtmp)
              rtmp.color = strarr(rnum) + 'R'
              rtmp.slitfile = strarr(rnum) + finstr[rdex[0]].slitfile
              bstr = [finstr[bdex], btmp]
              rstr = [rtmp, finstr[rdex]]
; now sort these arrays to keep only the unique elements.
              if isdeep then $
                select_uniq, bstr, rstr, bfin=bfin, rfin=rfin, /isdeep $
              else select_uniq, bstr, rstr, bfin=bfin, rfin=rfin
          endif
; examine cases where serendips are only found on one side of slit.
          if bcnt gt 0 and rcnt eq 0 and slitcnt gt 0 then begin
              rfin = finstr[bdex]
              rfin.color = strarr(bcnt) + 'R'   
              rfind = where(finstr.slitno eq slitnum and $
                            finstr.color eq 'R', rfindcnt)
              rfin.slitfile = finstr[rfind[0]].slitfile
              bfin = finstr[bdex]
          endif
          if rcnt gt 0 and bcnt eq 0 and slitcnt gt 0 then begin
              bfin = finstr[rdex]
              bfin.color = strarr(rcnt) + 'B'
              bfind = where(finstr.slitno eq slitnum and $
                            finstr.color eq 'B', bfindcnt)
              bfin.slitfile = finstr[bfind[0]].slitfile
              rfin = finstr[rdex]
          endif
; combine this with the final structure array.
          if bcnt gt 0 or rcnt gt 0 and slitcnt gt 0 then begin
              othdex = where(finstr.slitno ne slitnum or $
                             finstr.objtype ne 'Q', othcnt)
              if othcnt gt 0 and slitcnt gt 0 then $
                finstr = [finstr[othdex], bfin, rfin]
              if othcnt gt 0 and slitcnt eq 0 then begin
                  if rcnt eq 0 then finstr=[finstr[othdex], bfin]
                  if bcnt eq 0 then finstr=[finstr[othdex], rfin]
              endif
              if othcnt eq 0 then message, '(do_extract.pro) ERROR: ' + $
                'error in indexing structure!'
          endif
      endif
  endfor ;i-index
; make sure that the cat_fwhm is greater than or equal to the fwhm as
; given by the pcat seeing.
  if isdeep then begin
      pixscl = 0.117371 ;arcseconds/pixel
      cfht_pixscl = 0.207 ;arcseconds/pixel
      finstr.cat_fwhm = (finstr.cat_fwhm) > $
        (2.35 * psee * cfht_pixscl / pixscl)
      print, 'Seeing fwhm value from pcat: ' + $
        string((2.35 * psee * cfht_pixscl / pixscl), format='(F5.2)') + $
        ' (DEIMOS pixels).'
  endif

; TBD:not done by this stage...comparison of the positions on the red and
; blue side for the known objects. 


;--------------------------------
; step 2: do a comparison of the fwhm as measured from the spatial
; profile and the fwhm from the photometry catalog (aka
; pcat). calculate the seeing difference between the spectroscopic
; observations and the photometry observations and apply it to the pcat
; fwhm values. we will then use these seeing-corrected fwhm values for
; the extraction widths. do this comparison separately for the blue
; and red sides of the slits since the seeing differs with
; wavelength. recall that the seeing at keck is typically better in
; the blue! also in this step of the routine, write the finstr
; structure to a file (obj_info .fits file).

; separate the red and blue sides of slits.
  bdex = where(finstr.color eq 'B', bcnt)
  rdex = where(finstr.color eq 'R', rcnt)
; check where the fwhm (as determined from the spatial profile) was
; unmeasured by the peakinfo routine.
  flagdex = where(finstr.fwhm le 0., flagcnt)
  print, '(do_extract.pro): ' + strcompress(string(flagcnt), /rem) + $
    ' objects had spatial profiles for which a FWHM was unmeasured.'
; check for cases where either fwhm or cat_fwhm are bogus. we only
; want to compare the seeing difference where we have valid data from
; both the design specifications and from the spatial profile.
  if bcnt gt 0 then begin
      dex = where(finstr[bdex].fwhm gt 0. and $
                  finstr[bdex].cat_fwhm gt 0., cnt)
; check that there are enough valid slits over which to do the seeing
; comparison. if so, compare fwhm to cat_fwhm for the blue side.
      if cnt gt 75 then begin
          print, '(do_extract.pro): Doing seeing comparison on blue side...'
; calculate seeing difference for only the valid points.
          seediff = (finstr[bdex].fwhm^2 - finstr[bdex].cat_fwhm^2)
          seediff = seediff[dex]
; set the bin size for the histogram plot of the seeing difference.
          binsz = 5
; bin the seeing diff data.
          plothist, seediff, xhist, yhist, bin=binsz, /noplot
; fit a gaussian to the binned data.
          mcc_gauss1dfit, xhist, yhist, a=a, yfit=yfit
; check that the fit parameters are reasonable before applying them.
          if finite(a[1]) then corr = a[1] else corr = 0.
; also determine seeing difference by simply taking mean and median of
; the distribution.
          djs_iterstat, seediff, mean=avg_seediff, median=med_seediff, $
            sigma=unc_seediff, sigrej=3.0
          print, 'Seeing correction factor (blue) from Gaussian ' + $
            'fitting to histogram = ' + string(corr, format='(F6.2)')
          print, 'Mean seeing correction factor (blue) after outlier ' + $
            'rejection = ' + string(avg_seediff, format='(F6.2)')
          print, 'Median seeing correction factor (blue) after outlier ' + $
            'rejection = ' + string(med_seediff, format='(F6.2)')
; define an array of the various correction factors.
          corrB = [corr, avg_seediff, med_seediff]
; use the median seeing difference value.
          corr = med_seediff
; define a variable to hold the seeing difference between pcat and
; spectroscopic data for the blue side (in units of arcseconds).
          pixscl = 0.117371 ;deimos pixel scale in arcsec/pixel
          if corr lt 0 then bsd = -sqrt(abs(corr)) * pixscl $
          else bsd = sqrt(corr) * pixscl
; catch possible instances where cat_fwhm^2 + corr < 0 such that we
; don't get -NaNs when we take the sqrt. in these cases, just use the
; pcat fwhm value.
          if corr lt 0 then begin
              badpts = where( (finstr[bdex].cat_fwhm^2 + corr) lt 0., badwh)
              fwhm_tmp = finstr[bdex].cat_fwhm
              blu_fwhm = sqrt(finstr[bdex].cat_fwhm^2 + corr)
              if badwh gt 0 then blu_fwhm[badpts] = fwhm_tmp[badpts]
          endif else blu_fwhm = sqrt(finstr[bdex].cat_fwhm^2 + corr)
          finstr[bdex].corr_fwhm = blu_fwhm
; make a post-script plot of the seeing difference distribution.
          set_plot, 'ps'
          device, file='ps/seediff_blue.ps', /landscape, /color   
          lilx = (min(xhist) - abs(0.1 * min(xhist))) > (-750)
          bigx = (max(xhist) + abs(0.1 * max(xhist))) < 750
          plothist, seediff, bin=binsz, thick=2, xthick=1.5, ythick=1.5, $
            charsize=1.5, xtitle='(fwhm!U2!N!X - cat_fwhm!U2!N!X)', $
            ytitle='Number (Total = ' + string(cnt, format='(1I0)') + ')', $
          title=mask + ' (blue side)', /xstyle, xr=[lilx,bigx]
          oplot, xhist, yfit, thick=2, line=2, color=2
          legend, ['Amplitude = ' + string(a[0], format='(F6.2)'), $
                   'x0 = ' + string(a[1], format='(F6.2)'), $
                   'sigma = ' + string(a[2], format='(F6.2)'), $
                   'mean value = ' + string(avg_seediff, format='(F6.2)'), $
                   'median value = ' + string(med_seediff, format='(F6.2)')], $
            /right_legend
          device, /close
          set_plot, 'x' 
      endif else begin
          print, '(do_extract.pro) No pcat info or not enough ' + $
            'slitfiles to do seeing comparison (on blue side)!'
          finstr[bdex].corr_fwhm = finstr[bdex].fwhm
      endelse
  endif
  
; okay, now do the same analysis for the red side of each slit.
  if rcnt gt 0 then begin
      dex = where(finstr[rdex].fwhm gt 0. and $
                  finstr[rdex].cat_fwhm gt 0., cnt)
; check that there are enough valid slits over which to do the seeing
; comparison. if so, compare fwhm to cat_fwhm for the blue side.
      if cnt gt 75 then begin
          print, '(do_extract.pro): Doing seeing comparison on red side...'
; calculate seeing difference for only the valid points.
          seediff = (finstr[rdex].fwhm^2 - finstr[rdex].cat_fwhm^2)
          seediff = seediff[dex]
; set the bin size for the histogram plot of the seeing difference.
          binsz = 5
; bin the seeing diff data.
          plothist, seediff, xhist, yhist, bin=binsz, /noplot
; fit a gaussian to the binned data.
          mcc_gauss1dfit, xhist, yhist, a=a, yfit=yfit
; check that the fit parameters are reasonable before applying them.
          if finite(a[1]) then corr = a[1] else corr = 0.
; also determine seeing difference by simply taking mean and median of
; the distribution.
          djs_iterstat, seediff, mean=avg_seediff, median=med_seediff, $
            sigma=unc_seediff, sigrej=3.0
          print, 'Seeing correction factor (red) from Gaussian ' + $
            'fitting to histogram = ' + string(corr, format='(F6.2)')
          print, 'Mean seeing correction factor (red) after outlier ' + $
            'rejection = ' + string(avg_seediff, format='(F6.2)')
          print, 'Median seeing correction factor (red) after outlier ' + $
            'rejection = ' + string(med_seediff, format='(F6.2)')
; define an array of the various correction factors.
          corrR = [corr, avg_seediff, med_seediff]
; use the median seeing difference value.
          corr = med_seediff
; define a variable to hold the seeing difference between pcat and
; spectroscopic data for the red side (in units of arcseconds).
          pixscl = 0.117371 ;deimos pixel scale in arcsec/pixel
          if corr lt 0 then rsd = -sqrt(abs(corr)) * pixscl $
          else rsd = sqrt(corr) * pixscl
; catch possible instances where cat_fwhm^2 + corr < 0 such that we
; don't get -NaNs when we take the sqrt. in these cases, just use the
; pcat fwhm value.
          if corr lt 0 then begin
              badpts = where( (finstr[rdex].cat_fwhm^2 + corr) lt 0., badwh)
              fwhm_tmp = finstr[rdex].cat_fwhm
              red_fwhm = sqrt(finstr[rdex].cat_fwhm^2 + corr)
              if badwh gt 0 then red_fwhm[badpts] = fwhm_tmp[badpts]
          endif else red_fwhm = sqrt(finstr[rdex].cat_fwhm^2 + corr)
          finstr[rdex].corr_fwhm = red_fwhm
; make a post-script plot of the seeing difference distribution.
          set_plot, 'ps'
          device, file='ps/seediff_red.ps', /landscape, /color   
          lilx = (min(xhist) - abs(0.1 * min(xhist))) > (-750)
          bigx = (max(xhist) + abs(0.1 * max(xhist))) < 750
          plothist, seediff, bin=binsz, thick=2, xthick=1.5, ythick=1.5, $
            charsize=1.5, xtitle='(fwhm!U2!N!X - cat_fwhm!U2!N!X)', $
            ytitle='Number (Total = ' + string(cnt, format='(1I0)') + ')', $
            title=mask + ' (red side)', /xstyle, xr=[lilx,bigx]
          oplot, xhist, yfit, thick=2, line=2, color=2
          legend, ['Amplitude = ' + string(a[0], format='(F6.2)'), $
                   'x0 = ' + string(a[1], format='(F6.2)'), $
                   'sigma = ' + string(a[2], format='(F6.2)'), $
                   'mean value = ' + string(avg_seediff, format='(F6.2)'), $
                   'median value = ' + string(med_seediff, format='(F6.2)')], $
            /right_legend
          device, /close
          set_plot, 'x'
      endif else begin
          print, '(do_extract.pro) No pcat info or not enough ' + $
            'slitfiles to do seeing comparison (on red side)!'
          finstr[rdex].corr_fwhm = finstr[rdex].fwhm
      endelse
  endif

; while we are making plots, let's make a post-script plot of the
; objpos - cat_objpos distribution for the set of objects for which
; both were measured.
; set bin size for histogram.
  binsz = 1.0
  if bcnt gt 0 then begin
      hist_dat1 = (finstr.objpos - finstr.cat_objpos)[bdex]
      dex = where(finstr[bdex].objpos gt 0. and $
                  finstr[bdex].cat_objpos gt 0., cnt1)
      if cnt1 gt 0 then hist_dat1 = hist_dat1[dex]
      if cnt1 gt 50 then $
        plothist, hist_dat1, xhist1, yhist1, bin=binsz, /noplot $
      else print, 'Not enough objects to create objpos ' + $
        'histogram on blue side!'
  endif else cnt1 = 0
  if rcnt gt 0 then begin
      hist_dat2 = (finstr.objpos - finstr.cat_objpos)[rdex] 
      dex = where(finstr[rdex].objpos gt 0. and $
                  finstr[rdex].cat_objpos gt 0., cnt2)
      if cnt2 gt 0 then hist_dat2 = hist_dat2[dex]
      if cnt2 gt 50 then $
        plothist, hist_dat2, xhist2, yhist2, bin=binsz, /noplot $
      else print, 'Not enough objects to create objpos ' + $
        'histogram on red side!'
  endif else cnt2 = 0
  nslits = (cnt1 + cnt2) / 2
; make plot.
  if cnt1 gt 50 and cnt2 gt 50 then begin
      set_plot, 'ps'
      device, file='ps/objpos_hist.ps', /landscape
      plothist, hist_dat1, bin=binsz, thick=2, $
        xtitle='observed object position - design object position (pixels)', $
        ytitle='Number (total = ' + string(nslits, format='(I3.0)') + $
        ' slits)', /xsty, /ysty, $
        xr=[min([xhist1,xhist2])-binsz, max([xhist1,xhist2])+binsz], $
        yr=[min([yhist1,yhist2]), max([yhist1,yhist2])]
      plothist, hist_dat2, bin=binsz, /overplot, linesty=2, thick=2
      legend, ['Solid Line = Blue Slits', 'Dashed Line = Red Slits']
      device, /close
      set_plot, 'x'
  endif

; finally, save and analyze the info that was tabulated in step 1.
; create a file name in which to wrote the object info data.
  objfile = 'obj_info.' + mask + '.fits'
; copy the header from one of the slit files. remove the fields
; referring directly to the particular slit.
  hdr = copy_header(hdr, 'ObjInfo')
  sxdelpar, hdr, 'SLITNO'
  sxdelpar, hdr, 'SLITX0'
  sxdelpar, hdr, 'SLITX1'
  if keyword_set(corrB) then $
    sxaddpar, hdr, 'corrB', corrB[2], 'median seeing correction on Blue side'
  if keyword_set(corrR) then $
    sxaddpar, hdr, 'corrR', corrR[2], 'median seeing correction on Red side'
; add header entries for pcat seeing and for the seeing difference
; (both blue and red). 
  if keyword_set(psee) then $
    sxaddpar, hdr, 'pcat_see', psee, 'pcat seeing (cfht pixels)'
  if keyword_set(bsd) then $
    sxaddpar, hdr, 'SeeDiffB', bsd, 'Seeing Diff (arcsec)'
  if keyword_set(rsd) then $
    sxaddpar, hdr, 'SeeDiffR', rsd, 'Seeing Diff (arcsec)'

; write the object information structure into a fits file.
  mwrfits, finstr, objfile, hdr, /silent, /create


; if we are extracting from the non-local-sky-subtracted extension,
; then read back in the obj_info data.
  jump_nonlocal: if keyword_set(nonlocal) then begin
      objfile = findfile('*obj_info*.fits', count=nobjfiles)
      if nobjfiles eq 0 then $
        message, '(do_extract.pro) No obj_info file found!'
      finstr = mrdfits(objfile[0], 1, objhdr, /silent)
      mask = strcompress(sxpar(objhdr, 'SLMSKNAM'), /rem)
      deimos_isdeep, isdeep, mask
; make sure we define the keylist.
      keylist = ['FRAMENO', 'OUTFILE', 'EXPTIME', 'DARKTIME', 'OBSERVER', $
                 'OBJECT', 'OBSTYPE', 'ROTATVAL', 'DATE-OBS', 'UT', 'AIRMASS', $
                 'TARGNAME', 'EPOCH', 'EQUINOX', 'DEC',  'RA', 'AZ', 'EL', $    
                 'HA', 'ST', 'MJD-OBS', 'PARANG',  'SYNOPSIS', 'DWFILNAM', $
                 'SLMSKNAM', 'GRATEPOS', 'HPLOGTIM',  $
                 'SP2DVERS', 'AUTHOR', 'CHIPNO', 'GRATING', $
                 'GRTLTWAV', 'SLITNO', 'SLITXO', 'SLITX1', 'SLITWID', $
                 'SKYSIGMA','WAVETYPE']
  endif
;--------------------------------
; step3: iterate through the list of slit files (as stored in the finstr
; structure). for each slit, extract the 1-d spectra (blue and red
; portions), then store the 1-d spectra as structures in a new spec1d
; .fits file. 

; initiliaze the objdone array. we will redefine it later; presently,
; it just needs to be define for the first pass thru the for/do loop.
  objdone = [-1]
; determine the number of objects which we need to extract.
  nfiles = n_elements(finstr)
  for j=0,nfiles-1 do begin
; get all entries matching the jth object number and slit
; number. really, object number should be enough here!
      objdex = where(finstr.objno eq finstr[j].objno and $
                     finstr.slitno eq finstr[j].slitno, objcnt)
; make sure that this object wasn't extracted already.
      done = where(objdone eq j, done_num)
; if not done previously and 2 entries in finstr are found, then do
; the following...
      if done_num eq 0 and objcnt gt 0 then begin
          if objcnt gt 2 then $
            print, '(do_extract.pro) ERROR: more than 2 objects found ' + $
            'with same object number!'
; figure out which entry is the blue portion of the slit and which is
; the red portion.
          bdex = where(finstr[objdex].color eq 'B', bcnt)
          if bcnt gt 1 then $
            print, '(do_extract.pro) ERROR: multiple objects ' + $
            '(w/ same objno) found in blue portion of slit ' + $
            strcompress(string(finstr[j].slitno), /rem)
          if bcnt GT 0 then bdex = objdex[bdex[0]]
          rdex = where(finstr[objdex].color eq 'R', rcnt)
          if rcnt gt 1 then $
            print, '(do_extract.pro) ERROR: multiple objects ' + $
            '(w/ same objno) found in red portion of slit ' + $
            strcompress(string(finstr[j].slitno), /remove_all)
          if rcnt gt 0 then rdex = objdex[rdex[0]]
      endif else begin
          bcnt = 0
          rcnt = 0
      endelse

; get position(s) of object.
;TBD: perhaps add a keyword to do_extract so that this s/n value can
;be set by the user....or at least so that the user can set it to a
;very small value to force the measured objpos and fwhm to be used in
;the extractions.
      sng = 1.3
      if bcnt gt 0 then begin
; if the objpos value is non-zero and the signal-to-noise is decent,
; then use it. otherwise, use the cat_objpos value. 
          if finstr[bdex].objpos gt 0 and $
            finstr[bdex].s2n_fwhm gt sng then posB = finstr[bdex].objpos $
          else begin
;              print, '(do_extract.pro) Positions for object ' + $
;                finstr[bdex].objno + ' in slit ' + $
;                string(finstr[bdex].slitno, format='(I3.3)') + $
;                ' disagree by more than 5 pixels on blue side...' + $
;                'using cat_objpos!'
              if finstr[bdex].objtype ne 'Q' then $
                posB = finstr[bdex].cat_objpos $
              else posB = finstr[bdex].objpos
          endelse
      endif
      if rcnt gt 0 then begin
          if finstr[rdex].objpos gt 0 and $
            finstr[rdex].s2n_fwhm gt sng then posR = finstr[rdex].objpos $
          else begin
;              print, '(do_extract.pro) Positions for object ' + $
;                finstr[rdex].objno + ' in slit ' + $
;                string(finstr[rdex].slitno, format='(I3.3)') + $
;                ' disagree by more than 5 pixels on red side...' + $
;                'using cat_objpos!'
              if finstr[rdex].objtype ne 'Q' then $
                posR = finstr[rdex].cat_objpos $ 
              else posR = finstr[rdex].objpos
          endelse
      endif
      
; get the extraction widths. if the signal-to-noise is poor, then use
; the catalog width.
      if bcnt gt 0 then begin
          if (finstr[bdex].s2n_fwhm gt sng) or $
             (finstr[bdex].cat_fwhm eq 0.0) or $
             (finstr[bdex].objtype eq 'Q') then $
            fwhmB = finstr[bdex].fwhm $
          else fwhmB = finstr[bdex].cat_fwhm
      endif
      if rcnt gt 0 then begin
          if (finstr[rdex].s2n_fwhm gt sng) or $
             (finstr[rdex].cat_fwhm eq 0.0) or $
             (finstr[rdex].objtype eq 'Q') then $
            fwhmR = finstr[rdex].fwhm $
          else fwhmR = finstr[rdex].cat_fwhm
      endif
; match the extraction widths for the tophat extraction so that the
; continuum levels for the blue and red side of the spectrum will be
; equal (only for the tophat extraction).
      if rcnt gt 0 and bcnt gt 0 then begin
          avgfwhmB = (fwhmB + fwhmR)/2.
          avgfwhmR = avgfwhmB
      endif else begin
          if bcnt gt 0 then avgfwhmB = fwhmB
          if rcnt gt 0 then avgfwhmR = fwhmR
      endelse
      
; get the slit files.
      if bcnt gt 0 then bfile = finstr[bdex].slitfile
      if rcnt gt 0 then rfile = finstr[rdex].slitfile

; specify the extraction width parameters (nsigma).
      if n_elements(nsigma_optimal) gt 0 then $
        nsig_opt = nsigma_optimal[0] else nsig_opt = 1.5
      if n_elements(nsigma_boxcar) gt 0 then $
        nsig_box = nsigma_boxcar[0] else nsig_box = 1.1 ;/ 2.35482

; extract the blue spectrum via the optimal extraction and via the
; tophat (or boxcar) extraction algorithm.
      if bcnt gt 0 then begin
          if keyword_set(nonlocal) then begin
              blu_opt = extract1d(bfile, posB, fwhmB, /horne, $
                                  /nonlocal, nsigma=nsig_opt)
              blu_box = extract1d(bfile, posB, avgfwhmB, $
                                  /nonlocal, nsigma=nsig_box)
          endif else begin
              blu_opt = extract1d(bfile, posB, fwhmB, /horne, $
                                  nsigma=nsig_opt)
              blu_box = extract1d(bfile, posB, avgfwhmB, $
                                  /boxsprof, nsigma=nsig_box)
          endelse
      endif
; similarly, extract the red spectrum.
      if rcnt gt 0 then begin
          if keyword_set(nonlocal) then begin
              red_opt = extract1d(rfile, posR, fwhmR, /horne, $
                                  /nonlocal, nsigma=nsig_opt)
              red_box = extract1d(rfile, posR, avgfwhmR, /boxsprof, $
                                  /nonlocal, nsigma=nsig_box)
          endif else begin
              red_opt = extract1d(rfile, posR, fwhmR, /horne, $
                                  nsigma=nsig_opt)
              red_box = extract1d(rfile, posR, avgfwhmR, $
                                  /boxsprof, nsigma=nsig_box)
          endelse
      endif
      
; define the name of the 1-d output file name.
      specfile = 'spec1d.' + mask + '.' + $
        string(finstr[j].slitno, format='(I3.3)') + '.' + $
        strcompress(finstr[j].objno, /rem) + '.fits'

; write the 1-d spectra (blue and red / optimal and tophat) to a
; single fits file. first, grab the header from one of the slitfiles
; and carry it over to the 1-d file.
      hdr = headfits(finstr[j].slitfile, ext=1, /silent)
; write the blue portion of the tophat/boxcar extraction.
      if bcnt gt 0 then begin
          if keyword_set(nonlocal) then $
;            hdrB = copy_header(hdr, 'Boxcar-NL') $
            hdrB = copy_header(hdr, 'Bxspf-NL-B') $
;          else hdrB = copy_header(hdr, 'Boxcar')
          else hdrB = copy_header(hdr, 'Bxspf-B')
          sxaddpar, hdrB, 'objno', finstr[bdex].objno, $
            'Object number', after='DATE'
          sxaddpar, hdrB, 'objpos', finstr[bdex].objpos, $
            'Object Position from Spatial Profile', after='objno'
          sxaddpar, hdrB, 'fwhm', finstr[bdex].fwhm, $
            'FWHM from Spatial Profile', after='objpos'
          sxaddpar, hdrB, 'cat_objpos', finstr[bdex].cat_objpos, $
            'Design Object Position', after='fwhm'
          sxaddpar, hdrB, 'cat_fwhm', finstr[bdex].cat_fwhm, $
            'FWHM from PCAT', after='cat_objp'
          sxaddpar, hdrB, 'cor_fwhm', finstr[bdex].corr_fwhm, $
            'FWHM from PCAT corrected for seeing diff', after='cat_fwhm'
          sxaddpar, hdrB, 'ext_fwhm', avgfwhmB, $
            'FWHM employed in object extraction', after='cor_fwhm'
; add RA and DEC info is object is a serendip.
          if finstr[bdex].objtype eq 'Q' then begin
              sxaddpar, hdrB, 'RA_obj', finstr[bdex].ra, $
                'estimate of serendip RA', after='ext_fwhm'
              sxaddpar, hdrB, 'DEC_obj', finstr[bdex].dec, $
                'estimate of serendip DEC', after='RA_obj'
          endif else begin
              sxaddpar, hdrB, 'RA_obj', finstr[bdex].ra, $
                'RA of object', after='ext_fwhm'
              sxaddpar, hdrB, 'DEC_obj', finstr[bdex].dec, $
                'DEC of object', after='RA_obj'
          endelse  
          sxaddpar, hdrB, 'Xmm', finstr[bdex].xmm, $
            'X position of slit', after='DEC_obj'
          sxaddpar, hdrB, 'Ymm', finstr[bdex].ymm, $
            'Y position of slit', after='Xmm'
          sxaddpar, hdrB, 'magB', finstr[bdex].magb, $
            'B magnitude', after='Ymm'
          sxaddpar, hdrB, 'magR', finstr[bdex].magr, $
            'R magnitude', after='magB'
          sxaddpar, hdrB, 'magI', finstr[bdex].magi, $
            'I magnitude', after='magR'
          sxaddpar, hdrB, 'OBJPA', finstr[bdex].objpa, $
            'Object PA on sky', before='SLITPA'

          if bcnt gt 0 and size(blu_box, /tname) eq 'STRUCT' then begin
              print, '(do_extract.pro) Writing spec1d file: ' + $
                specfile + ' ......'
              if keyword_set(nonlocal) then $
                mwrfits, blu_box, specfile, hdrB, /silent $
              else mwrfits, blu_box, specfile, hdrB, /silent, /create
          endif
      endif
; write the red portion of tophat/boxcar extraction.
      if rcnt gt 0 then begin
          if keyword_set(nonlocal) then $
;            hdrR = copy_header(hdr, 'Boxcar-NL') $
            hdrR = copy_header(hdr, 'Bxspf-NL-R') $
;          else hdrR = copy_header(hdr, 'Boxcar') 
          else hdrR = copy_header(hdr, 'Bxspf-R')
          sxaddpar, hdrR, 'objno', finstr[rdex].objno, $
            'Object number', after='DATE'
          sxaddpar, hdrR, 'objpos', finstr[rdex].objpos, $
            'Object Position from Spatial Profile', after='objno'
          sxaddpar, hdrR, 'fwhm', finstr[rdex].fwhm, $
            'FWHM from Spatial Profile', after='objpos'
          sxaddpar, hdrR, 'cat_objpos', finstr[rdex].cat_objpos, $
            'Design Object Position', after='fwhm'
          sxaddpar, hdrR, 'cat_fwhm', finstr[rdex].cat_fwhm, $
            'FWHM from PCAT', after='cat_objp'
          sxaddpar, hdrR, 'cor_fwhm', finstr[rdex].corr_fwhm, $
            'FWHM from PCAT corrected for seeing diff', after='cat_fwhm'
          sxaddpar, hdrR, 'ext_fwhm', avgfwhmR, $
            'FWHM employed in object extraction', after='cor_fwhm'
; add RA and DEC info is object is a serendip.
          if finstr[rdex].objtype eq 'Q' then begin
              sxaddpar, hdrR, 'RA_obj', finstr[rdex].ra, $
                'estimate of serendip RA', after='ext_fwhm'
              sxaddpar, hdrR, 'DEC_obj', finstr[rdex].dec, $
                'estimate of serendip DEC', after='RA_obj'
          endif else begin
              sxaddpar, hdrR, 'RA_obj', finstr[rdex].ra, $
                'RA of object', after='ext_fwhm'
              sxaddpar, hdrR, 'DEC_obj', finstr[rdex].dec, $
                'DEC of object', after='RA_obj'
          endelse
          sxaddpar, hdrR, 'Xmm', finstr[rdex].xmm, $
            'X position of slit', after='DEC_obj'
          sxaddpar, hdrR, 'Ymm', finstr[rdex].ymm, $
            'Y position of slit', after='Xmm'
          sxaddpar, hdrR, 'magB', finstr[rdex].magb, $
            'B magnitude', after='Ymm'
          sxaddpar, hdrR, 'magR', finstr[rdex].magr, $
            'R magnitude', after='magB'
          sxaddpar, hdrR, 'magI', finstr[rdex].magi, $
            'I magnitude', after='magR'
          sxaddpar, hdrR, 'OBJPA', finstr[rdex].objpa, $
            'Object PA on sky', before='SLITPA'
          if size(red_box, /tname) eq 'STRUCT' then begin
              if bcnt gt 0 then $
                mwrfits, red_box, specfile, hdrR, /silent $
              else begin 
                  print, '(do_extract.pro) Writing spec1d file: ' + $
                    specfile + ' ......'
                  if keyword_set(nonlocal) then $
                    mwrfits, red_box, specfile, hdrR, /silent $
                  else mwrfits, red_box, specfile, hdrR, /silent, /create
              endelse
          endif
      endif
; write the blue portion of the optimal extraction.
      if bcnt gt 0 then begin
          if keyword_set(nonlocal) then $
;            sxaddpar, hdrB, 'EXTNAME', 'Optimal-NL', 'Extension Name' $
            sxaddpar, hdrB, 'EXTNAME', 'Horne-NL-B', 'Extension Name' $
;          else sxaddpar, hdrB, 'EXTNAME', 'Optimal', 'Extension Name'
          else sxaddpar, hdrB, 'EXTNAME', 'Horne-B', 'Extension Name'
          sxaddpar, hdrB, 'ext_fwhm', fwhmB, $
            'FWHM employed in object extraction'
          if size(blu_opt, /tname) eq 'STRUCT' then $
            mwrfits, blu_opt, specfile, hdrB, /silent
      endif
; write the red spectrum of optimal extraction.
      if rcnt gt 0 then begin
          if keyword_set(nonlocal) then $
;            sxaddpar, hdrR, 'EXTNAME', 'Optimal-NL', 'Extension Name' $
            sxaddpar, hdrR, 'EXTNAME', 'Horne-NL-R', 'Extension Name' $
;          else sxaddpar, hdrR, 'EXTNAME', 'Optimal', 'Extension Name' 
          else sxaddpar, hdrR, 'EXTNAME', 'Horne-R', 'Extension Name'
          sxaddpar, hdrR, 'ext_fwhm', fwhmR, $
            'FWHM employed in object extraction'
          if size(red_opt, /tname) eq 'STRUCT' then $
            mwrfits, red_opt, specfile, hdrR, /silent
      endif
      
; construct/add to the list of object for which extraction has been
; completed.
      if j eq 0 then begin
          if bcnt gt 0 and rcnt gt 0 then $
            objdone = [bdex, rdex] $
          else begin
              if bcnt gt 0 then objdone = [bdex] 
              if rcnt gt 0 then objdone = [rdex]
          endelse
      endif else begin
          if bcnt gt 0 then objdone = [objdone, bdex]
          if rcnt gt 0 then objdone = [objdone, rdex]
      endelse
      
  endfor


  
end








