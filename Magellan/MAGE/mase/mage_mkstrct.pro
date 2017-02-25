pro mage_mkstrct, mage, rawpath=rawpath

  cmd = "\ls "+strtrim(rawpath,2)+"mage*.fit*"

  spawn, cmd, files
  ;; remove the mage.fits file which causes problems
  all_files = fileandpath(files)
  ifiles = WHERE(strmatch(all_files, 'mage.fits') EQ 0, nmage)
  If nmage GT 0 THEN files = files[ifiles] $
  ELSE message, 'Problem with directory, files not found'

  nfiles = n_elements(files)

   tmp = {use:     1,  $
         rawpath: ' ', $
         frame:    0L, $
         fitsfile: ' ', $
         exptype:  ' ', $
         object:   ' ', $
         airmass:  -1.0, $
         slit:     ' ', $
         exptime:  0.0, $
         jd:       0d, $
         ra_deg:   0.0, $
         dec_deg:  0.0, $
         obj_id:   -1, $
         arcfile:  ' ', $
         pixflatfile: ' ', $
         illumflatfile: ' ', $
         orderfile: ' ', $
         slitfile: ' ', $
         stdfile:  ' ', $
         sensfunc: ' ' }
          ;;hdr:      strarr(500), $

;  tmp = {magestrct, $
;         use:     1,  $
;         rawpath: rawpath, $
;         fitsfile: '', $
;         arcfile:  '', $
;         exptype:  '', $
;         object:   '', $
;         airmass:  -1.0, $
;         slit:     '', $
;         exptime:  0.0, $
;         jd:       0d, $
;         ra_deg:   0.0, $
;         dec_deg:  0.0, $
;         obj_id:   -1 $
;         }

  ;;tmp = { magestrct }

   mage = replicate(tmp, nfiles)

  for ifile=0, nfiles-1 do begin

     hdr = xheadfits(strtrim(files[ifile],2))

     mage[ifile].rawpath   = rawpath
     ;;mage[ifile].hdr       = hdr
     mage[ifile].fitsfile  = sxpar(hdr, "FILENAME")+'.fits'
     temp = strsplit(mage[ifile].fitsfile, 'mage*.fits', /extract)
     mage[ifile].frame    = long(temp[0])
     mage[ifile].exptime   = sxpar(hdr, "EXPTIME")
     mage[ifile].object    = sxpar(hdr, "OBJECT")
     mage[ifile].ra_deg    = sxpar(hdr, "RA-D")
     mage[ifile].dec_deg   = sxpar(hdr, "DEC-D")
     mage[ifile].slit      = sxpar(hdr, "SLITNAME")
     mage[ifile].airmass      = sxpar(hdr, "AIRMASS")

;    Exposure types may be the following: 
;    TWIPIX   = Twilight pixel gain correction flat 
;    TWITRC   = Twilight trace/illum flat with the science slit
;    DOMEFLT  = Dome flat (used for both pixel and trace in the red
;    XE-FLASH = Xe flash lamp flat
;    STD      = standard star observation
;    SCIENCE  = science frame
;    ARC      = ThAr arc lamp

     type     = sxpar(hdr, "EXPTYPE")
     airmass  = sxpar(hdr, "AIRMASS")
     xe_flash = sxpar(hdr, "XE-FLASH")
     thar     = sxpar(hdr, "HC_LAMP")
     slit     = sxpar(hdr, "SLITNAME")
     object   = sxpar(hdr, "OBJECT")
     exptime  = sxpar(hdr, "EXPTIME")

     ; Tag the structure with the reduce JD to match Arcs 
     ; to the closest time science frames (taken at same RA/DEC)

     date_day = float(strsplit(sxpar(hdr,"UT-DATE"),"-", /extract))
     date_time = float(strsplit(sxpar(hdr,"UT-TIME"),":", /extract))
     date = [date_day[0],date_day[1],date_day[2],date_time[0],$
             date_time[1],date_time[2]]
     juldate, date, jd          
     mage[ifile].jd = jd

     ; Set the exposure type using crude set of rules.

     if (type EQ 'Flat') THEN $ ;;AND airmass EQ '1.0') then $
        mage[ifile].exptype = "DOMEFLT"

     if (slit NE '5.00' AND airmass NE '1.00'  $
         AND type NE 'ThAr-Lamp' AND type NE 'Flat' AND type NE 'Xe-Flash') $
     then mage[ifile].exptype = "TWITRC"

     if (type EQ "ThAr-Lamp") then $
        mage[ifile].exptype = "ARC"

     if (type EQ "Xe-Flash") then $
        mage[ifile].exptype = "XE-FLASH"

     if (slit EQ '5.00' AND airmass NE '1.00') then $
        mage[ifile].exptype = "TWIPIX"


     if (type EQ "Object") then begin
        if (strpos(strlowcase(object), "gd") EQ -1 AND (strpos(strlowcase(object), "GD") EQ -1 AND $
            strpos(strlowcase(object), "feige") EQ -1) AND strpos(strlowcase(object), "Feige") EQ -1) then begin
           mage[ifile].exptype = "SCIENCE"
        endif else begin
           mage[ifile].exptype = "STD"
        endelse
     endif

     if (type EQ "Object") AND exptime LE 450 then begin
        mage[ifile].exptype = "BRIGHT"
        
     endif


     
  endfor


; Try to assign an object ID number to each unique object observed.

  objframes = where(mage.exptype EQ 'SCIENCE' OR mage.exptype EQ 'STD'  OR mage.exptype EQ 'BRIGHT')
  objectnames = mage[objframes[uniq(mage[objframes].object)]].object

  for ifile=0, nfiles-1 do begin

     if (mage[ifile].exptype NE "ARC") then begin

        match = where(mage.exptype EQ "ARC" AND $ ;Look for arcs
                      $
                      $  ; Look for either same sky coords or zenith (for domes)
                      ((abs(mage.ra_deg-mage[ifile].ra_deg) LT 0.1 AND $
                        abs(mage.dec_deg-mage[ifile].dec_deg) LT 0.1) OR $
                       (float(mage[ifile].airmass) EQ 1.0)) $ 
                       $
                       AND $
                       $
                       $        ; Match Slits
                       (mage.slit EQ mage[ifile].slit OR $
                        mage[ifile].exptype EQ "TWIPIX"), narcs)

        if (mage[ifile].exptype EQ "SCIENCE" AND narcs EQ 0) then begin
           match = where(mage.slit EQ mage[ifile].slit AND mage.exptype EQ "ARC", narcs)
           print, "WARNING: Arc and science frames widely spaced: ", mage[ifile].fitsfile
        endif
        
        if (narcs EQ 0) then begin
           print, "No Arcs found for this setup: ", mage[ifile].fitsfile

           continue
        endif else begin
           arcs = mage[match]
        endelse

        ;; Closest spaced arc in time.
        bestarc = arcs[where(abs(arcs.jd-mage[ifile].jd) EQ min(abs(arcs.jd-mage[ifile].jd)))]
        ;;????? I would add something here that computes the angle on the
        ;; sky also so that we don't use arcs from very different elevations
        

        mage[ifile].arcfile = bestarc[0].fitsfile

     endif


;  Done matching arcs, now look to match an object ID to this frame
     
     if (mage[ifile].exptype EQ 'SCIENCE' OR $
         mage[ifile].exptype EQ 'STD' OR $
         mage[ifile].exptype EQ 'BRIGHT' ) then begin


        for iobj=0, n_elements(objectnames)-1 do begin
           if (mage[ifile].object EQ objectnames[iobj]) then begin
              mage[ifile].obj_id = iobj+1   ;adding an offset here so that the first object is picked up by mage_pipe

              break
           endif
        endfor
     endif

  endfor


  ; Print out a summary.
  for ifile=0, nfiles-1 do begin
     if (mage[ifile].obj_id NE -1) then begin
        print, mage[ifile].fitsfile, mage[ifile].object, mage[ifile].slit, float(mage[ifile].airmass), mage[ifile].exptype, mage[ifile].arcfile, mage[ifile].obj_id, format='(%"%15s %-17s%10s%7.3f%10s  %-15s %3d")'
     endif else begin
        print, mage[ifile].fitsfile, mage[ifile].object, mage[ifile].slit, float(mage[ifile].airmass), mage[ifile].exptype, mage[ifile].arcfile, format='(%"%15s %-17s%10s%7.3f%10s  %-15s")'
     endelse
  endfor
  ;;mwrfits, mage, "magestrct.fits", /create
RETURN
end
