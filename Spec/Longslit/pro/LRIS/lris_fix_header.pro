pro lris_fix_header, planfile
;+
; LRIS sometimes fails to write the correct headers to the FITS
; files. This causes, in particular, LONG_REDUCE to crash when
; LONG_HELIO tries to parse the RA and Dec of the target.
;
; This routine goes through the plan.par file in the current
; directory, and ensures that each 'Science' file has the appropriate
; FITS headers. If not, it will write it into the errant file
; (assuming that all Science files in the directory are on the same target.
;-

  if not keyword_set(planfile) then planfile='plan.par'

  if not file_test(planfile) then $
     print, 'Error: plan.par file needs to be present. Try running LONG_PLAN'


  planstr = yanny_readone(planfile, hdr=planhdr, /anonymous)
   if keyword_set(in_planstr) then $  ; Input??
      planstr = in_planstr
   planhdr = strcompress(planhdr)
   if (NOT keyword_set(planstr)) then begin
       splog, 'Empty plan file ', planfile
       return
   endif

   ;; yanny_readone reads all the lines in plan.par, while the header
   ;; stores information at the top of the file
   planstr = yanny_readone(planfile, hdr=planhdr, /anonymous)
   indir = yanny_par(planhdr, 'indir')

   scifiles = planstr[where(strmatch(planstr.flavor, 'science', /fold_case) EQ 1 OR $
                           strmatch(planstr.flavor, 'std', /fold_case) EQ 1)].filename
   nsci = n_elements(scifiles)

   ;; loop over the Science files to see which don't have the
   ;; headers we want
   ;;
   ;; these arrays are to store indices of files that don't
   ;; have the corresponding headers
   no_tel = []
   no_mjd = []
   no_equinox = []
   no_ra = []
   no_dec = []
   
   for ii=0, nsci-1 do begin

      print, 'Fixing '+indir+scifiles[ii]
      
      hdr_tmp = headfits(indir+scifiles[ii],exten=0)
      tel_tmp       = sxpar(hdr_tmp, 'TELESCOP', count=telcount)
      mjd_tmp       = sxpar(hdr_tmp,'MJD-OBS', count=mjdcount)
      equinox_tmp   = sxpar(hdr_tmp, 'EQUINOX', count=eqcount)
      ra_tmp        = sxpar(hdr_tmp, 'RA', count=racount)
      dec_tmp       = sxpar(hdr_tmp, 'DEC',  count=deccount)

      if telcount EQ 0 then no_tel = [no_tel, ii]
      if mjdcount EQ 0     then no_mjd=[no_mjd, ii]
      if eqcount EQ 0 then no_equinox = [no_equinox, ii] else $
         equinox_all=equinox_tmp
      if racount EQ 0 then no_ra = [no_ra, ii] else $
         ra_all = ra_tmp
      if deccount EQ 0 then no_dec = [no_dec, ii] else $
         dec_all = dec_tmp
      
   endfor

   sci_nohdr = [no_tel, no_mjd, no_equinox, no_ra, no_dec]
   nohdr_flg = n_elements(sci_nohdr)

   if nohdr_flg EQ 0 then begin
      print, 'All Science files have FITS headers necessary for ' + $
             'heliocentric corrections'
      return
   endif else begin
      sci_nohdr = sci_nohdr[uniq(sci_nohdr[sort(sci_nohdr)])]
      n_nohdr = n_elements(sci_nohdr)

      for ii=0, n_nohdr-1 do begin
         itmp = sci_nohdr[ii]

         hdr_tmp = headfits(indir+scifiles[itmp],exten=0)

         tel_tmp       = sxpar(hdr_tmp, 'TELESCOP', count=telcount)
         mjd_tmp       = sxpar(hdr_tmp,'MJD-OBS', count=mjdcount)
         equinox_tmp   = sxpar(hdr_tmp, 'EQUINOX', count=eqcount)
         ra_tmp        = sxpar(hdr_tmp, 'RA', count=racount)
         dec_tmp       = sxpar(hdr_tmp, 'DEC',  count=deccount)

         if telcount EQ 0 then begin
            print, 'Adding TELESCOP to '+scifiles[itmp]
            sxaddpar, hdr_tmp, 'TELESCOP', 'Keck-I'
         endif
         
         if mjdcount EQ 0 then begin
            ;; Parse MJD from the file prefix
            yr = 2000. + float(strmid(scifiles[itmp],1,2))
            mth = float(strmid(scifiles[itmp],3,2))
            day = float(strmid(scifiles[itmp],5,2))
            print, 'Adding MJD-OBS to '+scifiles[itmp]
            datestr = strtrim(fix(yr),2)+'-'+string(fix(mth),'(i02)')+ '-'+$
                      string(fix(day),'(i02)')+' 08:00:00.00'

            mjdstr = date_conv(datestr, 'MODIFIED')
            sxaddpar, hdr_tmp, 'MJD-OBS', string(mjdstr, '(f9.3)')
         endif

         if eqcount EQ 0 then begin 
            print, 'Adding EQUINOX to '+scifiles[itmp]
            sxaddpar, hdr_tmp, 'EQUINOX', equinox_all
         endif

         if racount EQ 0 then begin
            print, 'Adding RA to '+scifiles[itmp]
            sxaddpar, hdr_tmp, 'RA', ra_all
         endif

         if deccount EQ 0 then begin
            print, 'Adding DEC to '+scifiles[itmp]
            sxaddpar, hdr_tmp, 'DEC', dec_all
         endif

         modfits, indir+scifiles[itmp], 0, hdr_tmp, exten_no=0
         
      endfor
      
   endelse
   
 
end
