pro mage_vachelio, obj_strct, head=head

   hdr = xheadfits(obj_strct[0].img_fil)

   date_day = float(strsplit(sxpar(hdr,"UT-DATE"),"-", /extract))
   date_time = float(strsplit(sxpar(hdr,"UT-TIME"),":", /extract))
;   end_time = strsplit(sxpar(hdr,"UT-END"),":", /extract)

   date = [date_day[0],date_day[1],date_day[2],date_time[0],date_time[1],date_time[2]]

   ra_deg = sxpar(hdr, "RA-D")
   dec_deg = sxpar(hdr, "DEC-D")

   juldate, date, jd  ; Geocentric JD-2400000
   jdhelio = helio_jd(jd, ra_deg, dec_deg) ;Returns heliocentric JD
   
   jdhelio += 2400000.d

   ; Heliocentric velocity correction in km/s
   helio_corr = (-1.0)*x_keckhelio(ra_deg, dec_deg, jd=jdhelio, obs='lco')

   gamma = sqrt( (1.d + helio_corr/299792.458d) / $
                 (1.d - helio_corr/299792.458d))

   for i=0, n_elements(obj_strct)-1 do begin

      ; Air to vacuum correct
      linwv = obj_strct[i].wave
      airtovac, linwv
      obj_strct[i].wave = linwv
      ; Geocentric to heliocentric correct
      obj_strct[i].wave *= gamma

   endfor

   if (keyword_set(head)) then begin
      arr = strsplit(obj_strct[0].img_fil,"/", /extract)
      filenm = arr[n_elements(arr)-1]
      msg = 'mage_pipe: Vac/helio applied for '+filenm+', helio = '+strtrim(helio_corr,2)+' km/s'
      sxaddpar, head, "COMMENT", msg
   endif

end
