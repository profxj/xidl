FUNCTION TSPEC_GET_JD, hdr, REDUCED = REDUCED, MODIFIED = MODIFIED
  
  utshut = sxpar(hdr, 'UTSHUT')
  temp = strsplit(utshut, 'T', /extract)
  utdate = temp[0]
  uttime = temp[1]

  date_day = float(strsplit(utdate, "-", /extract))
  date_time = float(strsplit(uttime, ":", /extract))
  date = [date_day[0], date_day[1], date_day[2], date_time[0] $
          , date_time[1], date_time[2]]

  if NOT KEYWORD_SET(REDUCED) AND NOT KEYWORD_SET(MODIFIED) then begin
     jd = jd + 2400000.0
  endif else if KEYWORD_SET(MODIFIED) then begin
     jd = jd - 0.5
  endif
  
  RETURN, jd
  
END


  RETURN
END
