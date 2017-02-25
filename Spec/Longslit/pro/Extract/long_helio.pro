PRO long_helio, scihdr, final_struct, UNDO=undo

nstr = n_elements(final_struct)
FOR ii=0,nstr-1 DO BEGIN
    wvo = final_struct[ii].WAVE_OPT
    wvb = final_struct[ii].WAVE_BOX
    ;; HELIO
    IF NOT KEYWORD_SET(NOHELIO) THEN BEGIN
        telescope = strcompress(strmid(sxpar(scihdr[*, 0], 'TELESCOP'), 0, 3) $
                                , /rem)
        IF NOT KEYWORD_SET(TELESCOPE) THEN TELESCOPE = $
          strcompress(strmid(sxpar(scihdr[*, 0], 'TELID'), 0, 3), /rem)

        if telescope eq '' and strcompress(sxpar(scihdr[*,0], 'INSTRUME'),/rem) eq 'OSIRIS' then telescope='GTC'
        IF telescope EQ 'ESO' THEN BEGIN
           IF strmatch(sxpar(scihdr[*, 0], 'TELESCOP'), '*VLT*') THEN telescope = 'VLT' $
           ELSE IF strmatch(sxpar(scihdr[*, 0], 'TELESCOP'), '*NTT*') THEN telescope = 'NTT' 
        ENDIF

        case telescope OF
            'Kec': begin  ;; Keck/LRIS
                mjd = double(sxpar(scihdr, 'MJD-OBS'))  + 2400000.5D
                equinox = double(sxpar(scihdr, 'EQUINOX'))
                ra = sxpar(scihdr, 'RA')
                dec = sxpar(scihdr, 'DEC')
                x_radec, ra, dec, radeg, decdeg
            end
            'Gem': BEGIN ;; Gemini GMOS
                if n_elements(size(hdr,/dime)) EQ 1 then idx = 0 else idx = 1
                mjd = double(sxpar(scihdr[*, idx], 'MJD-OBS'))  + 2400000.5D
                equinox = double(sxpar(scihdr[*, idx], 'EQUINOX'))
                radeg = double(sxpar(scihdr[*, 0], 'RA'))
                decdeg = double(sxpar(scihdr[*, 0], 'DEC'))
            END
            'mmt': begin ;; MMT
                date = strtrim(sxpar(scihdr, 'DATE-OBS'), 2)
                ut = sxpar(scihdr, 'UT')
                mjd = x_setjdate(date, ut)
                equinox = double(sxpar(scihdr, 'EPOCH'))
                obs = 'kpno'
                ra = sxpar(scihdr, 'RA')
                dec = sxpar(scihdr, 'DEC')
                x_radec, ra, dec, radeg, decdeg
             end
            'kp4': begin ;; KPNO 4-m
                date = strtrim(sxpar(scihdr, 'DATE-OBS'), 2)
                ut = sxpar(scihdr, 'UT')
                mjd = x_setjdate(date, ut)
                equinox = double(sxpar(scihdr, 'EPOCH'))
                obs = 'kpno'
                ra = sxpar(scihdr, 'RA')
                dec = sxpar(scihdr, 'DEC')
                x_radec, ra, dec, radeg, decdeg
            end
            'Baa': begin ;; Magellan
                date = strtrim(sxpar(scihdr, 'DATE-OBS'), 2)
                ut = sxpar(scihdr, 'UT-TIME')
                mjd = x_setjdate(date, ut)
                equinox = double(sxpar(scihdr, 'EQUINOX'))
                obs = 'lco'
                ra = sxpar(scihdr, 'RA')
                dec = sxpar(scihdr, 'DEC')
                x_radec, ra, dec, radeg, decdeg
            end
            'DuP': begin ;; DuPont
                date = strtrim(sxpar(scihdr, 'DATE-OBS'), 2)
                ut = sxpar(scihdr, 'UT-TIME')
                mjd = x_setjdate(date, ut)
                equinox = double(sxpar(scihdr, 'EQUINOX'))
                obs = 'lco'
                ra = sxpar(scihdr, 'RA')
                dec = sxpar(scihdr, 'DEC')
                x_radec, ra, dec, radeg, decdeg
            end
            '200': begin ;; P200
                date1 = strtrim(sxpar(scihdr, 'DATE'), 2)
                date = strmid(date1, 0, 10)
                ut = sxpar(scihdr, 'UT')
                mjd = x_setjdate(date, ut)
                equinox = 2000.0
                obs = 'Palomar'
                ra = sxpar(scihdr, 'RA')
                dec = sxpar(scihdr, 'DEC')
                x_radec, ra, dec, radeg, decdeg
             end
            'CA-': begin ;; P200
               mjd = sxpar(scihdr, 'JUL-DATE')
               equinox = 2000.0
               obs = 'ca'
               radeg = sxpar(scihdr, 'RA')
               decdeg = sxpar(scihdr, 'DEC')
            end
            '3.5': begin ;; APO 3.5
               date1 = strtrim(sxpar(scihdr, 'DATE-OBS'), 2)
               date = strmid(date1, 0, 10)
               ut = sxpar(scihdr, 'UTC-OBS')
               mjd = x_setjdate(date, ut)
               equinox = 2000.0
               obs = 'apo'
               ra = sxpar(scihdr, 'RA')
               dec = sxpar(scihdr, 'DEC')
               x_radec, ra, dec, radeg, decdeg
            end
            'GTC': begin ;; GTC
               mjd = sxpar(scihdr, 'MJD-OBS')+2400000.5D
               equinox = 2000.0
               obs = 'lapalma'
               ra = sxpar(scihdr, 'RA')
               dec = sxpar(scihdr, 'DEC')
               x_radec, ra, dec, radeg, decdeg
            end
            'LBT': begin
               mjd = sxpar(scihdr, 'MJD-OBS')+2400000.5D
               equinox = 2000.0
               obs = 'lbt'
               instrument = strcompress(sxpar(scihdr, 'INSTRUME'), /rem)
               CASE instrument OF
                  'LUCI2': BEGIN
                     radeg =  sxpar(scihdr, 'RA_DEG')
                     decdeg =  sxpar(scihdr, 'DEC_DEG')
                  END
                  ELSE: BEGIN 
                     ra = sxpar(scihdr, 'OBJRA')
                     dec = sxpar(scihdr, 'OBJDEC')
                     IF strmatch(ra, '*:*') EQ 0 OR strmatch(dec, '*:*') EQ 0 THEN BEGIN
                        ra = repstr(ra, ' ', ':')
                        dec = repstr(dec, ' ', ':')
                     ENDIF
                     x_radec, ra, dec, radeg, decdeg
                  END
               ENDCASE
            end
            'VLT': begin 
               mjd = sxpar(scihdr, 'MJD-OBS')+2400000.5D
               equinox = 2000.0
               obs = 'vlt'
               radeg = sxpar(scihdr, 'RA')
               decdeg = sxpar(scihdr, 'DEC')
               ;;x_radec, ra, dec, radeg, decdeg
            end
            'NTT': begin 
               mjd = sxpar(scihdr, 'MJD-OBS')+2400000.5D
               equinox = 2000.0
               obs = 'eso'
               radeg = sxpar(scihdr, 'RA')
               decdeg = sxpar(scihdr, 'DEC')
            end
            'WHT': begin 
               mjd = sxpar(scihdr, 'MJD-OBS')+2400000.5D
               equinox = 2000.0
               obs = 'lapalma'
               radh = sxpar(scihdr, 'RA')
               decd = sxpar(scihdr, 'DEC')
               x_radec, radh, decd, radeg, decdeg
            end
            else: begin
                ;;
                if (stregex(sxpar(scihdr,'INSTRUME'),'.*kast.*',$
                            /boolean,/fold_case) eq 1) or $
                  (stregex(sxpar(scihdr,'VERSION'),'kast*', $
                           /bool, /fold_case) EQ 1) then begin
                    totdate = strtrim(sxpar(scihdr, 'DATE-OBS'), 2)
                    date = strmid(totdate, 0, 10)
                    if (stregex(sxpar(scihdr,'VERSION'),'kast*', /bool, $
                                /fold_case) EQ 1) then $
                      ut = strmid(totdate,11) $
                    else ut = strtrim(sxpar(scihdr, 'TIME'), 2)

                    mjd = x_setjdate(date, ut)
                    obs = 'lick'
                    ra = sxpar(scihdr, 'RA')
                    dec = sxpar(scihdr, 'DEC')
                    x_radec, ra, dec, radeg, decdeg
                endif else stop
            end
         endcase
        ;; Correct
        if not keyword_set(UNDO) then SIGN = -1. else SIGN = 1.
        helio = (SIGN)*x_keckhelio(radeg, decdeg, equinox, jd = mjd, OBS = OBS)
        print, 'long_reduce: heliocentric correction :', $
               helio, ' km/s', format = '(a,f8.3,a)'
        hel_corr = sqrt( (1.d + helio/299792.458d) / $
                         (1.d - helio/299792.458d) )
        wvo = wvo * hel_corr
        wvb = wvb * hel_corr
    ENDIF
    ;; Save
    final_struct[ii].WAVE_OPT = wvo
    final_struct[ii].WAVE_BOX = wvb
ENDFOR

RETURN
END
