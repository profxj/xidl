;+ 
; NAME:
; xdimg_strct   
;     Version 1.2
;
; PURPOSE:
;    Creates and outputs a structure for a series of direct images.
;    This is a key routine in organizing a night of imaging data.
;
; CALLING SEQUENCE:
;   
;  xdimg_strct, struct, ccd, tel, LIST=, /MKDIR, /NOFILE
;
; INPUTS:
;   CCD        - CCD/INSTRUMENT;  Options: Tek5, LRISR, SITe1, SITe3, NIRC2
;   tel        - Telescope; Options: (Keck, LCO-40, LCO-100)
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  Creates an IDL structure for direct images 
;               +  an ASCII file summarizing the structure
;               +  a FITS file saving the IDL structure
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list [default is Raw/*.fits]
;   /NOFILE     - Do not create the ASCII file
;   /MKDIR      - Create a set of directories which facilitate the rest
;                    of the reduction process (e.g. Flats/ Bias/
;                    Photo/)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   xdimg_strct, nght1_strct, 'Tek5', /MKDIR
;
; PROCEDURES/FUNCTIONS CALLED:
;   MRDFITS
;   SXPAR
;   X_SETJDATE
;   MWRFITS
;   WRITE_DIMGSTR
;
; REVISION HISTORY:
;   18-July-2001 Written by JXP
;   22-Aug-2001  Allows for LRISR
;   24-Dec-2001  LRISR modifications
;   14-June-2007 LKP updated to allow for NIRC2 data
;-
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro xdimg_strct, struct, ccd, tel, LIST=list, MKDIR=mkdir, NOFILE=nofile, $
                 IMG=img, NOSTATS=nostats

;
  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'xdimg_strct, struct, ccd, tel, LIST=, MKDIR=, NOFILE= (v1.2)'
      return
  endif 

  COMMON SITE, lat, lng, tzone
  DRADEG = 180.d0/!dpi
  
;  Optional Keywords
  
  if not keyword_set( IMG ) then begin
      if not keyword_set( LIST ) then img = findfile('Raw/*.fits*') $
      else readcol, list, img, format='A'
  endif
  nimg = n_elements(img)

  if keyword_set( MKDIR ) then begin
     a = findfile('Masks/..', count=count)
     if count EQ 0 then file_mkdir, 'Masks', 'Masks/Sky'
     a = findfile('Final/..', count=count)
     if count EQ 0 then file_mkdir, 'Final'
     a = findfile('Lists/..', count=count)
     if count EQ 0 then file_mkdir, 'Lists'
     a = findfile('Photo/..', count=count)
     if count EQ 0 then file_mkdir, 'Photo'
     a = findfile('Flats/..', count=count)
     if count EQ 0 then file_mkdir, 'Flats'
     
     if ccd ne 'NIRC2' then begin
        a = findfile('OV/..', count=count)
        if count EQ 0 then file_mkdir, 'OV'
        a = findfile('Bias/..', count=count)
        if count EQ 0 then file_mkdir, 'Bias'
     endif else begin
        a = findfile('Darks/..', count=count)
        if count EQ 0 then file_mkdir, 'Darks', 'Darks/Sub'
        a = findfile('PixMask/..', count=count)
        if count EQ 0 then file_mkdir, 'PixMask'
     endelse
     
  endif
  
;  Header Keywords

  case tel of 
      'Keck': begin
          case ccd of 
              'LRISR' : begin
                  expcrd = 'ELAPTIME'
                  frmcrd = 'FRAMENO'
                  racrd = 'RA'
                  deccrd = 'DEC'
                  eqxcrd = 'EQUINOX'
                  filtcrd = 'REDFILT'
                  utcrd = 'UT'
                  objcrd = 'OBJECT'
                  amcrd = 'AIRMASS'
                  datecrd = 'DATE'
               end
              'NIRC2' : begin
                  expcrd = 'ITIME'
                  frmcrd = 'FRAMENO'
                  racrd = 'RA'  ;For NIRC2 these RA and DEC are given in decimal degrees.
                  deccrd = 'DEC'
                  eqxcrd = 'EQUINOX'
                  filtcrd = 'FWINAME'
                  utcrd = 'UTC'
                  objcrd = 'OBJECT'
                  amcrd = 'AIRMASS'
                  datecrd = 'DATE-OBS'
                  gaincrd = 'GAIN'
                  coaddscrd = 'COADDS'
                  sampmodecrd = 'SAMPMODE'
                   ; NOTE: The sampmode can either be 1, 2, or 3.  1 is never used.
                   ; 2 means correlated readouts at the beginning and the end of the exposure.
                   ; 3 means an another type of correlated readouts, even better than sampmode 2 if you're dominated by readnoise.
                  multisamcrd = 'MULTISAM'
                   ; NOTE: for a given sampmode type, you can have N number of reads.  This is the multisam number.
                  camnamecrd = 'CAMNAME'
                  naxis1crd = 'NAXIS1'
                  naxis2crd = 'NAXIS2'
                  shutter = 'SHRNAME'  ;*** lindsey
               end
              else : begin
                  print, 'Not prepared for this ccd ', ccd
                  return
              end
          endcase
      end
      'KPNO-4m': begin
          case ccd of 
              'MOSA' : begin
                  expcrd = 'EXPTIME'
                  frmcrd = 'FRAMENO'
                  racrd = 'RA'
                  deccrd = 'DEC'
                  eqxcrd = 'RADECEQ'
                  filtcrd = 'FILPOS'
                  utcrd = 'UT'
                  objcrd = 'OBJECT'
                  amcrd = 'AIRMASS'
                  datecrd = 'DATE-OBS'
              end
              else : begin
                  print, 'Not prepared for this ccd ', ccd
                  return
              end
          endcase
      end
      'LCO-100': begin
          case ccd of 
              'Tek5' : begin
                  expcrd = 'EXPTIME'
                  frmcrd = 'CCDPICNO'
                  gaincrd = 'GAIN'
                  racrd = 'RA'
                  deccrd = 'DEC'
                  eqxcrd = 'EQUINOX'
                  filtcrd = 'FILTER'
                  utcrd = 'UTSTART'
                  objcrd = 'OBJECT'
                  amcrd = 'AIRMASS'
              end
              else : begin
                  print, 'Not prepared for this ccd ', ccd
                  return
              end
          endcase
      end
      'LCO-40': begin
          case ccd of 
              'SITe1' : begin
                  expcrd = 'EXPTIME'
                  frmcrd = 'CCDPICNO'
                  gaincrd = 'GAIN'
                  filtcrd = 'FILTER'
                  utcrd = 'UTSTART'
                  objcrd = 'OBJECT'
                  datecrd = 'DATE-OBS'
                  racrd = 'RA'
                  deccrd = 'DEC'
              end   
              'SITe3' : begin
                  racrd = 'RA'
                  deccrd = 'DEC'
                  expcrd = 'EXPTIME'
                  frmcrd = 'CCDPICNO'
                  gaincrd = 'GAIN'
                  filtcrd = 'FILTER'
                  utcrd = 'UTSTART'
                  objcrd = 'OBJECT'
                  datecrd = 'DATE-OBS'
                  amcrd = 'AIRMASS'
              end   
              else : begin
                  print, 'CCD not defined!'
                  return
              end
          endcase
      end
      else : begin
          print, 'Telescope not defined!'
          return
      end
  endcase
          
  
;    
;  Create the Structure
  tmp = { distruct }
  struct = replicate(tmp,nimg) 

  struct.Obj = ' '
  struct.UT = ' '
  struct.RA = ' '
  struct.DEC = ' '
  struct.TEL = ' '
  struct.CCD = ' '
  struct.type = ' '
  struct.rootpth = ' '
  struct.img_root = ' '
  struct.img_ov = ' '
  struct.img_drk = ' '
  struct.img_msk = ' '
  struct.img_skymsk = ' '
  struct.img_final = ' '
  struct.img_drksub = ' '

;  Loop on Indv images

  for q=0,nimg-1 do begin
      print, 'Reading ', img[q]
      if not keyword_set( NOSTATS ) then $
      raw = xmrdfits(img[q], 0, head, /fscale, /silent) $
      else head = headfits(img[q])
      
; Default all files to be analysed

      struct[q].flg_anly = 1

;        Image stats

      if not keyword_set( NOSTATS ) then begin
          struct[q].med_raw = median(raw)
          mmx = minmax(raw)
          struct[q].min_raw = mmx[0]
          struct[q].max_raw = mmx[1]
      endif

;        Parse the Header

      if keyword_set( expcrd ) then struct[q].exp = sxpar(head, expcrd)
      if keyword_set( frmcrd ) then struct[q].Frame = sxpar(head, frmcrd)
      if keyword_set( gaincrd ) then struct[q].gain = sxpar(head, gaincrd)
      if keyword_set( racrd ) then begin
         if ccd eq 'NIRC2' then ra_deg=sxpar(head, racrd) else struct[q].RA = sxpar(head, racrd)
      endif
      if keyword_set( deccrd ) then begin
         if ccd eq 'NIRC2' then dec_deg=sxpar(head, deccrd) else struct[q].DEC = sxpar(head, deccrd)
      endif
      if ccd eq 'NIRC2' then begin
         coords=adstring(ra_deg, dec_deg, 3)
         coords=coords[0]
         ra_sex=strmid(coords,1,2)+':'+strmid(coords,4,2)+':'+strmid(coords,7,7)
         dec_sex=strmid(coords,16,3)+':'+strmid(coords,20,2)+':'+strmid(coords,23,7)
         struct[q].RA = ra_sex
         struct[q].DEC = dec_sex         
      endif
      if keyword_set( eqxcrd ) then struct[q].Equinox = sxpar(head, eqxcrd)
      if keyword_set( filtcrd ) then $
        struct[q].filter = strtrim(sxpar(head, filtcrd),2)
      if keyword_set( utcrd ) then struct[q].UT = strtrim(sxpar(head, utcrd),2)
      if keyword_set( objcrd ) then struct[q].Obj = strtrim(sxpar(head, objcrd),2)
      if keyword_set( amcrd ) then struct[q].AM = sxpar(head, amcrd)
      if keyword_set( datecrd ) then begin
          date = sxpar(head, datecrd)
          struct[q].date = x_setddate(strtrim(date,2))
      endif
      if keyword_set( coaddscrd ) then struct[q].coadds = sxpar(head, coaddscrd)
      if keyword_set( sampmodecrd ) then struct[q].sampmode = sxpar(head, sampmodecrd)
      if keyword_set( multisamcrd ) then struct[q].multisam = sxpar(head, multisamcrd)
      if keyword_set( naxis1crd ) then struct[q].naxis1 = sxpar(head, naxis1crd)
      if keyword_set( naxis2crd ) then struct[q].naxis2 = sxpar(head, naxis2crd)

;        Telescope

      struct[q].tel = tel

;        CCD specfic (readnoise)
      
;  NOTE: Before this point we have not distinguished between NIRC2
;  wide and narrow cameras because their header keywords are pretty
;  much the same.  But now we'll get into some specific stuff like
;  readnoise, so we'd better differentiate here.
      if ccd eq 'NIRC2' then begin
         if keyword_set( camnamecrd ) then begin
            camname = strtrim(sxpar(head, camnamecrd),2)
            if camname eq 'wide' then struct[q].ccd = 'NIRC2W'
            if camname eq 'narrow' then struct[q].ccd = 'NIRC2N'
            if camname eq 'medium' then struct[q].ccd = 'NIRC2M'
         endif
      endif else begin
         ccd = strtrim(ccd,2)
         struct[q].ccd = ccd
      endelse

      case struct[q].ccd of 
          'MOSA' : begin
              struct[q].gain = 2.5   ; MOSA
              struct[q].readno = 5.0
              struct[q].satur = 30000.
              struct[q].statsec = [0, 845, 0, 2047]
          end
          'LRISR' : begin
              struct[q].gain = 1.9   ; LHS
              struct[q].readno = 7.0
              struct[q].satur = 60000.
              struct[q].statsec = [0, 845, 0, 2047]
          end
          'NIRC2W' : begin
             struct[q].readno = 60.0  ;Lindsey - look up the right number here.
             struct[q].satur = 18000.  ;Lindsey - look up the right number here.
             struct[q].statsec = [534, 650, 394, 482]
          end
          'NIRC2N' : begin
             struct[q].readno = 60.0 ;Lindsey - look up the right number here.
             struct[q].satur = 18000. ;Lindsey - look up the right number here.
             struct[q].statsec = [534, 650, 394, 482]
          end
          'Tek5' : begin
              struct[q].satur = 25000.
              struct[q].statsec = [0, 2047, 0, 2046]
              case struct[q].gain of
                  1 : begin
                      struct[q].gain = 1.0
                      struct[q].readno = 5.6
                  end
                  2 : begin
                      struct[q].gain = 2.0
                      struct[q].readno = 6.6
                  end
                  3 : begin
                      struct[q].gain = 2.7
                      struct[q].readno = 7.0
                  end
                  else :
              endcase
          end
          'SITe1' : begin
              struct[q].statsec = [0, 2047, 0, 2046]
              struct[q].satur = 25000.
              case struct[q].gain of
                  3 : begin
                      struct[q].gain = 2.2
                      struct[q].readno = 7.0
                  end
                  else :
              endcase
          end
          'SITe3' : begin
              struct[q].statsec = [0, 2047, 0, 3148]
              struct[q].satur = 25000.
              case struct[q].gain of
                  3 : begin
                      struct[q].gain = 2.5
                      struct[q].readno = 6.6
                  end
                  else :
              endcase
          end
          else : print, 'xdimg_strct: ', ccd, $
            ' not recognized, moving on..' ; Move along
      endcase

;        Image type

      case struct[q].ccd of 
          'LRISR' : begin
              trap = sxpar(head, 'TRAPDOOR')
              shut = sxpar(head, 'AUTOSHUT')
              rotmod = sxpar(head, 'ROTMODE')
              ; Shutter is closed
              if shut EQ 0 then begin
                  if struct[q].exp LE 1 then struct[q].type = 'ZRO' $
                    else struct[q].type = 'DRK' ; Dark
                  break
              endif
              ; Rotator is not rotating
              if strmid(strtrim(rotmod,2),0,4) EQ 'stat' OR $
                strmid(strtrim(rotmod,2),0,4) EQ 'STAT' then begin
                  struct[q].type = 'DFT' ; Dome Flat
                  break
              endif
              ; Trapdoor closed
              if trap EQ 'closed' then begin
                  struct[q].type = 'UNK' ; Unknown
                  break
              endif
              ; Twilight/Std
              if struct[q].exp LE 60 then begin
                  if struct[q].med_raw GT 10000L then begin
                      struct[q].type = 'TWI' ; Twighlight Flat
                      break
                  endif else begin
                      struct[q].type = 'STD' ; Standard
                      break
                  endelse
              endif
              struct[q].type = 'OBJ' ; Object
          end
          'NIRC2W' : begin
             shut = sxpar(head, 'SHRNAME')
             domeaz = sxpar(head, 'DOMEPOSN')  ;dome's azimuth position
             telaz = sxpar(head, 'AZ')  ; telescope's azimuth position
             telaz = (telaz + 360.) mod 360.
             domeaz = (domeaz + 360.) mod 360.
             observatory, 'keck', obs_struct
             tzone = obs_struct.tz
             lng = 360.d0 - obs_struct.longitude
             lat = obs_struct.latitude
             jd = 2400000.5D + sxpar(head, 'MJD-OBS')
             sunpos, jd, ra1, dec1           ; returns degrees
             zenpos, jd, ra2, dec2           ; returns radians
             ra2 = ra2 * DRADEG
             dec2 = dec2 * DRADEG
             sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
             ; Darks is shutter is closed
             if strtrim(shut,2) eq 'closed' then begin
                struct[q].type = 'DRK'
                break
             endif
             ; Domeflat if shutter's open and telescope isn't pointing out of dome slit.
             if abs(domeaz - telaz) GT 4. then begin
                struct[q].type = 'DFT'
                break
             endif
             ; Twilight if shutter's open, telescope points out of dome, 
             ; but sun angle is greater than -12 degrees.
             if sun_angle GT -12. then begin
                struct[q].type = 'TWI'
                break
             endif
             ; Standard frame if shutter's open, telescope points out of dome,
             ; and sun is lower than -12 deg below horizon, but exposure time is LT 3 seconds.
             if sxpar(head, 'ITIME') LT 3 then begin
                struct[q].type = 'STD'
             endif else begin
                struct[q].type = 'OBJ'
             endelse
          end
          'NIRC2N' : begin
             shut = sxpar(head, 'SHRNAME')
             domeaz = sxpar(head, 'DOMEPOSN') ;dome's azimuth position
             telaz = sxpar(head, 'AZ')        ; telescope's azimuth position
             telaz = (telaz + 360.) mod 360.
             domeaz = (domeaz + 360.) mod 360.
             observatory, 'keck', obs_struct
             tzone = obs_struct.tz
             lng = 360.d0 - obs_struct.longitude
             lat = obs_struct.latitude
             jd = 2400000.5D + sxpar(head, 'MJD-OBS')
             sunpos, jd, ra1, dec1           ; returns degrees
             zenpos, jd, ra2, dec2           ; returns radians
             ra2 = ra2 * DRADEG
             dec2 = dec2 * DRADEG
             sun_angle = 90. - djs_diff_angle(ra1, dec1, ra2, dec2)
             ; Darks is shutter is closed
             if strtrim(shut,2) eq 'closed' then begin
                struct[q].type = 'DRK'
                break
             endif
             ; Domeflat if shutter's open and telescope isn't pointing out of dome slit.
             if abs(domeaz - telaz) GT 4. then begin
                struct[q].type = 'DFT'
                break
             endif
             ; Twilight if shutter's open, telescope points out of dome, 
             ; but sun angle is greater than -12 degrees.
             if sun_angle GT -12. then begin
                struct[q].type = 'TWI'
                break
             endif
             ; Standard frame if shutter's open, telescope points out of dome,
             ; and sun is lower than -12 deg below horizon, but exposure time is LT 3 seconds.
             if sxpar(head, 'ITIME') LT 3 then begin
                struct[q].type = 'STD'
             endif else begin
                struct[q].type = 'OBJ'
             endelse
          end
          'Tek5' : begin
              imhead = sxpar(head, 'IMTYPE')
              case imhead of 
                  'dark' : struct[q].type = 'DRK' ; Dark
                  'flat' : struct[q].type = 'DFT' ; Dome Flat
                  'bias' : struct[q].type = 'ZRO' ; Bias
                  else : begin
                      if ( struct[q].exp GT 180. AND abs(float(struct[q].med_raw- $
                                                 struct[q].min_raw)) / $
                           float(struct[q].med_raw) LT 0.05) then begin
                          struct[q].type = 'DRK' ; Dark
                          break
                      endif
                      if ( struct[q].exp LT 0.1) then begin
                          struct[q].type = 'ZRO' ; Bias
                          break
                      endif
                      if ( struct[q].exp LT 60. AND struct[q].med_raw GT 6000 ) then begin
                          struct[q].type = 'TWI' ; Twighlight Flat
                          break
                      endif
                      if ( struct[q].exp LT 120 AND struct[q].med_raw LT 1200 ) then $
                        struct[q].type = 'STD' $ ; Standard star field 
                      else struct[q].type = 'OBJ'  ; Object
                  end
              endcase
          end
          'SITe1' : begin
              imhead = sxpar(head, 'IMTYPE')
              case imhead of 
                  'dark' : struct[q].type = 'DRK' ; Dark
                  'flat' : struct[q].type = 'DFT' ; Dome Flat
                  'bias' : struct[q].type = 'ZRO' ; Bias
                  else : begin
                      if ( struct[q].exp GT 180. AND abs(float(struct[q].med_raw- $
                                                 struct[q].min_raw)) / $
                           float(struct[q].med_raw) LT 0.05) then begin
                          struct[q].type = 'DRK' ; Dark
                          break
                      endif
                      if ( struct[q].exp LT 0.1) then begin
                          struct[q].type = 'ZRO' ; Bias
                          break
                      endif
                      if ( struct[q].exp LT 60. AND struct[q].med_raw GT 6000 ) then begin
                          struct[q].type = 'TWI' ; Twighlight Flat
                          break
                      endif
                      if ( struct[q].exp LT 90 AND struct[q].med_raw LT 3000 ) then $
                        struct[q].type = 'STD' $ ; Standard star field 
                      else struct[q].type = 'OBJ'  ; Object
                  end
              endcase
          end
          'SITe3' : begin
              imhead = sxpar(head, 'IMAGETYP')
              case imhead of 
                  'dark' : struct[q].type = 'DRK' ; Dark
                  'flat' : struct[q].type = 'DFT' ; Dome Flat
                  'bias' : struct[q].type = 'ZRO' ; Bias
                  else : begin
                      if ( struct[q].exp GT 180. AND abs(float(struct[q].med_raw- $
                                                 struct[q].min_raw)) / $
                           float(struct[q].med_raw) LT 0.05) then begin
                          struct[q].type = 'DRK' ; Dark
                          break
                      endif
                      if ( struct[q].exp LT 0.1) then begin
                          struct[q].type = 'ZRO' ; Bias
                          break
                      endif
                      if ( struct[q].exp LT 60. AND struct[q].med_raw GT 6000 ) then begin
                          struct[q].type = 'TWI' ; Twighlight Flat
                          break
                      endif
                      if ( struct[q].exp LT 90 AND struct[q].med_raw LT 3000 ) then $
                        struct[q].type = 'STD' $ ; Standard star field 
                      else struct[q].type = 'OBJ'  ; Object
                  end
              endcase
          end
          else : print, 'xdimg_strct (image type): ', ccd, $
            ' not recognized, moving on..' ; Move along
      endcase
      
;        Set image names

      lslsh = strpos(img[q],'/',/reverse_search)
      if lslsh NE -1 then begin
          struct[q].img_root = strmid(img[q], lslsh+1) 
          ;; Truncate gz if it is there
          len = strlen(struct[q].img_root)
          if strmid(struct[q].img_root, len-1) EQ 'z' then $
            struct[q].img_root = strmid(struct[q].img_root,0,len-3)
          ;; Path
          struct[q].rootpth = strmid(img[q], 0, lslsh+1)
      endif else struct[q].img_root = img[q]
      
;        Check for OV (unlikely to find it)
      
      ovfil = strjoin(['OV/','ov_',struct[q].img_root])
      a = findfile(ovfil,count=count)
      if count NE 0 then begin
          struct[q].flg_ov = 1
          struct[q].img_ov = ovfil
          ovd = xmrdfits(ovfil, 0, ohead, /fscale, /silent)
          struct[q].med_ov = median(ovd)
      endif else struct[q].img_ov = ' '
      
;        Check for Mask
      
      mskfil = strjoin(['Masks/','mk_',struct[q].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[q].flg_msk = 1
          struct[q].img_msk = a[0]
      endif else struct[q].img_msk = ' '

;        Check for Sky Mask
      
      mskfil = strjoin(['Masks/Sky/','sm_',struct[q].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[q].flg_skymsk = 1
          struct[q].img_skymsk = a[0]
      endif else struct[q].img_skymsk = ' '
      
;        Check for Final
      
      mskfil = strjoin(['Final/','f_',struct[q].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[q].flg_final = 1
          struct[q].img_final = a[0]
      endif else struct[q].img_final = ' '
              
;       No darks or coadd divided images have been created
      struct[q].img_drk = ' '
      struct[q].img_drksub = ' '

  endfor

; Close the ASCII file

  if not keyword_set( NOFILE ) then write_dimgstr, struct

; Write the structure to FITS

  mwrfits, struct, 'struct.fits', /create

end
