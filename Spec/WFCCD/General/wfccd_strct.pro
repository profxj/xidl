;+ 
; NAME:
; wfccd_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates an outputs a structure for a series of WFCCD
;    spectroscopic images
;
; CALLING SEQUENCE:
;   
;  wfccd_strct, struct, [ccd], [tel], LIST=list, MKDIR=mkdir, NOFILE=nofile, NOLIST=nolist
;
; INPUTS:
;   [CCD]        - Set specific CCD header keywords
;                   Options: WFTek5 (default)
;   [tel]        - Set telescope
;                   Options: LCO-100 (default)
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  Creates an IDL structure for direct images 
;         -  ASCII file summarizing the structure
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list
;   MKDIR      - Make directories
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_strct, nght1_strct, 'WFTek5', 'LCO-100', /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Nov-2001 Written by JXP
;   18-Feb-2001 Modified by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_strct, struct, ccd, tel, LIST=list, MKDIR=mkdir, NOFILE=nofile, NOLIST=nolist, OUTFIL=outfil

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'wfccd_strct, struct, [ccd, tel], LIST=, MKDIR=, NOFILE=, NOLIST=, OUTFIL= (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( CCD ) then ccd = 'WFTek5'
  if not keyword_set( OUTFIL ) then outfil = 'struct.fits'
  if not keyword_set( TEL ) then tel = 'LCO-100'
  if not keyword_set( LIST ) then img = findfile('Raw/ccd*.fits*') $
  else readcol, list, img, format='A'

  nimg = n_elements(img)

  ; Output directories
  if keyword_set( MKDIR ) then begin
      a = findfile('pro/..', count=count)
      if count EQ 0 then file_mkdir, 'pro'
      a = findfile('Maps/..', count=count)
      if count EQ 0 then file_mkdir, 'Maps'
      a = findfile('Slits/..', count=count)
      if count EQ 0 then file_mkdir, 'Slits'
      a = findfile('Masks/..', count=count)
      if count EQ 0 then file_mkdir, 'Masks'
      a = findfile('OV/..', count=count)
      if count EQ 0 then file_mkdir, 'OV'
      a = findfile('Final/..', count=count)
      if count EQ 0 then file_mkdir, 'Final'
      a = findfile('Lists/..', count=count)
      if count EQ 0 then file_mkdir, 'Lists'
      a = findfile('Flats/..', count=count)
      if count EQ 0 then file_mkdir, 'Flats'
      a = findfile('Bias/..', count=count)
      if count EQ 0 then file_mkdir, 'Bias'
      a = findfile('Arcs/..', count=count)
      if count EQ 0 then file_mkdir, 'Arcs'
      a = findfile('Logs/..', count=count)
      if count EQ 0 then file_mkdir, 'Logs'
      a = findfile('Extract/..', count=count)
      if count EQ 0 then file_mkdir, 'Extract'
  endif

;  Header Keywords

  case tel of 
      'LCO-100': begin
          case ccd of 
              'WFTek5' : begin
                  expcrd = 'EXPTIME'
                  frmcrd = 'CCDPICNO'
                  gaincrd = 'GAIN'
                  racrd = 'RA'
                  deccrd = 'DEC'
                  eqxcrd = 'EQUINOX'
                  filtcrd1 = 'FILTER'
                  filtcrd2 = 'FILTERP'
                  apercrd1 = 'APERTURE'
                  apercrd2 = 'APERPPOS'
                  utcrd = 'UTSTART'
                  objcrd = 'OBJECT'
                  amcrd = 'AIRMASS'
                  casscrd = 'CASSPOS'
              end
              else : begin
                  print, 'Not prepared for this ccd ', ccd
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
  tmp = { wfccdstr }
  struct = replicate(tmp,nimg) 

;  Loop on Indv images

  for q=0,nimg-1 do begin
      print, 'Reading ', img[q]
      head = xheadfits(img[q], /silent)
      
; Default all files to be analysed

      struct[q].flg_anly = 1

;        Image stats

;      struct[q].med_raw = median(raw)
;      mmx = minmax(raw)
;      struct[q].min_raw = mmx[0]
;      struct[q].max_raw = mmx[1]

;        Parse the Header

      if keyword_set( expcrd ) then struct[q].exp = sxpar(head, expcrd)
      if keyword_set( frmcrd ) then struct[q].Frame = sxpar(head, frmcrd)
      if keyword_set( gaincrd ) then struct[q].gain = sxpar(head, gaincrd)
      if keyword_set( racrd ) then struct[q].RA = sxpar(head, racrd)
      if keyword_set( deccrd ) then struct[q].DEC = sxpar(head, deccrd)
      if keyword_set( eqxcrd ) then struct[q].Equinox = sxpar(head, eqxcrd)
      if keyword_set( utcrd ) then struct[q].UT = sxpar(head, utcrd)
      if keyword_set( objcrd ) then struct[q].Obj = sxpar(head, objcrd)
      if keyword_set( amcrd ) then struct[q].AM = sxpar(head, amcrd)
      if keyword_set( datecrd ) then begin
          date = sxpar(head, datecrd)
          struct[q].date = x_setjdate(date)
      endif

;        Filter and Grating
      if keyword_set( filtcrd1 ) then begin
          struct[q].filter = sxpar(head, filtcrd1)
         case struct[q].filter of
             'O': begin
                 struct[q].filter = 'C'
                 struct[q].grism = 'NO'
             end
             'RGRSM': begin
                 struct[q].filter = 'C'
                 struct[q].grism = 'RG'
             end
             'BGRSM': begin
                 struct[q].filter = 'C'
                 struct[q].grism = 'BG'
             end
             else: struct[q].grism = 'NONE'
         endcase
     endif

      if keyword_set( filtcrd2 ) then struct[q].filtpos = sxpar(head, filtcrd2)
      
;        Aperture
      if keyword_set( apercrd1 ) then begin
          dums = sxpar(head, apercrd1)
          if dums EQ 'OPEN' then struct[q].masknm = 'NO'
      endif
      if keyword_set( apercrd2 ) then struct[q].aperpos = sxpar(head, apercrd2)

;        CASS
      if keyword_set( casscrd ) then struct[q].casspos = sxpar(head, casscrd)

;        Telescope

      struct[q].tel = tel

;        CCD specfic (readnoise)

      struct[q].ccd = ccd
      case ccd of 
          'WFTek5' : begin
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
          else : print, ccd, 'not recognized, moving on..' ; Move along
      endcase

;        Image type

      case ccd of 
          'WFTek5' : begin
              imhead = sxpar(head, 'IMTYPE', count=count)
              if count EQ 0 then $
                imhead = sxpar(head, 'IMAGETYP', count=count)
              case imhead of 
                  'dark' : struct[q].type = 'DRK' ; Dark
                  'flat' : struct[q].type = 'FLT' ; Dome Flat
                  'bias' : struct[q].type = 'ZRO' ; Bias
                  else : begin
                      if ( struct[q].exp GT 600) then begin
                          struct[q].type = 'OBJ' ; Object
                          break
                      endif
                      if ( struct[q].exp LT 0.1) then begin
                          struct[q].type = 'ZRO' ; Bias
                          break
                      endif
                      struct[q].type = 'UNK' ; Unknown
                  end
              endcase
              if ( struct[q].exp LT 0.1) then struct[q].type = 'ZRO' ; Bias
          end
          else : print, ccd, 'not recognized, moving on..' ; Move along
      endcase
      
;        Set image names

  
      lslsh = strpos(img[q],'/',/reverse_search)
      if lslsh NE -1 then begin
          struct[q].img_root = strmid(img[q], lslsh+1) 
          ;; Truncate gz if it is there
          len = strlen(struct[q].img_root)
          if strmid(struct[q].img_root, len-1) EQ 'z' then $
            struct[q].img_root = strmid(struct[q].img_root,0,len-3)
          struct[q].rootpth = strmid(img[q], 0, lslsh+1)
      endif else struct[q].img_root = img[q]
      
;        Check for OV (unlikely to find it)
      
      ovfil = strjoin(['OV/','ov_',struct[q].img_root])
      a = findfile(ovfil,count=count)
      if count NE 0 then begin
          struct[q].flg_ov = 1
          struct[q].img_ov = ovfil
;          ovd = mrdfits(ovfil, 0, ohead)
;          struct[q].med_ov = median(ovd)
      endif else struct[q].img_ov = ' '
      
;        Check for Mask
      
      mskfil = strjoin(['Masks/','mk_',struct[q].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[q].flg_msk = 1
          struct[q].img_msk = a[0]
      endif else struct[q].img_msk = ' '

;        Check for Final
      
      mskfil = strjoin(['Final/','f_',struct[q].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[q].flg_final = 1
          struct[q].img_final = a[0]
      endif else struct[q].img_final = ' '
        
;      Zero out files
      struct[q].msk_fil = ' '
      struct[q].slit_fil = ' '
      struct[q].arc_fil = ' '
      struct[q].map_fil = ' '
      struct[q].flat_fil = ' '
      struct[q].mask_id = -1L
      struct[q].obj_fil = ' '

      
  endfor

; Close the ASCII file

  if not keyword_set( NOFILE ) then write_wfccdstr, struct

; Write the structure to FITS

  mwrfits, struct, outfil, /create

; All done
  print, 'wfccd_strct: All done!'

end
