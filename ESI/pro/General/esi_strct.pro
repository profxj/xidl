;+ 
; NAME:
; esi_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of ESI frames
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the ESI package    
;
; CALLING SEQUENCE:
;   
;  esi_strct, struct, LIST=, /MKDIR, /NOFILE, OUTFIL=, /NOEDIT
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  IDL structure 
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list:  e.g.  'gd_files.lst'
;              Default is 'Raw/esi*.fits'
;   MKDIR      - Make directories
;   NOEDIT     - Do not edit the hand
;   OUTFIL     - Name of fits output file
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   esi_strct, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   11-Jul-2002 Written by JXP
;   29-Jan-2003 Polished by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro esi_strct, struct, LIST=list, MKDIR=mkdir, NOFILE=nofile, $
                   OUTFIL=outfil, NOEDIT=noedit, IMG=img

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'esi_strct, struct, LIST=, MKDIR=, NOFILE=, NOLIST=, /NOEDIT (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( CCD ) then ccd = 'MIT-LL'
  if not keyword_set( OUTFIL ) then outfil = 'esistrct.fits'
  if not keyword_set( TEL ) then tel = 'KeckII'

; List
  if not keyword_set( IMG ) then begin
      if not keyword_set( LIST ) then begin
          img = findfile('Raw/esi*.fits*',count=nimg) 
          if nimg EQ 0 then begin
              print, 'esi_strct: No images in Raw!'
              return
          endif
      endif else begin
          if x_chkfil(list[0]) NE 1 then begin
              print, 'esi_strct: Trouble with filename ', list
              return
          endif
          readcol, list, img, FORMAT='A'
          nimg = n_elements(img)
      endelse
  endif else nimg = n_elements(img)


; Make directories
  if keyword_set( MKDIR ) then begin
      a = findfile('Maps/..', count=count)
      if count EQ 0 then file_mkdir, 'Maps'
      a = findfile('OV/..', count=count)
      if count EQ 0 then file_mkdir, 'OV'
      a = findfile('Final/..', count=count)
      if count EQ 0 then file_mkdir, 'Final'
      a = findfile('Flats/..', count=count)
      if count EQ 0 then file_mkdir, 'Flats'
      a = findfile('Bias/..', count=count)
      if count EQ 0 then file_mkdir, 'Bias'
      a = findfile('Arcs/..', count=count)
      if count EQ 0 then file_mkdir, 'Arcs'
      a = findfile('Arcs/TRC/..', count=count)
      if count EQ 0 then file_mkdir, 'Arcs/TRC'
      a = findfile('Logs/..', count=count)
      if count EQ 0 then file_mkdir, 'Logs'
      a = findfile('pro/..', count=count)
      if count EQ 0 then file_mkdir, 'pro'
      a = findfile('Sky/..', count=count)
      if count EQ 0 then file_mkdir, 'Sky'
      a = findfile('Extract/..', count=count)
      if count EQ 0 then file_mkdir, 'Extract'
      a = findfile('FSpec/..', count=count)
      if count EQ 0 then file_mkdir, 'FSpec'
      a = findfile('FSpec/txt/..', count=count)
      if count EQ 0 then file_mkdir, 'FSpec/txt'
  endif

;  Header Keywords

  case tel of 
      'KeckII': begin
          case ccd of 
              'MIT-LL' : begin
                  expcrd = 'EXPOSURE'
                  frmcrd = 'FRAMENO'
                  utcrd = 'UT'
                  gaincrd = 'CCDGAIN'
                  ccdspeed = 'CCDSPEED'
                  racrd = 'RA'
                  deccrd = 'DEC'
                  eqxcrd = 'EQUINOX'
                  rotmode = 'ROTMODE'
                  nampcrd = 'NUMAMPS'
                  slitnm = 'SLMSKNAM'
                  objcrd = 'OBJECT'
                  binning = 'BINNING'
                  amcrd = 'AIRMASS'
                  datecrd = 'DATE-OBS'
                  prismcrd = 'PRISMNAM'
                  imflatcrd = 'IMFLTNAM'
                  ldflatcrd = 'LDFLTNAM'
                  dwfiltcrd = 'DWFILNAM'
                  qtzcrd = 'LAMPQTZ1'
                  cuarcrd = 'LAMPCU1'
                  xecrd = 'LAMPNE1'
                  hgnecrd = 'LAMPAR1'
                  bincrd = 'BINNING'
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
          
  
;  Create the Structure
  tmp = { esistrct}
  struct = replicate(tmp,nimg) 
  struct.ut = ' '
; Default all files to be analysed
  struct.flg_anly = 1

;  Loop on Indv images
  for q=0,nimg-1 do begin
      print, 'Reading ', img[q]
      head = xheadfits(img[q], /silent)
      
      ;; Parse the Header

      if keyword_set( expcrd ) then struct[q].exp = sxpar(head, expcrd)
      if keyword_set( frmcrd ) then struct[q].frame = sxpar(head, frmcrd)
      if keyword_set( gaincrd ) then ccdgain = sxpar(head, gaincrd)
      if keyword_set( racrd ) then struct[q].RA = sxpar(head, racrd)
      if keyword_set( deccrd ) then struct[q].DEC = sxpar(head, deccrd)
      if keyword_set( eqxcrd ) then struct[q].equinox = sxpar(head, eqxcrd)
      if keyword_set( utcrd ) then struct[q].UT = sxpar(head, utcrd)
      ;; Name
      if keyword_set( objcrd ) then begin
          struct[q].Obj = sxpar(head, objcrd)
          struct[q].Obj = repchr(struct[q].Obj, '"',' ')
          struct[q].Obj = strcompress(struct[q].Obj, /remove_all)
      endif
      
      if keyword_set( amcrd ) then struct[q].AM = sxpar(head, amcrd)
      if keyword_set( ccdspeed ) then struct[q].ccdspeed = sxpar(head, ccdspeed)
      if keyword_set( datecrd ) then begin
          date = sxpar(head, datecrd)
          struct[q].date = x_setjdate(date, struct[q].UT)
      endif
      ;; MODE
      prism = sxpar(head, prismcrd)
      if strtrim(prism,2) EQ 'out' then struct[q].mode = 0 $
      else begin
          ldflat = sxpar(head,ldflatcrd)
          if strtrim(ldflat,2) EQ 'in' then struct[q].mode = 1 $
          else struct[q].mode = 2
      endelse
      ;; Binning
      if keyword_set( bincrd ) then begin
          binc = sxpar(head, bincrd)
          case strtrim(binc, 2) of 
              '1, 1': begin
                  struct[q].cbin = 4
                  struct[q].rbin = 1
              end
              '2, 2': begin
                  struct[q].cbin = 2
                  struct[q].rbin = 2
              end
              '4, 4': begin
                  struct[q].cbin = 4
                  struct[q].rbin = 4
              end
              else: stop
          endcase
      endif
      ;; SLIT
      if keyword_set( slitnm ) then begin
          slit = sxpar(head, slitnm)
          case slit of 
              '0.30_arcsec': struct[q].slit = 0.30
              '0.50_arcsec': struct[q].slit = 0.50
              '0.75_arcsec': struct[q].slit = 0.75
              '1.00_arcsec': struct[q].slit = 1.00
              '1.25_arcsec': struct[q].slit = 1.25
              '1.50_arcsec': struct[q].slit = 1.50
              '6.00_arcsec': struct[q].slit = 6.00
              'LowD_1.0': struct[q].slit = 1.00
              'LowD_0.75': struct[q].slit = 0.75
              'LowD_0.50': struct[q].slit = 0.50
              'MultiHoles': struct[q].slit = 9.99
              else:
          endcase
      endif
      ;; FILTER 
      if keyword_set( dwfiltcrd ) then begin
          imfilt = sxpar(head, dwfiltcrd)
         case strtrim(imfilt,2) of
             'Clear_S': struct[q].imfilt = 'C'
             'B': struct[q].imfilt = 'B'
             'V': struct[q].imfilt = 'V'
             'R': struct[q].imfilt = 'R'
             'I': struct[q].imfilt = 'I'
             else: 
         endcase
     endif

; LAMPS
     ;; QTZ
     if strtrim(sxpar(head,qtzcrd),2) EQ 'on' then struct[q].qtzlamp = 1
     ;; ARCS
     if strtrim(sxpar(head,cuarcrd),2) EQ 'on' then $  ; CuAr
       struct[q].arclamp = struct[q].arclamp + 1
     if strtrim(sxpar(head,xecrd),2) EQ 'on' then $  ; Xe
       struct[q].arclamp = struct[q].arclamp + 2
     if strtrim(sxpar(head,hgnecrd),2) EQ 'on' then $  ; HgNe
       struct[q].arclamp = struct[q].arclamp + 4

; AMPS
     if keyword_set( nampcrd ) then numamp = sxpar(head, nampcrd)
     struct[q].namp = numamp

     ;; Telescope
     struct[q].tel = tel

; CCD specfic (gain, rn)

     struct[q].ccd = ccd
     case ccd of 
         'MIT-LL' : begin
             case strtrim(ccdgain,2) of
                  'low' : begin
                      struct[q].gain = 1.29
                      struct[q].readno = 2.7
                  end
                  'high' : begin
                      struct[q].gain = 0.5
                      struct[q].readno = 2.7
                  end
                  else : stop
              endcase
          end
          else : stop
      endcase

; Image type

      type = sxpar(head, 'OBSTYPE')
      case strtrim(type,2) of 
          'Line' : struct[q].type = 'ARC' ; Dark
          'Dark' : struct[q].type = 'ZRO' ; Dome Flat
          'Bias' : struct[q].type = 'ZRO' ; Dome Flat
          'SkyFlat' : struct[q].type = 'TWI' ; Dome Flat
          'IntFlat' : struct[q].type = 'IFLT' ; Dome Flat
          'DmFlat' : struct[q].type = 'DFLT' ; Dome Flat
          'Object' : struct[q].type = 'OBJ' ; Bias
          else : begin
              if ( struct[q].exp GT 600) then begin
                  struct[q].type = 'OBJ' ; Object
                  break
              endif
              struct[q].type = 'UNK' ; Unknown
          end
      endcase
      ;; Catch zero frame
      if ( struct[q].exp LT 0.1) then struct[q].type = 'ZRO' ; Bias
      
;        Set image names

  
      lslsh = strpos(img[q],'/',/reverse_search)
      if lslsh NE -1 then begin
          struct[q].img_root = strmid(img[q], lslsh+1) 
          ;; Truncate gz if it is there
          len = strlen(struct[q].img_root)
          if strmid(struct[q].img_root, len-1) EQ 'z' then $
            struct[q].img_root = strmid(struct[q].img_root,0,len-3)
          ;; Root path
          struct[q].rootpth = strmid(img[q], 0, lslsh+1)
      endif else struct[q].img_root = img[q]
      
;        Check for OV (unlikely to find it)
      
      ovfil = strjoin(['OV/','ov_',struct[q].img_root])
      a = findfile(ovfil,count=count)
      if count NE 0 then begin
          struct[q].flg_ov = 1
          struct[q].img_ov = ovfil
      endif else struct[q].img_ov = ' '
      
;        Check for Final
      
      mskfil = strjoin(['Final/','f_',struct[q].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[q].flg_final = 1
          struct[q].img_final = a[0]
      endif else struct[q].img_final = ' '
        
;      Zero out files
      struct[q].ccdspeed = ' '
      struct[q].arc_fil = ' '
      struct[q].map_fil = ' '
      struct[q].flat_fil = ' '
      struct[q].obj_id = -1L
      struct[q].obj_fil = ' '
  endfor

; Close the ASCII file
  if not keyword_set( NOFILE ) then esi_wrstrct, struct

; Edit
  if not keyword_set( NOEDIT ) then esi_editstrct, struct

; Write the structure to FITS
  mwrfits, struct, outfil, /create

; All done
  print, 'esi_strct: All done!  Fits file in ', outfil

end
