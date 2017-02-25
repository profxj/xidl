;+ 
; NAME:
; imacsls_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of IMACS long slit frames.
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the package    
;
; CALLING SEQUENCE:
;  imacsls_strct, struct, LIST=, /MKDIR, /NOFILE, OUTFIL=, /NOEDIT
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;   struct   -  IMACS long slit IDL structure 
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list:  e.g.  'gd_files.lst'
;              Default is 'Raw/*3.fits*', 'Raw/*8.fits*'
;   /MKDIR      - Make directories
;   /NOEDIT     - Do not edit the hand
;   OUTFIL=     - Name of fits output file
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   imacsls_strct, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   08-Dec-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro imacsls_strct, struct, LIST=list, MKDIR=mkdir, NOFILE=nofile, $
                   OUTFIL=outfil, NOEDIT=noedit

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'imacsls_strct, struct, LIST=, MKDIR=, NOFILE=, NOLIST=, /NOEDIT (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( CCD ) then ccd = 'Site1'
  if not keyword_set( OUTFIL ) then outfil = 'imacslsstrct.fits'

; List
  if not keyword_set( LIST ) then begin
      img3 = findfile('Raw/*3.fits*',count=nimg) 
      img8 = findfile('Raw/*8.fits*',count=nimg) 
      img = [img3,img8]
      nimg = n_elements(img)
      if nimg EQ 0 then begin
          print, 'imacsls_strct: No images in Raw!'
          return
      endif
  endif else begin
      if x_chkfil(list) NE 1 then begin
          print, 'imacsls_strct: Trouble with filename ', list
          return
      endif
      readcol, list, img, FORMAT='A'
      nimg = n_elements(img)
  endelse


; Make directories
  if keyword_set( MKDIR ) then begin
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

  case ccd of 
      'Site1' : begin
          expcrd = 'EXPTIME'
          frmcrd = 'OBSNUM'
          utcrd = 'TIME'
          telescope = 'TELESCOP'
          pa = 'ROTANGLE'
          racrd = 'RA'
          deccrd = 'DEC'
          eqxcrd = 'EQUINOX'
          chipid = 'CHIP'
          objcrd = 'OBJECT'
          gangl = 'G-ANGLE'
          binning = 'BINNING'
          datecrd = 'UT-DATE'
          ccdspeed = 'SPEED'
          airmass = 'AIRMASS'
      end
      else : begin
          print, 'Not prepared for this ccd ', ccd
          return
      end
  endcase
          
  
;  Create the Structure
  tmp = { imacslsstrct}
  struct = replicate(tmp,nimg) 
  struct.ut = ' '
  struct.Obj = ' '
  struct.tel = ' '
; Default all files to be analysed
  struct.flg_anly = 1

;  Loop on Indv images
  for q=0,nimg-1 do begin
      print, 'Reading ', img[q]
      data = xmrdfits(img[q], 0, head, /silent, /fscale)

      ;; Frame
      spos = strpos(img[q], '.fits')
      struct[q].frame = long(strmid(img[q], spos-6, 4))

      ;; Parse the Header
      if keyword_set( expcrd ) then struct[q].exp = sxpar(head, expcrd)
      if keyword_set( racrd ) then struct[q].RA = sxpar(head, racrd)
      if keyword_set( deccrd ) then struct[q].DEC = sxpar(head, deccrd)
      if keyword_set( eqxcrd ) then struct[q].equinox = sxpar(head, eqxcrd)
      if keyword_set( utcrd ) then struct[q].UT = sxpar(head, utcrd)
;      if keyword_set( frmcrd ) then struct[q].frame = sxpar(head, frmcrd)
      if keyword_set( slitsz ) then struct[q].slit = sxpar(head, slitsz)
      if keyword_set( splitter ) then struct[q].splitter = sxpar(head, splitter)
      if keyword_set( airmass ) then struct[q].AM = sxpar(head, airmass)
      if keyword_set( objcrd ) then struct[q].Obj = strtrim(sxpar(head, objcrd),2)
      if keyword_set( telescope ) then $
        struct[q].TEL = strtrim(sxpar(head, telescope),2)
      if keyword_set( PA ) then struct[q].PA = strtrim(sxpar(head, PA),2)
      if keyword_set( GANGL ) then $
        struct[q].grangle = strtrim(sxpar(head, GANGL),2)

      ;; CCD
      struct[q].ccd = ccd
      struct[q].ccdspeed = sxpar(head, ccdspeed)
      case strtrim(struct[q].ccdspeed,2) of
          'Fast': begin
              struct[q].gain = 0.9
              struct[q].readno = 4.9
          end
          else: stop
      endcase

      case strtrim(sxpar(head,binning)) of 
          '1x1': begin
              struct[q].rbin = 1
              struct[q].cbin = 1
          end
          '2x2': begin
              struct[q].rbin = 2
              struct[q].cbin = 2
          end
          else: stop
      endcase
          
      ;; RED vs BLUE
      if sxpar(head,chipid) EQ 8 then struct[q].side = 1  $
      else struct[q].side = 2

      ;; Mode (set all to spec)
      struct[q].mode = 1

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
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
      struct[q].grising = ' '
  endfor

; Close the ASCII file
  if not keyword_set( NOFILE ) then imacsls_wrstrct, struct

; Edit
;  if not keyword_set( NOEDIT ) then imacsls_editstrct, struct

; Write the structure to FITS
  mwrfits, struct, outfil, /create

; All done
  print, 'imacsls_strct: All done!  Fits file in ', outfil

end
