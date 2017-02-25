;+ 
; NAME:
; kast_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of Kast frames.
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the Kast package    
;
; CALLING SEQUENCE:
;  kast_strct, struct, LIST=, /NOMKDIR, OUTFIL=, /NOEDIT
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;  struct  -- Kast IDL structure
;
; OPTIONAL KEYWORDS:
;   LIST=     - Image list:  e.g.  'gd_files.lst'
;               [Default is 'Raw/kast*.fits']
;   /NOMKDIR   - Do not make default directories
;   /NOEDIT    - Do not edit the hand
;   OUTFIL=    - Name of fits output file [default: 'kaststrct.fits']
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   kast_strct, nght1_strct, /MKDIR
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_strct, struct, LIST=list, NOMKDIR=nomkdir, $
                   OUTFIL=outfil, NOEDIT=noedit

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'kast_strct, struct, LIST=, MKDIR=, LIST=, OUTFIL=, /NOEDIT [v1.1]'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( CCD ) then ccd = 'Reticon1'
  if not keyword_set( OUTFIL ) then outfil = 'kaststrct.fits'

; List
  if not keyword_set( LIST ) then begin
      img = findfile('Raw/*.ccd*',count=nimg) 
      if nimg EQ 0 then begin
          print, 'kast_strct: No images in Raw!'
          return
      endif
  endif else begin
      if x_chkfil(list) NE 1 then begin
          print, 'kast_strct: Trouble with filename ', list
          return
      endif
      readcol, list, img, FORMAT='A'
      nimg = n_elements(img)
  endelse


; Make directories
  if not keyword_set( NOMKDIR ) then begin
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
      'Reticon1' : begin
          expcrd = 'EXPOSURE'
          frmcrd = 'OBSNUM'
          utcrd = 'TIME'
          ccdspeed = 'CLKSPEED'
          shutter = 'SHUTTER'
          side = 'SPSIDE'
          telescope = 'TELESCOP'
          pa = 'TUB'
          racrd = 'RA'
          deccrd = 'DEC'
          eqxcrd = 'EQUINOX'
          slitnm = 'SLITSIZE'
          objcrd = 'OBJECT'
          binning = 'BINNING'
          datecrd = 'UT-DATE'
          splitter = 'BSPLIT_N'
          slitsz = 'SLIT_P'
      end
      else : begin
          print, 'Not prepared for this ccd ', ccd
          return
      end
  endcase
          
  
;  Create the Structure
  tmp = { kaststrct}
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

      ;; Parse the Header
      if keyword_set( expcrd ) then struct[q].exp = sxpar(head, expcrd)
      if keyword_set( racrd ) then struct[q].RA = sxpar(head, racrd)
      if keyword_set( deccrd ) then struct[q].DEC = sxpar(head, deccrd)
      if keyword_set( eqxcrd ) then struct[q].equinox = sxpar(head, eqxcrd)
      if keyword_set( utcrd ) then struct[q].UT = sxpar(head, utcrd)
      if keyword_set( frmcrd ) then struct[q].frame = sxpar(head, frmcrd)
      if keyword_set( slitsz ) then struct[q].slit = sxpar(head, slitsz)
      if keyword_set( splitter ) then struct[q].splitter = sxpar(head, splitter)
      if keyword_set( objcrd ) then struct[q].Obj = strtrim(sxpar(head, objcrd),2)

      ;; CCD
      struct[q].ccd = ccd
      struct[q].ccdspeed = sxpar(head, ccdspeed)
      case strtrim(struct[q].ccdspeed,2) of
          'Slow': begin
              struct[q].gain = 3.8
              struct[q].readno = 6.
          end
          else: stop
      endcase

      ;; RED vs BLUE
      if strtrim(sxpar(head,side),2) EQ 'blue' then begin
          struct[q].side = 1 
          ;; MODE
          if strtrim(sxpar(head,'GRISM_N'),2) EQ 'open' then struct[q].mode = 0 $
          else begin
              struct[q].mode = 1
              struct[q].grising = strtrim(sxpar(head,'GRISM_N'),2)
          endelse
      endif else begin ;; RED
          struct[q].side = 2
          ;; MODE
          if strtrim(sxpar(head,'GRATNG_N'),2) EQ 'open' then struct[q].mode = 0 $
          else begin
              struct[q].mode = 1
              struct[q].grising = strtrim(sxpar(head,'GRATNG_N'),2)
          endelse
      endelse

      ;; TYPE
      if strtrim(sxpar(head,'SHUTTER'),2) NE 'OPEN' then begin
          if struct[q].exp LE 2 then struct[q].type = 'ZRO' $
          else struct[q].type = 'DRK'
      endif else begin
          if struct[q].mode EQ 0 then begin ;; IMAGE
              if strtrim(sxpar(head,'SLIT_N'),2) EQ 'open' then $
                struct[q].type = 'ACQ' $
              else struct[q].type = 'SLT'
          endif else begin
              if struct[q].exp LT 60. then begin
                  if median(data) GT 50. then struct[q].type = 'QTZ' $
                  else struct[q].type = 'STD'  ;; ARC is ignored!
              endif else struct[q].type = 'OBJ'  
          endelse
      endelse
      
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
  endfor

  ;; Create ASCII file
  kast_wrstrct, struct

; Edit
  if not keyword_set( NOEDIT ) then kast_editstrct, struct

; Write the structure to FITS
  mwrfits, struct, outfil, /create

; All done
  print, 'kast_strct: All done!  Fits file in ', outfil

end
