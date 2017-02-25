;+ 
; NAME:
; uves_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of MIKE images.
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the MIKE package.  It will
;    attempt to identify the type of object using the algorithm
;    (uves_autoid).  It also parses the header for information related
;    to exposure time, RA,DEC, etc. Finally, the code calculates the
;    gain for pairs of Milky flat exposures with the routine
;    uves_gain.  At present these gain values are not used.
;    Use uves_setgain to set the values [highly recommended].
;
;    See uvesstrct__define.pro for all of the tags of the uves structure.
;
;
; CALLING SEQUENCE:
;   
;  uves_strct, struct, LIST=, /NOMKDIR, /NOFILE, OUTFIL=, /EDIT, FILE_LIST=
;
; INPUTS:
;    By default the program examines all of the files with the form
;    *b####.fits in the Raw/ directory.
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  IDL structure based on the MIKE images
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list:  e.g.  'gd_files.lst'
;   NOMKDIR    - Suppress the creation of sub directories (not
;                recommended)
;   EDIT     - Edit the final structure with the GUI
;   OUTFIL     - Name of fits output file (default = uvesstrct.fits)
;   FILE_LIST  - String Array of all input filenames
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;   Need to do some consistency checks between red and blue side,
;         i.e. they should be the same type
;
; EXAMPLES:
;   uves_strct, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro uves_strct, struct, LIST=list, NOMKDIR=nomkdir, NOFILE=nofile, $
                OUTFIL=outfil, EDIT=edit, FILE_LIST=file_list, $
                SILENT=silent, GETRA=getra

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'uves_strct, struct, LIST=, IMG=, /NOMKDIR, /NOFILE, ' + $
             'OUTFIL=, /EDIT, FILE_LIST=, /SILENT (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( TEL ) then tel = 'VLT2'
  if not keyword_set( CCD ) then ccd = 'BLUE-RED'
  if not keyword_set( OUTFIL ) then outfil = 'uves_strct.fits'


; Resolve the programs
  uves_rslvall

  if NOT keyword_set( FILE_LIST ) then begin
     
; List
    if not keyword_set( LIST ) then begin
      uves_file = '[0-9][0-9][0-9][0-9].fits'

;      img_archive = findfile('Raw/UV*'+uves_file+'*', count=narc) 
;      img_sarac = findfile('Raw/c*.fits*',count=ntele) 
;      img_saras = findfile('Raw/s*.fits*',count=ntele) 
      img = findfile('Raw/*.fits*',count=ntele) 
      
      notempty = where(img NE '',nimg)
      
      if nimg EQ 0 then begin
          print, 'uves_strct: No images in Raw!'
          return
      endif
      img=img[notempty]

    endif else begin
      if x_chkfil(list) NE 1 then begin
          print, 'uves_strct: Trouble with filename ', list
          return
      endif
      readcol, list, img, FORMAT='A'
    endelse
  endif else img = file_list

  nimg = n_elements(img)


; Make directories
  if not keyword_set( NOMKDIR ) then begin
      a = findfile('QA/..', count=count)
      if count EQ 0 then file_mkdir, 'QA'
;      a = findfile('Maps/..', count=count)
;      if count EQ 0 then file_mkdir, 'Maps'
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
      a = findfile('Arcs/Fits/..', count=count)
      if count EQ 0 then file_mkdir, 'Arcs/Fits'
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

  keywd = { uveskey }

  case ccd of 
      'BLUE-RED' : begin
          keywd.expcrd = 'EXPTIME'
          keywd.frmcrd = 'FRAMENO'
          keywd.utcrd = 'MJD-OBS'
          keywd.ccdrn = 'CCDRN'
          keywd.ccdgain = 'CCDGAIN'
          keywd.racrd = 'RA'
          keywd.deccrd = 'DEC'
          keywd.eqxcrd = 'EQUINOX'
          keywd.slitwid = 'SLITWID'
          keywd.objcrd = 'TARGNAME'
          keywd.binx = 'BINX'
          keywd.biny = 'BINY'
          keywd.prescan = 'PRESCAN'
          keywd.oscan = 'OSCAN'
          keywd.block = 'FIL1NAME'
          keywd.XDISPERS = 'XDISPERS'
          keywd.xdangl = 'WLEN1'
          keywd.echangl = 'ECHANGL'
          keywd.amcrd = 'AIRMASS'
          keywd.obstyp = 'OBSTYPE'
          keywd.hatch = 'HATOPEN'
          keywd.lamp = 'LAMPNAME'
          keywd.lampfil = 'LFILNAME'
          keywd.pane = 'PANELIST'
          keywd.ampmod = 'AMPMODE'
          keywd.mosmod = 'MOSMODE'
      end
      else : begin
          print, 'Not prepared for this ccd ', ccd
          return
      end
  endcase
          
  
;  Create the Structure
  tmp = { uvesstrct }
  struct = replicate(tmp, nimg)

;  Loop on Indv images
  for q=0,nimg-1 do begin
      if not keyword_set(SILENT) then print, 'Reading ', img[q]
      if keyword_set(HEAD2) then delvarx, head2

      ;; Set Image type
      print, img[q]
      phead = xheadfits(img[q],/silent)  ;; Primary header

      ;; Side
      mt = where(strmid(phead,0,21) EQ 'HIERARCH ESO INS MODE',nmt1)
      if strpos(phead[mt],'BLUE') GE 0 then struct[q].side = 1 $
      else struct[q].side = 2

      ;; Bias frames
      mt = where(strmid(phead,0,27) EQ 'HIERARCH ESO DET OUT1 PRSCX',nmt1)
      if nmt1 EQ 0 then begin
          struct[q].side = 2
          head2 = xheadfits(img[q],exten=1)
          struct[q].exten = 1  ;; Only doing the redder data for now
      endif

      if keyword_set(SARA) then dumf = img[q] 
      guess = uves_headid(phead, struct[q].side, SILENT=silent, SARA=dumf, $
                         H2=head2)

      
      ;; Angle
      struct[q].xdangl = sxpar(phead,keywd.xdangl)
      if struct[q].xdangl LE 0. then struct[q].xdangl = sxpar(head2,'WLEN2')
      if struct[q].xdangl GT 440. then struct[q].side = 2 else struct[q].side = 1

      ;; Create one structure entry per CCD
;      idx = uves_parsehdu(struct, phead, keywd)

      ;; DATE
      struct[q].date = sxpar(phead,keywd.utcrd) + 2400000.5D

      ;; Setup 
      struct[q].type = guess
      struct[q].exp = sxpar(phead,keywd.expcrd)
      struct[q].slitwid = strtrim(sxpar(phead,keywd.slitwid),2)
;      struct[q].block = strtrim(sxpar(phead,keywd.block),2)
      struct[q].lamp = strtrim(sxpar(phead,keywd.lamp),2)
;      struct[q].lampfil = strtrim(sxpar(phead,keywd.lampfil),2)
;      struct[q].cross = strtrim(sxpar(phead,keywd.xdispers),2)
      struct[q].echangl = sxpar(phead,keywd.echangl)
;      struct[q].frame = sxpar(phead,keywd.frmcrd)


      ;; Gain
      struct[q].gain = float(strtrim(sxpar(phead, keywd.ccdgain),2))
      struct[q].readno = float(strtrim(sxpar(phead, keywd.ccdrn),2))
;      uves_setgainrn, struct, idx, ccdgain

      struct[q].AM = sxpar(phead, keywd.amcrd)

      ;; RA, DEC
      if struct[q].type EQ 'OBJ' or keyword_set(GETRA) then begin
          struct[q].RA = sxpar(phead,keywd.racrd)
          struct[q].DEC= sxpar(phead,keywd.deccrd)
          struct[q].EQUINOX= sxpar(phead,keywd.eqxcrd)
          ;; Name
          spltn = strsplit(img[q], '_', /extract)
          struct[q].Obj = sxpar(phead,'OBJECT')
;          ps2 = strpos(spltn[2],'.fits')
;          struct[q].frame = long(strmid(spltn[2],0,ps2))
      endif else struct[q].Obj = 'CALIB'
      struct[q].frame = q+1

      ;; Gain, readnoise
      struct[q].ccd = ccd
      struct[q].tel = tel
;      struct[idx].ccdspeed = sxpar(phead, keywd.ccdspeed)

      ;; Binning
      struct[q].colbin= sxpar(phead,keywd.binx)
      struct[q].rowbin= sxpar(phead,keywd.biny)
      struct[q].poscan= [sxpar(phead,keywd.prescan), sxpar(phead,keywd.oscan)] 

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Set image names
;      lslsh = strpos(img[q],'/',/reverse_search)
;      if lslsh NE -1 then begin
;          ;; Look for uniq frame number
;          mtc = where(struct.frame EQ struct[idx[0]].frame, nmt)
;          if nmt GT n_elements(idx) then begin
;              for qq=0L,6 do begin
;                  ;; Print
;                  print, 'uves_strct: Incrementing frame # to avoid match!'
;                  dumf = struct[idx].frame + 1000L*(qq+1)
;                  mtc = where(struct.frame EQ dumf[0], nmt)
;                  if nmt EQ 0 then begin
;                      struct[idx].frame = dumf[0]
;                      break
;                  endif 
;              endfor
;          endif
;          ;; Base on frame number
;          cfm = strtrim(struct[idx[0]].frame,2)
;          for zz=0L,3-strlen(cfm) do cfm = '0'+cfm
;          struct[idx].frm_nam = 'uves'+cfm+'.fits'
;          struct[idx].img_root = strmid(img[q], lslsh+1) 
;          ;; Truncate gz if it is there
;          len = strlen(struct[idx[0]].img_root)
;          if strmid(struct[idx[0]].img_root, len-1) EQ 'z' then $
;            struct[idx].img_root = strmid(struct[idx].img_root,0,len-3)
;          ;; Root path
;          struct[idx].rootpth = strmid(img[q], 0, lslsh+1)
;      endif else struct[idx].img_root = img[q]

      struct[q].rootpth = 'Raw/'
      pos = strpos(img[q], '.gz')
      if pos LE 0 then pos = strlen(img[q])
      struct[q].img_root = strmid(img[q],4,pos-4)

      
        
;      Zero out files
      struct[q].arc_fil = ' '
      struct[q].arc_img = ' '
      struct[q].flat_fil = ' '
      struct[q].obj_id = -1L
      struct[q].obj_fil = ' '
      struct[q].img_ov = ' '
      struct[q].img_final = ' '
      struct[q].block = ' '
      struct[q].cross = ' '
      struct[q].lampfil = ' '
      struct[q].frm_nam = ' '
      struct[q].ampmode = ' '

  endfor

  ;; Drop the first (null) entry
;  struct = struct[1:*]

  ;; 
  struct.flg_anly = 1

; Close the ASCII file
  if not keyword_set( NOFILE ) then uves_wrstrct, struct

; Edit
  if keyword_set( EDIT ) then uves_editstrct, struct

; Write the structure to FITS
  if not keyword_set(NOFILE) then mwrfits, struct, outfil, /create

; All done
  if not keyword_set(SILENT) then $
    print, 'uves_strct: All done!  Fits file in ', outfil

end
