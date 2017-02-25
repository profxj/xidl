;+ 
; NAME:
; hamspec_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of HIRES images.
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the HIRES package.  It will
;    attempt to identify the type of object using the algorithm
;    (hamspec_autoid).  It also parses the header for information related
;    to exposure time, RA,DEC, etc. Finally, the code calculates the
;    gain for pairs of Milky flat exposures with the routine
;    hamspec_gain.  At present these gain values are not used.
;    Use hamspec_setgain to set the values [highly recommended].
;
;    See hamspecstrct__define.pro for all of the tags of the hamspec structure.
;
;
; CALLING SEQUENCE:
;   
;  hamspec_strct, struct, LIST=, /NOMKDIR, /NOFILE, OUTFIL=, /EDIT, FILE_LIST=
;
; INPUTS:
;    By default the program examines all of the files with the form
;    *b####.fits in the Raw/ directory.
;
; RETURNS:
;
; OUTPUTS:
;   struct     -  IDL structure based on the HIRES images
;
; OPTIONAL KEYWORDS:
;   LIST       - Image list:  e.g.  'gd_files.lst'
;   NOMKDIR    - Suppress the creation of sub directories (not
;                recommended)
;   EDIT     - Edit the final structure with the GUI
;   OUTFIL     - Name of fits output file (default = hamspecstrct.fits)
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
;   hamspec_strct, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hamspec_strct, struct, LIST=list, NOMKDIR=nomkdir, NOFILE=nofile, $
                   OUTFIL=outfil, EDIT=edit, FILE_LIST=file_list, $
                 SILENT=silent, ELAP=elap

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hamspec_strct, struct, LIST=, IMG=, /NOMKDIR, /NOFILE, ' + $
             'OUTFIL=, /EDIT, FILE_LIST=, /SILENT, /ELA	P (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( TEL ) then tel = 'Shane3m'
  if not keyword_set( OUTFIL ) then outfil = 'hamspec_strct.fits'


; Resolve the programs
  hamspec_rslvall

  if NOT keyword_set( FILE_LIST ) then begin
     
; List
    if not keyword_set( LIST ) then begin
      hamspec_file = '[0-9][0-9][0-9].fits'
;      hamspec_file = '[r,b][0-9][0-9][0-9][0-9].fits'

      img_archive = findfile('Raw/HI*'+hamspec_file+'*', count=narc) 
      img_tele = findfile('Raw/d*'+hamspec_file+'*',count=ntele) 
      img = [img_archive, img_tele]
      
      notempty = where(img NE '',nimg)
      
      if nimg EQ 0 then begin
          print, 'hamspec_strct: No images in Raw!'
          return
      endif
      img=img[notempty]

    endif else begin
      if x_chkfil(list) NE 1 then begin
          print, 'hamspec_strct: Trouble with filename ', list
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

  ;; Deal with dewar
  hd0 = xheadfits(strtrim(img[0],2))
  
  ccd = sxpar(hd0,'DSENSOR')
  if strlen(ccd) EQ 0 then ccd = sxpar(hd0,'CCDID')
  ccd = strcompress(ccd,/remove_all)

  ;;  Header Keywords
  keywd = { hamspeckey }

  case ccd of 
      'Loral2Kx2K' : begin
         keywd.expcrd = 'EXPTIME'
         keywd.frmcrd = 'OBSNUM'
         keywd.utcrd = 'DATE-OBS'
         keywd.mpmid = 'MP-MID'
;         keywd.ccdspeed = 'CCDSPEED'
;         keywd.ccdgain = 'CCDGAIN'
         keywd.racrd = 'RA'
         keywd.deccrd = 'DEC'
         keywd.eqxcrd = 'EQUINOX'
         keywd.decker = 'PLATENAM'
         keywd.dewarfoc = 'DFOCRAW'
         keywd.objcrd = 'OBJECT'
         keywd.rbin = 'RBIN'
         keywd.cbin = 'CBIN'
         keywd.block = 'DFILTNAM'
;         keywd.XDISPERS = 'XDISPERS'
         keywd.xdangl = 'GTILTRAW'
         keywd.echangl = 'DHEITRAW'
;         keywd.amcrd = 'AIRMASS'
         keywd.obstyp = 'OBSTYPE'
;         keywd.hatch = 'HATOPEN'
         keywd.lamp = 'LAMPPOS'
;         keywd.lampfil = 'LFILNAME'
;         keywd.pane = 'PANELIST'
;         keywd.ampmod = 'AMPMODE'
;         keywd.mosmod = 'MOSMODE'
      end
      'e2vCCD203-824kx4kthin' : begin
         keywd.expcrd = 'EXPTIME'
         keywd.frmcrd = 'OBSNUM'
         keywd.utcrd = 'DATE'
         keywd.mpmid = 'MP-MID'
;         keywd.ccdspeed = 'CCDSPEED'
;         keywd.ccdgain = 'CCDGAIN'
         keywd.racrd = 'RA'
         keywd.deccrd = 'DEC'
         keywd.eqxcrd = 'EQUINOX'
         keywd.decker = 'PLATENAM'
         keywd.dewarfoc = 'DFOCRAW'
         keywd.objcrd = 'OBJECT'
         keywd.rbin = 'RBIN'
         keywd.cbin = 'CBIN'
         keywd.block = 'DFILTNAM'
;         keywd.XDISPERS = 'XDISPERS'
         keywd.xdangl = 'GTILTRAW'
         keywd.echangl = 'DHEITRAW'
;         keywd.amcrd = 'AIRMASS'
         keywd.obstyp = 'OBSTYPE'
;         keywd.hatch = 'HATOPEN'
         keywd.lamp = 'LAMPPOS'
;         keywd.lampfil = 'LFILNAME'
;         keywd.pane = 'PANELIST'
;         keywd.ampmod = 'AMPMODE'
;         keywd.mosmod = 'MOSMODE'
      end
      else : begin
          print, 'Not prepared for this ccd ', ccd
          return
      end
   endcase

  
;  Create the Structure
  tmp = { hamspecstrct }
  struct = replicate(tmp,nimg)

;  Loop on Indv images
  for q=0,nimg-1 do begin
      if not keyword_set(SILENT) then print, 'Reading ', img[q]

      ;; Set Image type
      phead = xheadfits(img[q],/silent)  ;; Primary header
      guess = hamspec_headid(phead, keywd, SILENT=silent)

      ;; DATE
      date = (strsplit(sxpar(phead,keywd.utcrd),'T',/extract))[0]
      ut = sxpar(phead,keywd.mpmid)
      mjd = x_setjdate(date, ut)
      struct[q].date = mjd 

      ;; Setup 
      struct[q].type = guess
      struct[q].exp = sxpar(phead,keywd.expcrd)
      struct[q].block = strtrim(sxpar(phead,keywd.block),2)
      struct[q].lamp = strtrim(sxpar(phead,keywd.lamp),2)
;      struct[q].lampfil = strtrim(sxpar(phead,keywd.lampfil),2)
;      struct[q].cross = strtrim(sxpar(phead,keywd.xdispers),2)
      struct[q].xdangl = sxpar(phead,keywd.xdangl)
      struct[q].echangl = sxpar(phead,keywd.echangl)
      struct[q].frame = sxpar(phead,keywd.frmcrd)

      ;; Slit width and length
      decker = strsplit(strtrim(sxpar(phead,keywd.decker),2),':',/extract)
      struct[q].decker = strtrim(sxpar(phead,keywd.decker),2)
      struct[q].width = long(decker[0])
      struct[q].length = float(decker[1])

      ;; dewar focus
      struct[q].dewarfoc = long(sxpar(phead,keywd.dewarfoc))
      if guess eq 'TFLT' and  struct[q].dewarfoc GT 1800 then begin
         struct[q].type = 'MFLT'
         guess = 'MFLT'
      endif
      ;; Gain
;      ccdgain = strtrim(sxpar(phead, keywd.ccdgain),2)
;      hamspec_setgainrn, struct, idx, ccdgain
      ;; THESE NEED TO BE UPDATED WITH CCD, ETC.
      case ccd of 
         'Loral2Kx2K': begin
            struct[q].gain = 1.2
            struct[q].readno = 18.0
            struct[q].amp = 1 ;; N amps
            struct[q].ratio_gain = 1. ;; N amps
         end
         'e2vCCD203-824kx4kthin' : begin
            struct[q].gain = 0.9
            struct[q].readno = 2.5
            struct[q].amp = 2 ;; N amps
            struct[q].ratio_gain = 1. ;; N amps
         end
         else: stop
      endcase
      struct[q].ampmode =  ' '

      ;; RA, DEC
      struct[q].RA = sxpar(phead,keywd.racrd)
      struct[q].DEC= sxpar(phead,keywd.deccrd)
      struct[q].EQUINOX= sxpar(phead,keywd.eqxcrd)
      tobj = strtrim(sxpar(phead, keywd.objcrd),2)
      tobj = strjoin(strsplit(tobj, ' ', /extract))
      if strlen(tobj) EQ 0 then tobj = ' '
      struct[q].Obj = tobj

      ;; Calculate airmass
      ras = sxpar(phead, 'RA')
      decs = sxpar(phead, 'DEC')
      x_radec, ras, decs, rad, decd
      hangl = float(strsplit(sxpar(phead,'HA'),':',/extrac))
      if hangl[0] LT 0. then hangl = hangl[0]-hangl[1]/60. $
      else hangl = hangl[0]+hangl[1]/60. 
      airmass = airmass(40., decd, hangl)
      struct[q].AM = airmass

      ;; Gain, readnoise
;      struct[q].ccdspeed = sxpar(phead, keywd.ccdspeed)
      struct[q].ccd = ccd
      struct[q].tel = tel

      ;; Binning
;      binstring = strtrim(sxpar(phead,keywd.binning),2)
      case ccd of 
         'Loral2Kx2K' : begin
            struct[q].colbin = sxpar(phead,keywd.cbin)
            struct[q].rowbin = sxpar(phead,keywd.rbin)
         end
         'e2vCCD203-824kx4kthin' : begin
            struct[q].colbin = long(round(4160./sxpar(phead,'NAXIS1')))
            struct[q].rowbin = long(round(4096./sxpar(phead,'NAXIS2')))
         end
         else: stop
      endcase

      ;; Image name
      struct[q].img_root = img[q]

      
      ovfil = strjoin(['OV/','ov_',struct[q[0]].img_root])
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
      struct[q].arc_fil = ' '
      struct[q].arc_img = ' '
      struct[q].flat_fil = ' '
      struct[q].pflat_fil = ' '
      struct[q].obj_id = -1L
      struct[q].obj_fil = ' '
      struct[q].img_ov = ' '

      struct[q].rootpth = ' '
  endfor

  ;; 
  struct.flg_anly = 1

; Close the ASCII file
  if not keyword_set( NOFILE ) then hamspec_wrstrct, struct

; Edit
  if keyword_set( EDIT ) then hamspec_editstrct, struct

; Write the structure to FITS
  if not keyword_set(NOFILE) then mwrfits, struct, outfil, /create

; All done
  if not keyword_set(SILENT) then $
    print, 'hamspec_strct: All done!  Fits file in ', outfil

end
