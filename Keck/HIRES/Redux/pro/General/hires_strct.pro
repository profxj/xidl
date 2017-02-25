;+ 
; NAME:
; hires_strct   
;     Version 1.1
;
; PURPOSE:
;    Creates and outputs a structure for a series of HIRES images.
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the HIRES package.  It will
;    attempt to identify the type of object using the algorithm
;    (hires_autoid).  It also parses the header for information related
;    to exposure time, RA,DEC, etc. Finally, the code calculates the
;    gain for pairs of Milky flat exposures with the routine
;    hires_gain.  At present these gain values are not used.
;    Use hires_setgain to set the values [highly recommended].
;
;    See hiresstrct__define.pro for all of the tags of the hires structure.
;
;
; CALLING SEQUENCE:
;   
;  hires_strct, struct, LIST=, /NOMKDIR, /NOFILE, OUTFIL=, /EDIT, FILE_LIST=
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
;   OUTFIL     - Name of fits output file (default = hiresstrct.fits)
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
;   hires_strct, nght1_strct
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   24-Jan-2005 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_strct, struct, LIST=list, NOMKDIR=nomkdir, NOFILE=nofile, $
                 OUTFIL=outfil, EDIT=edit, FILE_LIST=file_list, $
                 SILENT=silent, ELAP=elap, CCD=ccd

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'hires_strct, struct, LIST=, IMG=, /NOMKDIR, /NOFILE, ' + $
             'OUTFIL=, /EDIT, FILE_LIST=, /SILENT, /ELA, CCD=	 (v1.1)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( TEL ) then tel = 'KeckI'
  if not keyword_set( CCD ) then ccd = 'MOSAIC'
  if not keyword_set( OUTFIL ) then outfil = 'hiresstrct.fits'


; Resolve the programs
  hires_rslvall

  if NOT keyword_set( FILE_LIST ) then begin
     
; List
    if not keyword_set( LIST ) then begin
      hires_file = '[0-9][0-9][0-9][0-9].fits'
;      hires_file = '[r,b][0-9][0-9][0-9][0-9].fits'

      img_marcy = findfile('Raw/j*'+hires_file+'*', count=narc) 
      img_archive = findfile('Raw/HI*'+hires_file+'*', count=narc) 
      img_tele = findfile('Raw/hires'+hires_file+'*',count=ntele) 
      img = [img_archive, img_tele, img_marcy]
      
      notempty = where(img NE '',nimg)
      
      if nimg EQ 0 then begin
          print, 'hires_strct: No images in Raw!'
          return
      endif
      img=img[notempty]

    endif else begin
      if x_chkfil(list) NE 1 then begin
          print, 'hires_strct: Trouble with filename ', list
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

  keywd = { hireskey }

  case ccd of 
      'MOSAIC' : begin
;          keywd.expcrd = 'EXPTIME'  ;; Changed from default on Jan 18, 2008
          keywd.expcrd = 'ELAPTIME'
          keywd.frmcrd = 'FRAMENO'
          keywd.utcrd = 'MJD'
          keywd.ccdspeed = 'CCDSPEED'
          keywd.ccdgain = 'CCDGAIN'
          keywd.racrd = 'RA'
          keywd.deccrd = 'DEC'
          keywd.eqxcrd = 'EQUINOX'
          keywd.decker = 'DECKNAME'
          keywd.objcrd = 'TARGNAME'
          keywd.binning = 'BINNING'
          keywd.block = 'FIL1NAME'
          keywd.XDISPERS = 'XDISPERS'
          keywd.xdangl = 'XDANGL'
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
      'SINGLE' : begin  ;; Old-school detector (Tektronix)
          keywd.expcrd = 'ELAPTIME'
          keywd.frmcrd = 'FRAMENO'
          keywd.utcrd = 'MJD'
          keywd.ccdspeed = 'CCDSPEED'
          keywd.ccdgain = 'CCDGAIN'
          keywd.racrd = 'RA'
          keywd.deccrd = 'DEC'
          keywd.eqxcrd = 'EQUINOX'
          keywd.decker = 'DECKNAME'
          keywd.objcrd = 'TARGNAME'
          keywd.binning = 'BINNING'
          keywd.block = 'FIL1NAME'
          keywd.XDISPERS = 'XDISPERS'
          keywd.xdangl = 'XDANGL'
          keywd.echangl = 'ECHANGL'
          keywd.amcrd = 'AIRMASS'
          keywd.obstyp = 'OBSTYPE'
          keywd.hatch = 'HATOPEN'
          keywd.lamp = 'LAMPNAME'
          keywd.lampfil = 'LFILNAME'
          keywd.ampmod = 'AMPMODE'
          ;keywd.mosmod = 'MOSMODE'
      end
      else : begin
          print, 'Not prepared for this ccd ', ccd
          return
      end
  endcase
          
  
;  Create the Structure
  tmp = { hiresstrct }
  struct = [tmp]

;  Loop on Indv images
  for q=0,nimg-1 do begin
      if not keyword_set(SILENT) then print, 'Reading ', img[q]

      ;; Set Image type
      phead = xheadfits(img[q],/silent)  ;; Primary header
      guess = hires_headid(phead, keywd, SILENT=silent)

      ;; Create one structure entry per CCD
      if strmatch(ccd, 'MOSAIC') then $ 
         idx = hires_parsehdu(struct, phead, keywd) $
      else begin ;; Original CCD
         struct = [struct, {hiresstrct}]
         struct.ampmode = ' '
         idx = q+1
         struct.chip = -1
      endelse

      ;; DATE
      struct[idx].date = sxpar(phead,keywd.utcrd) + 2400000.5D

      ;; Setup 
      struct[idx].type = guess
      struct[idx].exp = sxpar(phead,keywd.expcrd)
      struct[idx].decker = strtrim(sxpar(phead,keywd.decker),2)
      struct[idx].block = strtrim(sxpar(phead,keywd.block),2)
      struct[idx].lamp = strtrim(sxpar(phead,keywd.lamp),2)
      struct[idx].lampfil = strtrim(sxpar(phead,keywd.lampfil),2)
      struct[idx].cross = strtrim(sxpar(phead,keywd.xdispers),2)
      struct[idx].xdangl = sxpar(phead,keywd.xdangl)
      struct[idx].echangl = sxpar(phead,keywd.echangl)
      struct[idx].frame = sxpar(phead,keywd.frmcrd)

      ;; Gain
      ccdgain = strtrim(sxpar(phead, keywd.ccdgain),2)
      if strmatch(ccd,'SINGLE') then ccdgain = 'SINGLE'
      hires_setgainrn, struct, idx, ccdgain

      ;; RA, DEC
      struct[idx].RA = sxpar(phead,keywd.racrd)
      struct[idx].DEC= sxpar(phead,keywd.deccrd)
      struct[idx].EQUINOX= sxpar(phead,keywd.eqxcrd)
      tobj = strtrim(sxpar(phead, keywd.objcrd),2)
      tobj = strjoin(strsplit(tobj, ' ', /extract))
      if strlen(tobj) EQ 0 then tobj = ' '
      struct[idx].Obj = tobj
      struct[idx].AM = sxpar(phead, keywd.amcrd)

      ;; Gain, readnoise
      struct[idx].ccd = ccd
      struct[idx].tel = tel
;      struct[idx].ccdspeed = sxpar(phead, keywd.ccdspeed)

      ;; Binning
      binstring = strtrim(sxpar(phead,keywd.binning),2)
      case binstring of
          '1,1': begin
              struct[idx].colbin = 1
              struct[idx].rowbin = 1
          end
          '2,1': begin
              struct[idx].colbin = 2
              struct[idx].rowbin = 1
          end
          '1,2': begin
              struct[idx].colbin = 1
              struct[idx].rowbin = 2
          end
          '2,2': begin
              struct[idx].colbin = 2
              struct[idx].rowbin = 2
          end
          '3,1': begin
              struct[idx].colbin = 3
              struct[idx].rowbin = 1
          end
          else: stop ;begin
;            struct[q].colbin = long(strmid(binstring,0,1))
;            struct[q].rowbin = long(strmid(binstring,strpos(binstring,',')+1,1))
;            print, 'Grabbing binning from header.. ', $
;               strcompress(string(struct[q].colbin, ' x', struct[q].rowbin)) 
;          end
      endcase

      
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;        Set image names
      lslsh = strpos(img[q],'/',/reverse_search)
      if lslsh NE -1 then begin
          ;; Look for uniq frame number
          mtc = where(struct.frame EQ struct[idx[0]].frame, nmt)
          if nmt GT n_elements(idx) then begin
              for qq=0L,6 do begin
                  ;; Print
                  print, 'hires_strct: Incrementing frame # to avoid match!'
                  dumf = struct[idx].frame + 1000L*(qq+1)
                  mtc = where(struct.frame EQ dumf[0], nmt)
                  if nmt EQ 0 then begin
                      struct[idx].frame = dumf[0]
                      break
                  endif 
              endfor
          endif
          ;; Base on frame number
          cfm = strtrim(struct[idx[0]].frame,2)
          for zz=0L,3-strlen(cfm) do cfm = '0'+cfm
          struct[idx].frm_nam = 'hires'+cfm+'.fits'
          struct[idx].img_root = strmid(img[q], lslsh+1) 
          ;; Truncate gz if it is there
          len = strlen(struct[idx[0]].img_root)
          if strmid(struct[idx[0]].img_root, len-1) EQ 'z' then $
            struct[idx].img_root = strmid(struct[idx].img_root,0,len-3)
          ;; Root path
          struct[idx].rootpth = strmid(img[q], 0, lslsh+1)
      endif else struct[idx].img_root = img[q]

      
      ovfil = strjoin(['OV/','ov_',struct[idx[0]].img_root])
      a = findfile(ovfil,count=count)
      if count NE 0 then begin
          struct[idx].flg_ov = 1
          struct[idx].img_ov = ovfil
      endif else struct[idx].img_ov = ' '
      
;        Check for Final
      
      mskfil = strjoin(['Final/','f_',struct[idx].img_root,'*'])
      a = findfile(mskfil,count=count)
      if count NE 0 then begin
          struct[idx].flg_final = 1
          struct[idx].img_final = a[0]
      endif else struct[idx].img_final = ' '
        
;      Zero out files
      struct[idx].arc_fil = ' '
      struct[idx].arc_img = ' '
      struct[idx].flat_fil = ' '
      struct[idx].obj_id = -1L
      struct[idx].obj_fil = ' '
      struct[idx].img_ov = ' '

;;      Check gain for sequential Flats
;      if (struct[q].type EQ 'MFLT') then begin
;          if q GT 1 then begin
;              if struct[q-1].type EQ 'MFLT' then begin
;;                  if ((struct[q].side EQ struct[q-1].side) AND $
;                      (struct[q].colbin EQ struct[q-1].colbin) AND $
;                      (struct[q].rowbin EQ struct[q-1].rowbin)) then begin
;                      gain = hires_gain(data_prev, data)
;                      print, 'The gain measurement from the ratio is ', gain
;                  endif
;              endif  
;          endif
;          data_prev = data
;      endif
  endfor

  ;; Drop the first (null) entry
  struct = struct[1:*]

  ;; 
  struct.flg_anly = 1

; Close the ASCII file
  if not keyword_set( NOFILE ) then hires_wrstrct, struct

; Edit
  if keyword_set( EDIT ) then hires_editstrct, struct

; Write the structure to FITS
  if not keyword_set(NOFILE) then mwrfits, struct, outfil, /create

; All done
  if not keyword_set(SILENT) then $
    print, 'hires_strct: All done!  Fits file in ', outfil

end
