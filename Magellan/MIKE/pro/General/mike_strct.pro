;+ 
; NAME:
; mike_strct   
;     Version 2.0
;
; PURPOSE:
;    Creates and outputs a structure for a series of MIKE images.
;    This structure organizes the data for the night and is used 
;    to run most of the programs in the MIKE package.  It will
;    attempt to identify the type of object using the algorithm
;    (mike_autoid).  It also parses the header for information related
;    to exposure time, RA,DEC, etc. Finally, the code calculates the
;    gain for pairs of Milky flat exposures with the routine
;    mike_gain.  At present these gain values are not used.
;    Use mike_setgain to set the values [highly recommended].
;
;    See mikestrct__define.pro for all of the tags of the mike structure.
;
;
; CALLING SEQUENCE:
;   
;  mike_strct, struct, LIST=, /NOMKDIR, /NOFILE, OUTFIL=, /NOEDIT, FILE_LIST=
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
;   NOEDIT     - Do not edit the final structure
;   OUTFIL     - Name of fits output file (default = mikestrct.fits)
;   FILE_LIST  - String Array of all input filenames
;   NOFILE     - Do not write out the final structure (not
;                recommended)
;   
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   mike_strct, mike
;
; PROCEDURES/FUNCTIONS CALLED:
;    mike_rslvall
;    mike_autoid
;    mike_gain  (x_calcgain)
;    mike_wrstrct
;    mike_editstrct
;
; REVISION HISTORY:
;   27-Feb-2003 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro mike_strct, struct, LIST=list, NOMKDIR=nomkdir, NOFILE=nofile, $
                   OUTFIL=outfil, NOEDIT=noedit, FILE_LIST=file_list, $
                CHK_GAIN=chk_gain

;
  if  N_params() LT 1  then begin 
      print,'Syntax - ' + $
        'mike_strct, struct, LIST=, /NOMKDIR, /NOFILE, ' + $
             'OUTFIL=, /NOEDIT, FILE_LIST=file_list, /CHK_GAIN (v2.0)'
      return
  endif 
  
;  Optional Keywords
  
  if not keyword_set( TEL ) then tel = 'Magellan'
  if not keyword_set( CCD ) then ccd = 'Site#1'
  if not keyword_set( OUTFIL ) then outfil = 'mikestrct.fits'


; Resolve the programs
  mike_rslvall

  if NOT keyword_set( FILE_LIST ) then begin
     
; List
    if not keyword_set( LIST ) then begin
      mike_file = '[r,b][0-9][0-9][0-9][0-9].fits'

      img_old = findfile('Raw/m'+mike_file+'*') 
      img_new = findfile('Raw/'+mike_file+'*') 
      img = [img_old,img_new]
      notempty = where(img NE '',nimg)
      
      if nimg EQ 0 then begin
          print, 'mike_strct: No images in Raw!'
          return
      endif
      img=img[notempty]

    endif else begin
      if x_chkfil(list) NE 1 then begin
          print, 'mike_strct: Trouble with filename ', list
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
      a = findfile('Maps/..', count=count)
      if count EQ 0 then file_mkdir, 'Maps'
      a = findfile('OV/..', count=count)
      if count EQ 0 then file_mkdir, 'OV'
      a = findfile('Final/..', count=count)
      if count EQ 0 then file_mkdir, 'Final'
      a = findfile('Flats/..', count=count)
      if count EQ 0 then file_mkdir, 'Flats'
      a = findfile('Lists/..', count=count)
      if count EQ 0 then file_mkdir, 'Lists'
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

  case ccd of 
      'Site#1' : begin
          expcrd = 'EXPTIME'
          frmcrd = 'FRAMENO'
          utcrd = 'UT-TIME'
          ccdspeed = 'SPEED'
          racrd = 'RA-D'
          deccrd = 'DEC-D'
          eqxcrd = 'EQUINOX'
          slitnm = 'SLITSIZE'
          objcrd = 'OBJECT'
          binning = 'BINNING'
          amcrd = 'AIRMASS'
          datecrd = 'UT-DATE'
      end
      else : begin
          print, 'Not prepared for this ccd ', ccd
          return
      end
  endcase
          
  
;  Create the Structure
  tmp = { mikestrct }
  struct = replicate(tmp,nimg) 
  struct.ut = -1L
; Default all files to be analysed
  struct.flg_anly = 1

;  Loop on Indv images
  for q=0,nimg-1 do begin
      print, 'Reading ', img[q], format='(a,a,$)'
      data = 0
      head = ''
      guess = mike_autoid(data, xbinguess, ybinguess, $
                       filename=img[q], hdr=head, /silent)
;
;      Let's check start time of exposure, end time would be nice too
;       see if telescope is tracking or not...
;
      uttime = sxpar(head, 'UT-TIME')
      if size(uttime,/tname) EQ 'LONG' then begin
        uthour = uttime / 3600.
        if uthour GT 11.5 AND uthour LT 22.0 then begin
          if guess EQ 'OBJ' then begin
            guess = 'UNK'
            print, ' ... Exposure between sunrise and sunset, no OBJ allowed'
          endif
        endif
      endif


;; Set JD
      utd = sxpar(head, datecrd)
      utt = sxpar(head, 'UT-START')
      struct[q].date = x_setjdate(utd, utt)
      
      ll = strmid(guess,0,1)
      print, " ...it's ", (ll EQ 'A') OR (ll EQ 'O') ? " an " : " a  ",  guess
      
      ;; Frame+Side
      ;;  Changed to deal without 'm'
      name = sxpar(head,'FILENAME')
      instr = sxpar(head,'INSTRUME')
;      mhere = strpos(name, 'm') 
;      case strmid(name,mhere+1,1) of
      posb = strpos(name, 'b')
      posr = strpos(name, 'r')
      if posb NE -1 then begin
          if strmid(instr,5,3) NE 'Blu' then stop
          struct[q].side = 1
      endif else begin
          if strmid(instr,5,3) NE 'Red' then stop
          struct[q].side = 2
      endelse
      struct[q].frame = long(strmid(name,(posr>posb)+1))

      ;; Parse the Header
      if keyword_set( expcrd ) then struct[q].exp = sxpar(head, expcrd)
      if keyword_set( gaincrd ) then ccdgain = sxpar(head, gaincrd)
      if keyword_set( racrd ) then begin
           temp_ra = sxpar(head, racrd)
           if (size(temp_ra,/tname) EQ 'INT') then temp_ra = sxpar(head, 'RA')
           if (size(temp_ra,/tname) EQ 'INT') then begin 
               print, 'Cannot find RA Header Card, setting to UNK'
               temp_ra = 'UNK' 
           endif
           struct[q].RA = temp_ra
      endif

      if keyword_set( deccrd ) then begin
           temp_dec = sxpar(head, deccrd)
           if (size(temp_dec,/tname) EQ 'INT') then temp_dec = sxpar(head, 'DEC')
           if (size(temp_dec,/tname) EQ 'INT') then begin 
               print, 'Cannot find DEC Header Card, setting to UNK'
               temp_dec = 'UNK' 
           endif
           struct[q].DEC = temp_dec
      endif

      if keyword_set( eqxcrd ) then struct[q].equinox = sxpar(head, eqxcrd)
      if keyword_set( utcrd ) then struct[q].UT = sxpar(head, utcrd)
      if keyword_set( objcrd ) then struct[q].Obj = sxpar(head, objcrd)
      if keyword_set( amcrd ) then struct[q].AM = sxpar(head, amcrd)

; CCD
      ;; Gain, readnoise
      struct[q].ccd = ccd
      struct[q].tel = tel
      struct[q].ccdspeed = sxpar(head, ccdspeed)
      case strmid(struct[q].ccdspeed,0,4) of
          'Slow': begin
              struct[q].readno = 3.70
           if struct[q].side EQ 2 then struct[q].gain = 0.93  $
           else begin
              if (struct[q].date LT 2453132.) then struct[q].gain = 1.01 $
              else begin
                  struct[q].gain = 0.47
                  struct[q].readno = 2.5
             endelse
           endelse
          end
          else: begin
             print, 'Fast readout? setting FLG_ANLY to 0'
             struct[q].flg_anly = 0
             struct[q].readno = 5.0  ; this is a guess
             if struct[q].side EQ 2 then struct[q].gain = 0.93 $
             else if (struct[q].date LT 2453132.) then struct[q].gain = 1.01 $
             else struct[q].gain = 0.47
          end
            
      endcase
      ;; Binning
      binstring = sxpar(head,binning)
      case binstring of
          '2x1': begin
              struct[q].colbin = 2
              struct[q].rowbin = 1
          end
          '2x2': begin
              struct[q].colbin = 2
              struct[q].rowbin = 2
          end
          '3x2': begin
              struct[q].colbin = 3
              struct[q].rowbin = 2
          end
          '3x3': begin
              struct[q].colbin = 3
              struct[q].rowbin = 3
          end
          '4x2': begin
              struct[q].colbin = 4
              struct[q].rowbin = 2
          end
          else: begin
;            print, 'Unusual binning, but might as well give it a shot'
            struct[q].colbin = long(strmid(binstring,0,1))
            struct[q].rowbin = long(strmid(binstring,strpos(binstring,'x')+1,1))
            print, 'Grabbing binning from header.. ', $
               strcompress(string(struct[q].colbin, ' x', struct[q].rowbin)) 
               

          end
      endcase
      if struct[q].colbin EQ 0 then continue
          
      ;; SLIT
      if keyword_set( slitnm ) then begin
          slit = sxpar(head, slitnm)
          case slit of 
              else:
          endcase
      endif

;     IMAGE TYPE
      ;; Bias
      if struct[q].exp LT 1. then begin
          struct[q].type = 'ZRO'
      endif else struct[q].type = guess
      
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

      ;; Strip off the 'm'
;      if strmid(struct[q].img_root,0,1) EQ 'm' then $
;        struct[q].img_root = strmid(struct[q].img_root,1)
      
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
      struct[q].arc_img = ' '
      struct[q].map_fil = ' '
      struct[q].flat_fil = ' '
      struct[q].obj_id = -1L
      struct[q].obj_fil = ' '
      struct[q].img_ov = ' '


;;      Check gain for sequential Milky Flats

      if keyword_set(CHK_GAIN) then begin
          if (struct[q].type EQ 'MFLT') then begin
              if q GT 1 then begin
                  if struct[q-1].type EQ 'MFLT' then begin
                      if ((struct[q].side EQ struct[q-1].side) AND $
                          (struct[q].colbin EQ struct[q-1].colbin) AND $
                          (struct[q].rowbin EQ struct[q-1].rowbin)) then begin
                          gain = x_calcgain(data_prev, data)
                          print, 'The gain measurement from the ratio is ', gain
                      endif
                  endif  
              endif
              data_prev = data
          endif
      endif
  endfor

; Close the ASCII file
  if not keyword_set( NOFILE ) then mike_wrstrct, struct

; Edit
  if not keyword_set( NOEDIT ) then mike_editstrct, struct

; Write the structure to FITS
  mwrfits, struct, outfil, /create

; All done
  print, 'mike_strct: All done!  Fits file in ', outfil

end
