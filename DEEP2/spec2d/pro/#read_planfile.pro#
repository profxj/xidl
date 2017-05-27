pro read_planfile,  planfile, maskno, rawdatadir, outdatadir, $
        flatnames, arcnames, sciencenames, slitnums, nbright, ndim, nsimp, nlong, polyflag=polyflag,chips=chips, bluearc = bluearc, redarc = redarc, linelist = linelist, degree = degree
;+
; NAME:
;   read_planfile
;
; PURPOSE:
;   Reads in an ASCII file detailing plan on how to reduce a given
;   DEIMOS mask
;
; CALLING SEQUENCE:
;   read_planfile, planfile, maskno, rawdatadir, outdatadir, $
;          flatnames, arcnames, sciencenames $
;         [, slitnums, nbright, ndim, nsimp, nlong]
;   
; INPUTS:
;   planfile -- file name, (in current directory, or specify full path)
; 
; OUTPUTS:
;   maskno  -- string name of mask
;   rawdatadir  -- subdirectory (under DEIMOS_DATA) where raw data is stored
;   outdatadir  -- subdirectory (relative to current directory, or absolute) 
;                  for reduced data output. (if not specified or if
;                  '.' then output goes to current directory
;   flatnames -- string list of suitable flat files
;   arcnames  -- string list of arcs for this mask
;   sciencenames -- string list of science frames for this mask
;   slitnums -- integer list of slit numbers for this mask (if
;               quicklook)
;   nbright  -- int # of brightest slits for this mask (if quicklook)
;   ndim     -- int # of dimmest slits for this mask (if quicklook)
;   nsimp    -- int # of simplest slits for this mask (if quicklook)
;   nlong    -- int # of longest slits for this mask (if quicklook)
;   bluearc -- file names for blue-side arc
;   redarc -- file names for red-side arc
;   linelist -- file name for line list, in spec2d/etc
;   degree -- degree for arc line fit
;   
; KEYWORDS:
;   polyflag -- this passes back the value of polyflag specified in
;               the file.  'POLY' in the plan file would make
;               polyflag=1 (default), 'TRACE' polyflag=0 .
;   chips    -- this passes back a list of chips to reduce 
;  
;
; COMMENTS:
;   there is no restriction on the length of the file lists.  returns
;   full pathname for each file.  Plan file can read SUBDIRECTORY's  and
;   if included, it will insert into filename.  Be certain this
;   keyword preceeds the files that are to use it.  Multiple
;   subdirectories are fine.  Note that the code looks for 3 letter
;   keywords: FLAt, ARC, SCIence, MASk, RAWdirectory OUTdirectory.  
;   Lower or upper,
;   or mixed case is fine. File name must be on same line as keyword,
;   separated by tab or spaces, followed by more tabs or spaces
;   (comments OK on same line)
;
; REVISION HISTORY:
;   17may02 MD
;   05jun02 MD add OUTdirectory
;   14jun02 DPF add reality check at end
;----------------------------------------------------------------------


  deimos_data = getenv('DEIMOS_DATA')+'/'
;  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'
  if deimos_data eq '/' then deimos_data = '' ;not always set

  foo = findfile(planfile, count=filect)
  if filect eq 0 then begin 
     print, 'Cannot find plan file: ', planfile
     retall
  endif 
  
  openr, unit, planfile, /get_lun
 
  line =  ''
  flatnames = ''
  arcnames = ''
  sciencenames = ''
  rawdatadir = ''
  outdatadir = ''
  bluearc = ''
  redarc = ''
  linelist = ''
  slitnums = 0
  chips=''
  nbright = 0
  nlong = 0
  ndim = 0
  nsimp = 0
  polyflag = 1
  degree = 6

  while NOT EOF(unit) do begin

    readf, unit, line

;    substring = strsplit(line, /extract) ;split into two or more 
    substring = strsplit(strcompress(line), $
                         '[^A-Za-z0-9._#/]+', /REGEX, /extract)

    case strmid(strupcase(substring[0]), 0, 3) of ;chars 0-2, upper case
      'SCI': if sciencenames[0] eq '' then sciencenames = $
                rawdatadir + strcompress(substring[1], /remove_all) else $
                sciencenames= [sciencenames, rawdatadir + $
                  strcompress(substring[1], /rem)]

      'FLA': if flatnames[0] eq '' then flatnames = $
                rawdatadir + strcompress(substring[1], /remove_all) else $
                flatnames= [flatnames, rawdatadir + $
                    strcompress(substring[1], /rem)]

      'ARC': if arcnames[0] eq '' then arcnames = $
                rawdatadir + strcompress(substring[1], /remove_all) else $
                arcnames= [arcnames, rawdatadir + $
                  strcompress(substring[1], /remove_all)]

      'MAS': maskno =  strcompress(substring[1], /remove_all)

      'RAW': rawdatadir =  strcompress(substring[1], /remove_all) + '/'

      'OUT': outdatadir = strcompress(substring[1], /remove_all) + '/'

      'SLI': slitnums = fix(substring[1:*])

      'CHI': if chips eq '' then chips = strcompress(strjoin(substring[1:*],','), /remove_all) else chips=chips+','+strcompress(strjoin(substring[1:*],','), /remove_all)

      'BRI': nbright =  fix(strcompress(substring[1], /remove_all))

      'DIM': ndim =  fix(strcompress(substring[1], /remove_all))

      'LON': nlong =  fix(strcompress(substring[1], /remove_all))

      'SIM': nsimp =  fix(strcompress(substring[1], /remove_all))

      'POL': polyflag = 1

      'TRA': polyflag = 0

      'BLU': if bluearc[0] eq '' then bluearc = $
                rawdatadir + strcompress(substring[1], /remove_all) else $
                bluearc= [bluearc, rawdatadir + $
                  strcompress(substring[1], /remove_all)]

      'RED': if redarc[0] eq '' then redarc = $
                rawdatadir + strcompress(substring[1], /remove_all) else $
                redarc= [redarc, rawdatadir + $
                  strcompress(substring[1], /remove_all)]

      'LIN': linelist = strcompress(substring[1], /remove_all)

      'DEG': degree =  fix(strcompress(substring[1], /remove_all))


       else: a = 0 ;junk statement
      
    endcase

  endwhile

  quickfile = total(slitnums) gt 0 OR (nbright+ndim+nlong+nsimp) gt 0


  if chips eq "" then chips='1,2,3,4,5,6,7,8'

; adding capability to use different arcs for red/blue sides
  if arcnames[0] eq '' then  arcnames = bluearc
  if bluearc[0] eq '' then bluearc = arcnames
  if redarc[0] eq '' then redarc = arcnames
  if linelist eq '' then linelist = 'lamp_NIST.dat'
  
  if NOT quickfile then begin
     sciencenames = deimos_data  + sciencenames ;prepend path
     flatnames    = deimos_data  + flatnames 
     arcnames     = deimos_data  + arcnames
     bluearc      = deimos_data  + bluearc
     redarc       = deimos_data  + redarc
	
  foo = findfile(sciencenames[0], count=ct)
  if ct EQ 0 then message, 'Could not find first science image: ' $
      +sciencenames[0], /cont

  endif

; just do a reality check that these files exist!
  foo = findfile(arcnames[0], count=ct)
  if ct EQ 0 then message, 'Could not find first arc: '+arcnames[0], /cont

  foo = findfile(flatnames[0], count=ct)
  if ct EQ 0 then message, 'Could not find first flat: '+flatnames[0], /cont

  close, unit
  return
end 






