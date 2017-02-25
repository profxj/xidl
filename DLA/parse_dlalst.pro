;+ 
; NAME:
; parse_dlalst
;  V1.2
;
; PURPOSE:
;    Given a list of DLA base files, fill up the DLA structure.  THis
;    program also defines the DLA structure {sDLA}.
;
; CALLING SEQUENCE:
;   
;   parse_DLAlst, dla, [list]
;
; INPUTS:
;  [list] -- List of DLA base files.  [default:
;            '/u/xavier/DLA/Lists/tot_dla.lst']
;
; RETURNS:
;
; OUTPUTS:
;  dla  -- IDL DLA structure
;
; OPTIONAL KEYWORDS:
;  /ION - Input ionic column densities
;  /NOELM - Supress elemental and ion structures (saves memory)
;  /NORELM - Supress inputting Elemental [X/H] values
;  FILE=  -- The input is the dat filename not a list of dat files
;  /NOHIS -- Suppress HI error in [X/H] values
;  /EW    -- Fill up the ion arrays with EW values instead of column
;            densities.  This reads the .EW files instead of .ion
;  /FINE  -- This fills up arrays for fine-structure states
;            (e.g. SiII*).  By default CII* is already considered 
;  ROOT=  Path to the DLA tree
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_dlalst, dla, '/u/xavier/DLA/Lists/tot_dla.lst'
;   parse_dlalst, dla, '/u/xavier/DLA/Data/GRB051111.z154.dat', /FILE, /ION
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;   02-Jan-2003 Added metallicity structure
;-
;------------------------------------------------------------------------------
pro parse_dlalst, supstrc, list, NOELM=noelm, ION=ion, NOHIS=nohis, FILE=file,$
                  NORELM=norelm, ROOT=root, FINE=fine, EW=ew

; parse_dlalst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
          'parse_dlalst, struct, [list], /NOELM/, /ION, ' + $
          '/NOHIS, /FILE, ROOT=, /EW [v1.2]' 
    return
  endif 

;

  if not keyword_set(ROOT) then root = ''

  if not keyword_set( list ) and not keyword_set(FILE) then begin
     print, 'Using all_mtl.lst as the list'
     if strlen(getenv('DLA')) GT 0 then begin
        print, 'Using getnev(DLA)'
        list = getenv('DLA')+'/Lists/all_mtl.lst'
        root = getenv('DLA')
     endif else begin
        list = '~/DLA/Lists/all_mtl.lst'
        root = '~/DLA/'
     endelse
  endif

  fil=''
  dumc = ''
  dumr = 0.0d
  dumd = double(0.0)
  dumr2 = 0.0d
  dumi = 0
  dumi2 = 0
  dumi3 = 0L


  if not keyword_set( NOELM ) then dlas = {dlastruct} $
  else dlas = {dlanoestruct} 

; Parse the Base File
  if not keyword_set( FILE ) then begin
      readcol, list, listnms, format='A', /silent
  endif else begin
      listnms = list
  endelse
  ndla = n_elements(listnms)

  ;; Create Big struct
  supstrc = replicate(dlas,ndla)
;
  close, 1
  for i=0,ndla-1 do begin
      fil = ROOT+listnms[i]
      nlin = file_lines(fil)
      openr, 1, fil
;		openr, 1, strmid(fil,0,xlc(fil))
      supstrc[i].dlafil = strtrim(fil,2)
      readf, 1,  format='(a60)', dumc
      supstrc[i].qso = strtrim(dumc,2)
      readf, 1,  format='(a19)', dumc
      supstrc[i].qso_ra = strtrim(dumc,2)
      readf, 1,  format='(a19)', dumc
      supstrc[i].qso_dec = strtrim(dumc,2)
      readf, 1,  format='(f9.6)', dumr
      supstrc[i].qso_zem = dumr
      readf, 1,  format='(i2)', dumi
      supstrc[i].flg_QSOmag = dumi
      readf, 1,  format='(f9.6)', dumr
      supstrc[i].qso_mag = dumr
      readf, 1,  format='(f12.9)', dumr
      supstrc[i].zabs = dumr
      readf, 1,  format='(f6.3)', dumr
      supstrc[i].NHI = dumr
      readf, 1,  format='(2f6.3)', dumr, dumr2
      supstrc[i].sigNHI[*] = [dumr,dumr2]
      readf, 1,  format='(a60)', dumc
      ;; Dont put the ROOT here.  It will get added to the output!
      ;supstrc[i].abndfil = ROOT+strtrim(dumc,2)
      supstrc[i].abndfil = strtrim(dumc,2)
      readf, 1,  format='(i3)', dumi
      supstrc[i].flgFe = dumi   ; Fe
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].FeH = dumr
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].sigFeH = dumr
      readf, 1,  format='(i3)', dumi
      supstrc[i].flgZn = dumi
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].ZnH = dumr
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].sigZnH = dumr
      readf, 1,  format='(i3)', dumi ; Alpha
      supstrc[i].flgAlpha = dumi  
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].Alpha = dumr
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].sigAlpha = dumr
      readf, 1,  format='(i2)', dumi
      ;; Low ions
      supstrc[i].flglw = dumi
      readf, 1,  format='(a60)', dumc
      supstrc[i].lwfil = strtrim(dumc,2)
      readf, 1,  format='(f9.4)', dumd
      supstrc[i].lwwav = dumd
      readf, 1,  format='(2f7.1)', dumr, dumr2
      supstrc[i].lwvmn = dumr
      supstrc[i].lwvmx = dumr2
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].lwfvel = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].lwfmm = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].lwfedg = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].lwftpk = dumr
      readf, 1,  format='(i2)', dumi
      supstrc[i].flgCII = dumi
      readf, 1,  format='(f6.3)', dumr
      supstrc[i].CII = dumr
      readf, 1,  format='(f5.3)', dumr
      supstrc[i].sigCII = dumr
      readf, 1,  format='(i2)', dumi
      supstrc[i].flgciv = dumi
      readf, 1,  format='(a60)', dumc
      supstrc[i].civfil = dumc
      readf, 1,  format='(f9.4)', dumd
      supstrc[i].civwav = dumd
      readf, 1,  format='(2f7.1)', dumr, dumr2
      supstrc[i].civvmn = dumr
      supstrc[i].civvmx = dumr2
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civfvel = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civfmm = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civfedg = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civftpk = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civlwfdv = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civlwfrto = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civlwfnmm = dumr
      readf, 1,  format='(f7.2)', dumr
      supstrc[i].civlwftvm = dumr
      readf, 1,  format='(f5.2)', dumr
      supstrc[i].qso_ebv = dumr
      readf, 1,  format='(i2)', dumi
      supstrc[i].ffilt = dumi
      readf, 1,  format='(i2)', dumi
      supstrc[i].fslit = dumi
      readf, 1,  format='(i2)', dumi  ; QSO Survey
      supstrc[i].srvy = dumi
      readf, 1,  format='(f5.2)', dumr
      supstrc[i].srvy_mag = dumr
      readf, 1,  format='(a60)', dumc ; Ref
      supstrc[i].ref = strtrim(dumc,2)
      readf, 1,  format='(i6)', dumi ; metallicity
      supstrc[i].flgmtl = dumi   
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].mtl = dumr
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].sigmtl = dumr
      readf, 1,  format='(i5,1x,i5,1x,i6)', dumi, dumi2, dumi3 ; plate, fibid,MJD
      supstrc[i].sdss_plate = dumi
      supstrc[i].sdss_fibid = dumi2
      supstrc[i].sdss_mjd = dumi3
      ;; H2 and CI
      if nlin GT 52 then begin  ;; Should eventually make this automatic
         ;; VPFIT file
         readf, 1,  format='(a60)', dumc
         supstrc[i].vpfit_fil = strtrim(dumc,2)
         ;; CI
         readf, 1,  format='(i2)', dumi
         supstrc[i].flg_CI = dumi
         readf, 1,  format='(f6.3)', dumr
         supstrc[i].CI = dumr
         readf, 1,  format='(f5.3)', dumr
         supstrc[i].sig_CI = dumr
         ;; H2
         readf, 1,  format='(i2)', dumi
         supstrc[i].flg_H2 = dumi
         readf, 1,  format='(f6.3)', dumr
         supstrc[i].H2 = dumr
         readf, 1,  format='(f5.3)', dumr
         supstrc[i].sig_H2 = dumr
      endif  
      close, 1
  endfor

  ;;  Read in Element Info
  if not keyword_set ( NOELM ) AND not keyword_set(NORELM) then $
    fill_elmxh, supstrc, NOHIS=nohis, ROOT=root

  ;;  Read in Ion info
  if keyword_set ( ION ) then fill_ion, supstrc, ROOT=root

  ;;  Read in EW info
  if keyword_set ( EW ) then dla_fillew, supstrc, ROOT=root

  ;; Fine structure
  if keyword_set ( FINE ) then x_fineabnd, supstrc

  return
end
