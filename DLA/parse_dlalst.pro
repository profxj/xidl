;+ 
; NAME:
; parse_dlalst
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   parse_DLAlst, stucture, filename
;
; INPUTS:
;
; RETURNS:
;   structure      - IDL structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LIST - File
;  ION - Input ionic column densities
;  NOELM - Supress inputting Elemental values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_dlalst, struct, '/u/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;   02-Jan-2003 Added metallicity sturcture
;-
;------------------------------------------------------------------------------
pro parse_dlalst, supstrc, list, NOELM=noelm, ION=ion

; parse_dlalst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'parse_dlalst, struct, [filename], /NOELM/, /ION, /MTL (v1.2)' 
    return
  endif 

;


  if not keyword_set( list ) then begin
      print, 'Using tot_dla.lst as the list'
      list = '/u/xavier/DLA/Lists/tot_dla.lst'
  endif

  fil=''
  dumc = ''
  dumr = 0.0
  dumd = double(0.0)
  dumr2 = 0.0
  dumi = 0

;  strctcolm = {strctclm, flgclm: intarr(100),$
;	clm: fltarr(100),$
;	sigclm: fltarr(100) }


;  Structure for Elements
  dumelm = {strctelm, $ 
               flgclm: 0,$
               clm: 0.0,$
               sigclm: 0.0 }

;  Structure for Ions
  dumion  = {strction, $
             state: replicate(dumelm,8) $
            }

;  Structure for Metallicity
;  dummtl  = {strctmtl, $
;             flg: 0, $   ; 1 = non-refract; 2=refract + 0.3
;             val: 0., $  ; 
;             sig: 0. $
;            }
	
  strctnm = {sDLA, qso: '',$
             qso_ra: '',$
             qso_dec: '',$
             qso_zem: 0.0,$
             flg_QSOmag: 0,$
             qso_mag: 0.0,$
             zabs: 0.0,$
             NHI: 0.0,$
             sigNHI: fltarr(2),$
             abndfil: '',$
             flgFe: 0,$
             FeH: 0.0,$
             sigFeH: 0.0,$
             flgZn: 0,$
             ZnH: 0.0,$
             sigZnH: 0.0,$
             flgAlpha: 0,$
             Alpha: 0.0,$
             sigAlpha: 0.0,$
             flglw: 0,$
             lwfil: '',$
             lwwav: 0.0d,$
             lwvmn: 0.0,$
             lwvmx: 0.0,$
             lwfvel: 0.0,$
             lwfmm: 0.0,$
             lwfedg: 0.0,$
             lwftpk: 0.0,$
             flgCII: 0,$
             CII: 0.0,$
             sigCII: 0.0,$
             flgciv: 0.0,$
             civfil: '',$ 
             civwav: 0.0d,$
             civvmn: 0.0,$
             civvmx: 0.0,$
             civfvel: 0.0,$
             civfmm: 0.0,$
             civfedg: 0.0,$
             civftpk: 0.0,$
             civlwfdv: 0.0,$
             civlwfrto: 0.0,$
             civlwfnmm: 0.0,$
             civlwftvm: 0.0,$
             qso_ebv: 0.0,$
             ffilt: 0,$
             fslit: 0,$
             srvy: 0, $
             srvy_mag: 0., $
             ref: '', $
             flgmtl: 0, $
             mtl: 0., $
             sigmtl: 0., $
             ndfil: 0,$
             Hfil: '',$
             Efil: '',$
             Ufil: '',$
             Xfil: '',$
             elm: replicate(dumelm,90),$
             XH: replicate(dumelm,90), $
             ion: replicate(dumion,90) $
            }

; Parse the Base File
  readcol, list, listnms, format='A'
  ndla = n_elements(listnms)
  supstrc = replicate(strctnm,ndla)
;
  close, 1
  for i=0,ndla-1 do begin
      fil = listnms[i]
      openr, 1, fil
;		openr, 1, strmid(fil,0,xlc(fil))
      readf, 1,  format='(a15)', dumc
      supstrc[i].qso = dumc
      readf, 1,  format='(a15)', dumc
      supstrc[i].qso_ra = dumc
      readf, 1,  format='(a15)', dumc
      supstrc[i].qso_dec = dumc
      readf, 1,  format='(f9.6)', dumr
      supstrc[i].qso_zem = dumr
      readf, 1,  format='(i2)', dumi
      supstrc[i].flg_QSOmag = dumi
      readf, 1,  format='(f9.6)', dumr
      supstrc[i].qso_mag = dumr
      readf, 1,  format='(f9.6)', dumr
      supstrc[i].zabs = dumr
      readf, 1,  format='(f6.3)', dumr
      supstrc[i].NHI = dumr
      readf, 1,  format='(2f6.3)', dumr, dumr2
      supstrc[i].sigNHI[*] = [dumr,dumr2]
      readf, 1,  format='(a60)', dumc
      supstrc[i].abndfil = dumc
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
      supstrc[i].flglw = dumi
      readf, 1,  format='(a60)', dumc
      supstrc[i].lwfil = dumc
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
      readf, 1,  format='(a15)', dumc ; Ref
      supstrc[i].ref = dumc
      readf, 1,  format='(i4)', dumi ; metallicity
      supstrc[i].flgmtl = dumi   
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].mtl = dumr
      readf, 1,  format='(f7.3)', dumr
      supstrc[i].sigmtl = dumr
      close, 1
  endfor

  ;;  Read in Element Info
  if not keyword_set ( NOELM ) then fill_elmxh, supstrc
  ;;  Read in Ion info
  if keyword_set ( ION ) then fill_ion, supstrc
  
  return
end
