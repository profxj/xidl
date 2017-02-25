;+ 
; NAME:
; parse_ismlst
;  V1.2
;
; PURPOSE:
;    Given a list of ISM base files, fill up the ISM structure.  THis
;    program also defines the DLA structure {sDLA}.  Similar to
;    parse_dlalst.pro
;
; CALLING SEQUENCE:
;   
;   parse_DLAlst, ism, [list]
;
; INPUTS:
;  [list] -- List of DLA base files.  [default:
;            '/u/xavier/DLA/Lists/tot_ism.lst']
;
; RETURNS:
;
; OUTPUTS:
;  ism  -- ISM structure
;
; OPTIONAL KEYWORDS:
;  /ION - Input ionic column densities
;  /NOELM - Supress inputting Elemental [X/H] values
;  /NOHIS -- Suppress HI error in [X/H] values
;  FLAG=  -- 0=Default; 1=GRB
;  /FINE  -- This fills up arrays for fine-structure states
;            (e.g. SiII*).  By default CII* is already considered 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   parse_ismlst, ism, '/u/xavier/DLA/Lists/tot_ism.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;   02-Jan-2003 Added metallicity structure
;-
;------------------------------------------------------------------------------
pro parse_ismlst, supstrc, list, NOELM=noelm, ION=ion, NOHIS=nohis, FILE=file,$
                  FLAG=flag, FINE=fine

; parse_ismlst -- Reads in DLA data to a structure

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'parse_ismlst, struct, [list], /FINE, /NOELM/, /ION, /NOHIS, FLAG= (v1.2)' 
    return
  endif 

  if not keyword_set( FLAG ) then flag = 0

  ;; 
  fil=''
  dumc = ''
  dumr = 0.0d
  dumd = double(0.0)
  dumr2 = 0.0d
  dumi = 0


  ;; Define structure
  case flag of
      1: isms = {grbastruct}
      else: isms = {ismstruct}
  endcase

  ;; Parse the Base File
  if not keyword_set( FILE ) then begin
      readcol, list, listnms, format='A', /silent
      nism = n_elements(listnms)
  endif else begin
      listnms = list
      nism = 1L
  endelse

  ;; Create Big struct
  supstrc = replicate(isms,nism)
;
  close, 13
  for i=0,nism-1 do begin
      fil = listnms[i]
      openr, 13, fil
;		openr, 1, strmid(fil,0,xlc(fil))
      supstrc[i].ismfil = strtrim(fil,2)
      readf, 13,  format='(a60)', dumc
      supstrc[i].target = dumc
      readf, 13,  format='(a60)', dumc
      supstrc[i].targ_ra = dumc
      readf, 13,  format='(a60)', dumc
      supstrc[i].targ_dec = dumc
      case flag of
          1:
          else: begin
              readf, 13,  format='(f11.6)', dumr
              supstrc[i].targ_vel = dumr
          end
      endcase
      readf, 13,  format='(i2)', dumi
      supstrc[i].flg_mag = dumi
      readf, 13,  format='(f9.6)', dumr
      supstrc[i].targ_mag = dumr
      readf, 13,  format='(f12.9)', dumr
      supstrc[i].zabs = dumr
      readf, 13,  format='(f6.3)', dumr
      supstrc[i].NHI = dumr
      readf, 13,  format='(2f6.3)', dumr, dumr2
      supstrc[i].sigNHI[*] = [dumr,dumr2]
      readf, 13,  format='(a60)', dumc
      supstrc[i].abndfil = dumc
      ;; D/H
      case flag of
          1:
          else: begin
              readf, 13,  format='(i3)', dumi
              supstrc[i].flgD = dumi ; D
              readf, 13,  format='(f7.3)', dumr
              supstrc[i].DH = dumr
              readf, 13,  format='(f7.3)', dumr
              supstrc[i].sigDH = dumr
          end
      endcase
      ;; Fe/H
      readf, 13,  format='(i3)', dumi
      supstrc[i].flgFe = dumi   
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].FeH = dumr
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].sigFeH = dumr
      readf, 13,  format='(i3)', dumi
      supstrc[i].flgZn = dumi
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].ZnH = dumr
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].sigZnH = dumr
      ;; Alpha
      readf, 13,  format='(i3)', dumi 
      supstrc[i].flgAlpha = dumi  
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].Alpha = dumr
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].sigAlpha = dumr
      readf, 13,  format='(i2)', dumi
      ;; Low ions
      supstrc[i].flglw = dumi
      readf, 13,  format='(a60)', dumc
      supstrc[i].lwfil = dumc
      readf, 13,  format='(f9.4)', dumd
      supstrc[i].lwwav = dumd
      readf, 13,  format='(2f7.1)', dumr, dumr2
      supstrc[i].lwvmn = dumr
      supstrc[i].lwvmx = dumr2
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].lwfvel = dumr
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].lwfmm = dumr
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].lwfedg = dumr
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].lwftpk = dumr
      readf, 13,  format='(i2)', dumi
      ;; CII*
      supstrc[i].flgCII = dumi
      readf, 13,  format='(f6.3)', dumr
      supstrc[i].CII = dumr
      readf, 13,  format='(f5.3)', dumr
      supstrc[i].sigCII = dumr
      readf, 13,  format='(i2)', dumi
      ;; CIV
      supstrc[i].flgciv = dumi
      readf, 13,  format='(a60)', dumc
      supstrc[i].civfil = dumc
      readf, 13,  format='(f9.4)', dumd
      supstrc[i].civwav = dumd
      readf, 13,  format='(2f7.1)', dumr, dumr2
      supstrc[i].civvmn = dumr
      supstrc[i].civvmx = dumr2
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].civfvel = dumr
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].civfmm = dumr
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].civfedg = dumr
      readf, 13,  format='(f7.2)', dumr
      supstrc[i].civftpk = dumr
      ;; CIV cross
      case flag of
          1:
          else: begin
              readf, 13,  format='(f7.2)', dumr
              supstrc[i].civlwfdv = dumr
              readf, 13,  format='(f7.2)', dumr
              supstrc[i].civlwfrto = dumr
              readf, 13,  format='(f7.2)', dumr
              supstrc[i].civlwfnmm = dumr
              readf, 13,  format='(f7.2)', dumr
              supstrc[i].civlwftvm = dumr
          end
      endcase
      readf, 13,  format='(f5.2)', dumr
      supstrc[i].ebv = dumr
      case flag of
          1:
          9: begin
              readf, 13,  format='(i2)', dumi
              supstrc[i].ffilt = dumi
              readf, 13,  format='(i2)', dumi
              supstrc[i].fslit = dumi
              readf, 13,  format='(i2)', dumi ; QSO Survey
              supstrc[i].srvy = dumi
              readf, 13,  format='(f5.2)', dumr
              supstrc[i].srvy_mag = dumr
          end
          else:
      endcase
      ;; References
      readf, 13,  format='(a60)', dumc 
      supstrc[i].ref = dumc
      readf, 13,  format='(i4)', dumi ; metallicity
      supstrc[i].flgmtl = dumi   
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].mtl = dumr
      readf, 13,  format='(f7.3)', dumr
      supstrc[i].sigmtl = dumr
      close, 13
  endfor

  ;;  Read in Element Info
  if not keyword_set ( NOELM ) then fill_elmxh, supstrc, NOHIS=nohis
  ;;  Read in Ion info
  if keyword_set ( ION ) then fill_ion, supstrc
  ;; Fine structure
  if keyword_set ( FINE ) then x_fineabnd, supstrc
  
  return
end
