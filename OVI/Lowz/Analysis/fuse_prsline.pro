;+ 
; NAME:
; fuse_prsline
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   lowzovi_prsdat, stucture, filename
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
;   lowzovi_prsdat, struct, '/u/xavier/DLA/Lists/tot_dla.lst'
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   31-May-2001 Written by JXP
;   02-Jan-2003 Added metallicity sturcture
;-
;------------------------------------------------------------------------------
pro fuse_prsline, files, strct, OUTFIL=outfil

; lowzovi_prsdat -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lowzovi_prsdat, files, strct, OUTFIL= (v1.0)' 
    return
  endif 

;

;  if not keyword_set( list ) then begin
;      print, 'Using tot_ovi.lst as the list'
;      list = '/u/xavier/OVI/Lists/tot_ovi.lst'
;  endif

  fil=''
  dumc = ''
  dumr = 0.0
  dumd = double(0.0)
  dumda = dblarr(3)
  dumr2 = 0.0
  dumi = 0
  flg = 0
  instr = 0
  tmp = { fuselinstrct }

; Loop on files
  readcol, files, fil_nm, format='A'
  nfil = n_elements(fil_nm)

  close, 1
  for qq=0L,nfil-1 do begin

      ;; Num abs lin
      nlin = numlines(strtrim(fil_nm[qq],2))
      nlin = (nlin-2)/2
      substr = replicate(tmp, nlin)

      ;; open file
      openr, 1, strtrim(fil_nm[qq],2)
      readf, 1,  format='(a15)', dumc
      readf, 1,  format='(a15)', dumc

      ;;
      for i=0L,nlin-1 do begin
          ;; Dummy line
          readf, 1,  format='(a15)', dumc
          ;; 
          readf, 1,  format='(f10.5,f9.4,2x,f9.4,1x,i3,1x,i2)', $
            dumda, instr, flg
          ;; Fill it up
          substr[i].wrest = dumda[0]
          substr[i].wv_lim[0] = dumda[1]
          substr[i].wv_lim[1] = dumda[2]
          substr[i].flg = flg
          substr[i].instr = instr
          ;;
      endfor
          
      ;; close
      close, 1
      if qq EQ 0 then  strct = substr else strct = [strct, substr]
  endfor

  ;; ION
  ntot = n_elements(strct)
  getfnam, strct.wrest, fval, nam
  strct.ion = nam

  ;; Output
  if keyword_set( OUTFIL ) then mwrfits, strct, outfil, /create

  return
end
