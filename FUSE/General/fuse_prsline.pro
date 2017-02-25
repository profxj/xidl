;+ 
; NAME:
; fuse_prsline
;  V1.1
;
; PURPOSE:
;    Create a structure for the absorption lines given a set
;  of files for the absorption lines.
;
; CALLING SEQUENCE:
;  fuse_prsline, files, strct, OUTFIL=
;
; INPUTS:
;   files -- List of files for absorption lines
;
; RETURNS:
;   strct  - FUSE absorption line structure
;
; OUTPUTS:
;  OUTFIL= -- Writes structure to OUTFIL
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   Nov-2003   Created by JXP
;-
;------------------------------------------------------------------------------
pro fuse_prsline, files, strct, OUTFIL=outfil, NOLIST=nolist, _extra=extra

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lowzovi_prsdat, files, strct, OUTFIL= [v1.1]' 
    return
  endif 

;
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
  if keyword_set(nolist) then fil_nm = files $
  else readcol, files, fil_nm, format='A'
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
          readf, 1,  format='(f10.5,f9.4,2x,f9.4,1x,i4,1x,i2)', $
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
  getfnam, strct.wrest, fval, nam, _extra=extra ; e.g. fil=getenv('XIDL_DIR')+'/Spec/Lines/Lists/lls.lst'
  strct.ion = nam

  ;; Trim
  gd = where(strtrim(strct.ion,2) ne '',ngd)
  if ngd eq 0 then stop,'fuse_prsline: no known lines'
  strct = strct[gd]

  ;; Output
  if keyword_set( OUTFIL ) then mwrfits, strct, outfil, /create

  return
end
