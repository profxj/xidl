;+ 
; NAME:
; fuse_printion
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
;   fuse_calccolm, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   17-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_printion, strct_fil, zabs, OUTFIL=outfil, ADDHI=addhi

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'fuse_printion, strct, zabs, NHI= (v1.0)' 
    return
  endif 

;  Open structure
  if size(strct_fil, /type) EQ 7 then $
    strct = xmrdfits(strct_fil, 1, /silent) $
  else strct = strct_fil

  ;; Grab key lines
  a = where(abs(strct.zabs-zabs) LT 0.001 AND $
           strct.flg NE 0, nsys)

  strct = strct[a]

  ;;  Get Z, ion
  atmval = intarr(nsys)
  ionval = intarr(nsys)
  for i=0L,nsys-1 do begin
      getion, strct[i].wrest, ion, celm
      getabnd, celm, atm, abnd
      atmval[i] = atm
      ionval[i] = ion
  endfor

  ;; Parse out HI
;  if not keyword_set( ADDHI ) then a = where(atmval NE 1) $
;  else a = lindgen(nsys)

  ;; Print
;  printcol, strct[a].ion, strct[a].Ncolm, strct[a].sigNcolm, $
;    strct[a].flg

  ;; OUTFIL
  if keyword_set( OUTFIL ) then begin
      close, /all
      openw, 11, outfil
  endif

  if keyword_set(ADDHI) then sZ = 1L else sz = 2L
  for nZ=sZ,100 do begin
      for nion=1L,10L do begin
          a = where(atmval EQ nZ AND ionval EQ nion, ngd)
          if ngd EQ 0 then continue
          printcol, strct[a].ion, strct[a].Ncolm, strct[a].sigNcolm, $
            strct[a].flg
          ;; Outfil
          if keyword_set( OUTFIL ) then $
            writecol, strct[a].ion, strct[a].Ncolm, strct[a].sigNcolm, $
            strct[a].flg, FILNUM=11
      endfor
  endfor
      
  return
end
