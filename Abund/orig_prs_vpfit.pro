;+ 
; NAME:
; prs_vpfit
;   Version 1.1
;
; PURPOSE:
;    Given a fort.18 file, parse for the abundances
;     and return a structure with one struct per line.
;
; CALLING SEQUENCE:
;   
;   prs_vpfit, fil, struct, /SILENT
;
; INPUTS:
;   fil - VPFIT output file (fort.18)
;
; RETURNS:
;
; OUTPUTS:
;   struct -- Structure containing the column densities and errors
;
; OPTIONAL KEYWORDS:
;   /SILENT -- Restrict output to the screen
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
;   06-Apr-2007 Written by JXP
;   06-Jul-2012 Discontinued by JXP  -- Use x_parse_vpfit
;-
;------------------------------------------------------------------------------

pro orig_prs_vpfit, fil, strct, SILENT=silent

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
      'prs_vpfit, fil, strct,  /silent [v1.0]'
    return
  endif 

  ;; Check file exists
  a = findfile(fil, count=nfil)
  if nfil EQ 0 then begin
      print, 'prs_vpfit: File ', fil, ' does not exist!'
      print, 'prs_vpfit: Please try again..'
      return
  endif

  nlin = numlines(fil)

  strct = {newabslinstrct}

  ;; Open file
  close, 94
  openr, 94, fil 

  lin = ''
  svzf = ['']
  svbf = ['']
  count = -1L
  ;; Loop
  for jj=0L,nlin-1 do begin
      ;; Read
      readf, 94, lin

      ;; Comment?
      if strmid(lin,0,1) EQ '%' then continue
      count = count + 1

      ;; Parse
      elmc = strtrim(strmid(lin,1,2),2)
      ionc = strtrim(strmid(lin,3,5),2)
      z = double(strmid(lin,9,8))
      zflg = strtrim(strmid(lin,17,2),2)
      zerr = double(strmid(lin,20,8))
      b = float(strmid(lin,31,5))
      bflg = strtrim(strmid(lin,36,2),2)
      berr = float(strmid(lin,41,5))
      N = float(strmid(lin,48,6))
      Nsig = float(strmid(lin,58,5))

      ;; Fill up
      strct[count].ion = elmc+ionc
      strct[count].zabs = z
      strct[count].zsig = zerr 
      strct[count].b = b
      strct[count].bsig = berr
      strct[count].N = N
      strct[count].Nsig = Nsig

      svzf = [svzf,zflg]
      svbf = [svbf,bflg]

      ;; Increment
      if jj LT (nlin-1) then strct = [strct,{newabslinstrct}]
  endfor
  svzf = svzf[1:*]
  svbf = svbf[1:*]

  ;; Fill in tied component errors
  tied = where(strct.zsig LE 0., ntied)
  for jj=0L,ntied-1 do begin
      ;; Find the match
      mt = where( strmatch(svzf, svzf[tied[jj]],/fold_case) AND $
                  strct.zsig GT 0., nmt)
      if nmt EQ 0 then begin
          if not keyword_set(SILENT) then $
            print, 'prs_vpfit: Warning, no redshift error for compoent ', $
                 svzf[tied[jj]]
      endif else strct[tied[jj]].zsig = strct[mt[0]].zsig
  endfor
  
  ;; Repeat for b values
  tied = where(strct.bsig LE 0., ntied)
  for jj=0L,ntied-1 do begin
      ;; Find the match
      mt = where( strmatch(svbf, svbf[tied[jj]],/fold_case) AND $
                  strct.bsig GT 0., nmt)
      if nmt EQ 0 then begin
          print, 'prs_vpfit: Warning, no Dopp error for compoent ', $
                 svbf[tied[jj]]
      endif else strct[tied[jj]].bsig = strct[mt[0]].bsig
  endfor
  close, 94
  
  return
end
  
     
  

