;+ 
; NAME:
; x_parse_vpfit
;   Version 1.0
;
; PURPOSE:
;    Given a fort.26 file outputted by VPFIT, parse for the abundances
;     and return a structure with one struct per line.
;
; CALLING SEQUENCE:
;   
;   x_parse_vpfit, fil, struct, /SILENT
;
; INPUTS:
;   fil - VPFIT output file (fort.26)
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
;   06-Jul-2012  Modified from orig_prs_vpfit
;-
;------------------------------------------------------------------------------

pro x_parse_vpfit, fil, strct, SILENT=silent

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
      'x_parse_vpfit, fil, strct,  /silent [v1.0]'
    return
  endif 

  ;; Check file exists
  a = findfile(fil, count=nfil)
  if nfil EQ 0 then begin
      print, 'x_parse_vpfit: File ', fil, ' does not exist!'
      print, 'x_parse_vpfit: Please try again..'
      return
  endif

  nlin = numlines(fil)

  strct = {newabslinstrct}

  ;; Open file
  close, 94
  openr, 94, fil 

  lin = ''
  count = -1L

  svsetcz = ['xxxxx']
  svsetnz = [0L]
  svsetcb = ['xxxxx']
  svsetnb = [0L]

  ;; Loop
  for qq=0L,nlin-1 do begin
      ;; Read
      readf, 94, lin

      ;; Comment?
      if strmid(lin,0,1) NE ' ' then continue

      ;; Is it a special charachter?
      if strmatch(strmid(lin,1,1), '[_<?>-]') then continue
      
      ;; Should be an element
      count = count + 1

      ;; Figure out the Formatting
      prs_lin = strsplit(lin, ' ', /EXTRACT)
      nstr = n_elements(prs_lin)
      for ii=0L,nstr-1 do begin
         if  strpos(prs_lin[ii], '.') GE 0 then begin
            z_COLM = ii
            break
         endif
      endfor
      b_COLM = z_COLM + 2
      N_COLM = b_COLM + 2
      if COUNT EQ 0 then begin
         ;; Get the precision for z
         ipos = strpos(prs_lin[z_COLM+1], '.')
         nzsig = strlen(prs_lin[z_COLM+1]) - ipos - 1
         ;; Get the precision for b
         ipos = strpos(prs_lin[b_COLM+1], '.')
         nbsig = strlen(prs_lin[b_COLM+1]) - ipos - 1
         ;; Get the precision for N
         ipos = strpos(prs_lin[N_COLM+1], '.')
         nNsig = strlen(prs_lin[N_COLM+1]) - ipos - 1
      endif

      ;; Ion/Molecule
      ion = ''
      for ii=0L,z_COLM-1 do ion = ion+strtrim(prs_lin[ii],2)
      strct[count].ion = ion

      ;; Redshift
      ipos = strpos(prs_lin[z_COLM], '.')
      z = double( strmid(prs_lin[z_COLM], 0, ipos+nzsig+1))
      ipos = strpos(prs_lin[z_COLM+1], '.')
      zerr = double( strmid(prs_lin[z_COLM+1], 0, ipos+nzsig+1))
      strct[count].zabs = z
      strct[count].zsig = zerr 

      ;; Set for z
      if strlen(prs_lin[z_COLM]) GT (ipos+nzsig+1) then begin
         tag = strmid(prs_lin[z_COLM], ipos+nzsig+1)
         mtset = where(strmatch(svsetcz, tag, /FOLD_CASE), nmt)
         case nmt of
            0: begin
               svsetcz = [svsetcz, tag]
               svsetnz = [svsetnz, max(svsetnz)+1]
               strct[count].set += max(svsetnz)
            end
            1: strct[count].set += svsetnz[mtset]
            else: stop
         endcase
      endif

      ;; b value
      ipos = strpos(prs_lin[b_COLM], '.')
      b = double( strmid(prs_lin[b_COLM], 0, ipos+nbsig+1))
      ipos = strpos(prs_lin[b_COLM+1], '.')
      berr = double( strmid(prs_lin[b_COLM+1], 0, ipos+nbsig+1))
      strct[count].b = b
      strct[count].bsig = berr 

      ;; Set for b
      if strlen(prs_lin[b_COLM]) GT (ipos+nbsig+1) then begin
         tag = strmid(prs_lin[b_COLM], ipos+nbsig+1)
         mtset = where(strmatch(svsetcb, tag, /FOLD_CASE), nmt)
         case nmt of
            0: begin
               svsetcb = [svsetcb, tag]
               svsetnb = [svsetnb, max(svsetnb)+100]
               strct[count].set += max(svsetnb)
            end
            1: strct[count].set += svsetnb[mtset]
            else: stop
         endcase
      endif

      ;; N value
      ipos = strpos(prs_lin[N_COLM], '.')
      N = double( strmid(prs_lin[N_COLM], 0, ipos+nbsig+1))
      ipos = strpos(prs_lin[N_COLM+1], '.')
      Nerr = double( strmid(prs_lin[N_COLM+1], 0, ipos+nbsig+1))
      strct[count].N = N
      strct[count].Nsig = Nerr 

      ;; Increment
      if qq LT (nlin-1) then strct = [strct,{newabslinstrct}]

   endfor
  strct = strct[0:count-1]
  close, 94
  
  ;; Fill in tied component errors
  tied = where(strct.zsig LE 0., ntied)
  for jj=0L,ntied-1 do begin
      ;; Find the match
     mt = where( strct.set EQ strct[tied[jj]].set AND $
                 strct.zsig GT 0., nmt)
     case nmt of 
        0: begin
           if not keyword_set(SILENT) then $
              print, 'x_parse_vpfit: Warning, no redshift error for component ', tied[jj]
           ;stop
        end
        1: strct[tied[jj]].zsig = strct[mt[0]].zsig
        else: if total(strct[mt].zsig-median(strct[mt].zsig)) GT 1e-5 then stop
     endcase
  endfor
  
  ;; Repeat for b values
  tied = where(strct.bsig LE 0., ntied)
  for jj=0L,ntied-1 do begin
      ;; Find the match
     mt = where( strct.set EQ strct[tied[jj]].set AND $
                 strct.bsig GT 0., nmt)
     case nmt of 
        0: begin
           if not keyword_set(SILENT) then begin
              mtbs = where( svsetnb EQ ((strct[tied[jj]].set / 100) * 100) )
              print, 'x_parse_vpfit: Warning, no b error for component ', tied[jj], svsetcb[mtbs]
           endif
        end
        1: strct[tied[jj]].bsig = strct[mt[0]].bsig
        else: if total(strct[mt].bsig-median(strct[mt].bsig)) GT 1e-5 then stop
     endcase
  endfor
  

  return
end
  
     
  

