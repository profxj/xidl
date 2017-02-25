;+
; NAME:
; cldy_prsmodels
;  V1.1
;
; PURPOSE:
;    Run a sweet of cloudy models
; CALLING SEQUENCE:
;
;  cldy_prsmodels, infil, fluxfil, OUTDIR=,INDIR=,CLOUDYDIR=,NHI=,METALS=,U=,HDEN=
;
;
;  INPUT:
;    infil:
;       input file with:
;           <Variable name>  <desired minimum value> <desired maximum value> <step size>
;           variables which can be input  are:  NHI, nH, U, metals
;           NOTE:  MUST use these specific variable tags
;
;     NHI, METALS, U, HDEN:  specify a single value for these variables to be used
;                            in all models run (HDEN is for NH)
;
;     OUTDIR:  output file directory
;     INDIR:  input file directory
;     FILENAME: input/output file name lead (default: cloudy)
;
;  OUTPUT:
;     supstrc:  output cloudy structures
;
;
; REVISION HISTORY:
;   12-Sep-2005 Written by GEP
;-
;------------------------------------------------------------------------------

pro cldy_prsmodels, infil, supstrc, OUTDIR=outdir, INDIR=indir, $
                    CLOUDYDIR=cloudydir, NHI=nhi, $
                    METALS=metals, U=u, HDEN=hden, FILNAME=filename, $
                    NEND=nend, SILENT=silent

  if (N_params() LT 1) then begin
      print, 'Syntax - cldy_prsmodels, infil'
      return
  endif

  if not keyword_set(CLOUDYDIR) then cloudydir='~/LLS/CLOUDY/source/'

  if not keyword_set(infil) then begin
    print, 'infil not set!'
    stop
  endif

  if not keyword_set(FILENAME) then filename = 'cloudy'


  readcol, infil, varp, minp, maxp, stepp, FORMAT='A,F,F,F'

  nvars = n_elements(minp)

  if not keyword_set(nhi) then nhi = 0
  if not keyword_set(nh) then nh = 0
  if not keyword_set(metals) then metals = 0
  if not keyword_set(u) then u = 0

  if keyword_set(HDEN) then nh = hden

  for i=0L, nvars-1 do begin
      g = long((maxp[i]-minp[i])/stepp[i])+1
      case varp[i] of
          'NHI': nhi = findgen(g)*stepp[i]+minp[i]
          'nH': nh = findgen(g)*stepp[i]+minp[i]
          'metals': metals = findgen(g)*stepp[i]+minp[i]
          'U': u = findgen(g)*stepp[i]+minp[i]
          else: print, 'cldy_mkmodels:  not a valid variable!: '+varp[i]
      endcase
  endfor

  if not keyword_set(nhi) then nhi = 17.
  if not keyword_set(nh) then nh = -2
  if not keyword_set(metals) then metals = -1
  if not keyword_set(u) then u = -4

  nnhi = n_elements(nhi)
  nmetals = n_elements(metals)
  nu = n_elements(u)
  nnh = n_elements(nh)
  
  nloops = nnhi*nmetals*nu*nnh

  tmp1 = {strctcldy}
  supstrc = replicate(tmp1,nloops)
  step=0L

  if not keyword_set(NEND) then nend = nnhi*nmetals*nu*nnh

  for i=0,nnhi-1 do begin
      for j=0, nmetals-1 do begin
          for k=0, nu-1 do begin
              for l=0, nnh-1 do begin
                  outfil=filename+'_'+strtrim(string(i),2)+'_'+strtrim(string(j),2)+'_' $
                             +strtrim(string(k),2)+'_'+strtrim(string(l),2)+'.out'
                  if keyword_set(OUTDIR) then outfil = outdir+outfil
                  if not x_chkfil(outfil) then begin
                      if not keyword_set(SILENT) then print, 'cldy_prsmodels: file not found!  '+outfil
                      continue
                  endif
                  cldy_prsout, outfil, cldystrct=cldystrct
                  supstrc[step] = cldystrct
                  step = step+1
                  if step GE nend then break
              endfor
              if step GE nend then break
          endfor
          if step GE nend then break
      endfor
      if step GE nend then break
  endfor

end
