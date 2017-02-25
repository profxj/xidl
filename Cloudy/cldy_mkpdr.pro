;+
; NAME:
; cldy_mkpdr
;  V1.1
;
; PURPOSE:
;    Run a suite of cloudy models
; CALLING SEQUENCE:
;
;  cldy_mkmodels, infil, fluxfil, z, OUTDIR=,INDIR=,CLOUDYDIR=,NHI=,METALS=,U=,HDEN=
;                                    /CLOBBER
;
;  INPUT:
;    infil:  
;       input file with:
;           <Variable name>  <desired minimum value> <desired maximum value> <step size>
;
;
;           variables which can be input  are:  NHI, nH, U, metals
;           NOTE:  MUST use these specific variable tags
;
;     fluxfil: cldy_cuba output file defining flux
;     
;     z:  redshift 
;
;     NHI, METALS, U, HDEN:  specify a single value for these variables to be used 
;                            in all models run (HDEN is for NH)
;
;     OUTDIR:  output file directory
;     INDIR:  input file directory
;     TITLE:  cloudy title (default: model)
;     FILENAME: input/output file name lead (default: cloudy)
;     /CLOBBER: clobber previous runs with same FILENAME
;
;  NHI 16.0 18.0 0.25
;  metals -1.0 0.0 1.0
;  U -6.0 0.0 0.25
;  nH -2.0 0.0 2.0
;
;
; REVISION HISTORY:
;   10-Sep-2005 Written by GEP
;-
;------------------------------------------------------------------------------

pro cldy_mkpdr, infil, template, root, OUTDIR=outdir, $
                   INDIR=indir, CLOUDYDIR=cloudydir, CIE=cie, $
                   NHI=nhi, METALS=metals, U=u, HDEN=hden, $
                   TITLE=title, FILNAME=filename, SHAPECMD=shapecmd, $
                   CLOBBER=clobber, NODOUB=nodoub, v0602a=v0602a

  if (N_params() LT 2) then begin
      print, 'Syntax - cldy_mkpdr, infil, template, root '
      return
  endif

  if not keyword_set(CLOUDYDIR) then cloudydir=getenv('CLDYDIR')
  if not keyword_set(OUTDIR) then outdir = ''
  if not keyword_set(INDIR) then indir = ''

  if not keyword_set(infil) then begin
    print, 'infil not set!'
    stop
 endif
  
  ;; Read template
  close, /all
  nlin = numlines(template)
  templ_lines = strarr(nlin)
  openr, 12, template
  tmpc = ''
  for ii=0L,nlin-1 do begin
     readf, 12, tmpc, format='(a)'
     templ_lines[ii] = tmpc
  endfor
  close, /all

  if not keyword_set(FILENAME) then filename = 'cloudy'

  ;; Read in file
  readcol, infil, varp, minp, maxp, stepp, FORMAT='A,D,D,D'

  nvars = n_elements(minp)
  sz_vars = lonarr(nvars)
  szc_vars = lonarr(nvars)
  idx = lonarr(nvars)

  ;; Initialize
  for i=0L, nvars-1 do begin
     g = long((maxp[i]-minp[i]+1d-7)/stepp[i])+1
     sz_vars[i] = g
     case varp[i] of
        'cosmics': crate = findgen(g)*stepp[i]+minp[i]
        'secondary': sec_crate = findgen(g)*stepp[i]+minp[i]
        'density': density = findgen(g)*stepp[i]+minp[i]
        'metals': metals = findgen(g)*stepp[i]+minp[i]
        'chi': begin  ;; Exponent approach
           g = round( (alog(maxp[i]) - alog(minp[i]))/alog(stepp[i]) ) + 1
           chi = minp[i] * stepp[i]^findgen(g)
        end
        else: begin
           print, 'cldy_mkmodels:  not a valid variable!: '+varp[i]
           stop
        end
     endcase
  endfor

  ;; Total number of models
  n_total = 1
  for i=0L,nvars-1 do n_total = n_total * sz_vars[i]

  for ss=0L,n_total-1 do begin
     ;; Update the index
     flg = -1
     for jj=0L,nvars-1 do begin
        if idx[jj] LT (sz_vars[jj]-1) then begin
           idx[jj] = idx[jj] + 1
           flg = jj
        endif
        if FLG GE 0 then break
     endfor
  
     ;; Names
     dumc = ''
     for jj=0L,nvars-1 do begin
        dumc = dumc+strtrim(idx[jj],2)
        if jj LT (nvars-1) then dumc = dumc+'_'
     endfor
     infil = root+dumc+'.in'
     outroot = root+dumc
     outfil = root+dumc+'.out'

     ;; Modify Output file
     input_lines = templ_lines

     ;; Loop on input variables

     ;; Find the match
     case varp[flg] of
        'density': begin
           str_mt = 'hden'
           len_mt = 4L
           param = density
        end
        'cosmics': begin
           str_mt = 'cosmic ray rate'
           len_mt = 15L
           param = crate
        end
        'secondary': begin
           str_mt = 'set csupra'
           len_mt = 10L
           param = sec_crate
        end
        'metals': begin
           str_mt = 'metals grains'
           len_mt = 13L
           param = metals
        end
        'chi': begin
           str_mt = 'table draine'
           len_mt = 12L
           param = chi
        end
        else: stop
     endcase
     mtc = where(strmatch(strmid(input_lines,0,len_mt), str_mt), nmt)
     if nmt NE 1 then stop

     ;; 
     input_lines[mtc] = str_mt+string(param[idx[flg]])


     ;; Loop on output files
     for kk=0L,nlin-1 do begin
        ipos = strpos(input_lines[kk], 'outfile')
        if ipos GE 0 then begin
           input_lines[kk] = strmid(input_lines[kk],0,ipos)+outroot+$
                             strmid(input_lines[kk],ipos+7)
        endif
     endfor
                  
     ;; Generate Input file
     close, 13
     openw, 13, infil
     for kk=0L,nlin-1 do begin
        printf, 13, input_lines[kk]
     endfor
     close, 13

     ;; Run Cloudy
     spawn, CLOUDYDIR+'/cloudy.exe '+'<'+infil+'>'+' '+outfil
     if keyword_set(OUTDIR) then $
        spawn, 'mv '+outfil+' '+outdir
     if keyword_set(INDIR) then $
        spawn, 'mv '+infil+' '+indir
  endfor

end
