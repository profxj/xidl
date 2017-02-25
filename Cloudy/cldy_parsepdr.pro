;+
; NAME:
; cldy_parsepdr
;  V1.1
;
; PURPOSE:
;    Run a suite of cloudy models
; CALLING SEQUENCE:
;
;  cldy_parsepdr, infil, fluxfil, z, OUTDIR=,INDIR=,CLOUDYDIR=,NHI=,METALS=,U=,HDEN=
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

pro cldy_parsepdr, root, all_strct, OUTDIR=outdir, $
                   INDIR=indir, CLOUDYDIR=cloudydir, CIE=cie, $
                   NHI=nhi, METALS=metals, U=u, HDEN=hden, $
                   TITLE=title, FILNAME=filename, SHAPECMD=shapecmd, $
                   CLOBBER=clobber, NODOUB=nodoub, v0602a=v0602a

  if (N_params() LT 1) then begin
      print, 'Syntax - cldy_mkpdr, infil, template, root '
      return
  endif

  ;; Find all files
  all_files = findfile(root+'*', count=nfile)
  if nfile EQ 0 then return

  tmp2 = { $
         metal: 0., $
         hden: 0., $    ; Hydrogen Density
         crate: 0., $   ; CR rate
         scrate: 0., $  ; Secondary CR rate
         init_file: '', $
         grains_type: '', $
         chi: 0., $
         em_line_J: fltarr(100), $
         em_line_name: strarr(100), $
         em_line_wrest: dblarr(100) $
        }

  ;; Parse Input files
  in_files = findfile(root+'*.in', count=nin)
  nruns = nin
  all_strct = replicate(tmp2, nruns)
  
  for qq=0L,nruns-1 do begin
     close, 11
     nlines = numlines(in_files[qq])
     ;; 
     openr, 11, in_files[qq]
     for ii=0L,nlines-1 do begin
        dumc = ''
        readf, 11, dumc
        if strmid(dumc,0,4) EQ 'hden' then all_strct[qq].hden = float(strmid(dumc,4))
        if strmid(dumc,0,13) EQ 'metals grains' then all_strct[qq].metal = float(strmid(dumc,13))
        if strmid(dumc,0,15) EQ 'cosmic ray rate' then all_strct[qq].crate = float(strmid(dumc,15))
        if strmid(dumc,0,10) EQ 'set csupra' then all_strct[qq].scrate = float(strmid(dumc,10))
     endfor
  endfor
  close, 11
           

  ;; Parse CII* 158 from .lin files
  lin_files = findfile(root+'*.lin', count=nlinf)
  if nlinf NE nin then stop

  for qq=0L,nruns-1 do begin
     ;; Read the lines
     readcol, lin_files[qq], elm, ion, wrest, jnu, $
              format='A,L,D,D', /silent
     nout = n_elements(wrest)
     ;; Search for duplications
     a = where(abs(wrest[0]-wrest) LT 1e-3, ndup)
     nlins = nout / ndup

     ;; Fill in
     for ii=0L,nlins-1 do begin
        idx = (ndup-1)*nlins+ii
        all_strct[qq].em_line_J[ii] = jnu[idx]
        all_strct[qq].em_line_name[ii] = elm[idx]+'_'+strtrim(ion[idx],2)
        all_strct[qq].em_line_wrest[ii] = wrest[idx]
     endfor
  endfor

  ;; Write
  
  return

end
