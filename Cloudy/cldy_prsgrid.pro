;+
; NAME:
; cldy_prsgrid
;  V1.1
;
; PURPOSE:
;    parse a set of CLOUDY output files and create a FITS file
; CALLING SEQUENCE:
;
;   cldy_prsout, outfil, cldystrct=cldystrct
;
; INPUT:
;    indir:  Directory containing CLOUDY output files
;    outfil: name of CLOUDY output file
;
; OUTPUT:
;    cldystrct: output cloudy strct
;
; REVISION HISTORY:
;   12-Sep-2005 Written by JXP
;-
;----------------------------------------------------------------------

pro cldy_prsgrid, indir, outfil, cldystrct=cldystrct, zabs=zabs, SPEC=spec, $
                  ROOT=root, OUTPUT=cldygrid, NOWRITE=nowrite

  if (N_params() LT 2) then begin
      print, 'Syntax - cldy_prsgrid, infil, outfil, /nowrite, OUTPUT=, SPEC=, ROOT= [v1.1]'
      return
  endif

  if not keyword_set(ROOT) then root = '*'

  ;; Read files
  fil = file_search(indir+'/'+ROOT+'*out*',count=nfil)
  ;fil = file_search(indir+'/'+ROOT+'28_10_*_4.out',count=nfil)
  if nfil EQ 0 then begin
      print, 'cldy_prsgrid: No files!!'
      stop
      return
  endif

  tmp={strctcldy}
  cldygrid = replicate(tmp,nfil)

  for qq=0L,nfil-1 do begin
      ;; Parse
     ;if strmatch(fil[qq], 'Output/cloudy_28_10_45_4.out') then stop
      cldy_prsout, fil[qq], SPEC=spec, cldystrct=cstr, zabs=zabs
      ;; Fill up
      cldygrid[qq] = cstr
  endfor
;  a = where(cldygrid.nh EQ -2.d and cldygrid.feh eq -1.5d and $
;            cldygrid.nhi eq 19.0d and  cldygrid.U eq -6.d)
;  stop

  ;; Sort ?

  ;; Write + compress
  if not keyword_set(NOWRITE) then begin
     mwrfits, cldygrid, outfil, /create
     spawn, 'gzip -f '+outfil
  endif

  ;; All done
  print, 'cldy_prsgrid: All done'

    return
end
