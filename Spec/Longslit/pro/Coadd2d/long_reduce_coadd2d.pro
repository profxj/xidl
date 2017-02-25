PRO LONG_REDUCE_COADD2D


coadd_plan = findfile('plan*coadd2d.par', count = nplan_coadd)
IF nplan_coadd EQ 0 THEN BEGIN 
   if (NOT keyword_set(planfile)) then $
      planfile = findfile('plan*.par', count = nplan)
   nplan = n_elements(planfile)
   coadd_file = findfile('coadd_2d*.fits*', count = ncoadd)
   IF nplan GT 1 OR ncoadd GT 1 THEN $
      message, 'This script currently only supports a single plan.par file and a single coadd_2d*.fits file'
   lines = djs_readlines(planfile)
   nlines = n_elements(lines)
;; find all science files in the plan.par and remove them
   qtrim = lonarr(nlines) + 1L
   isci = WHERE(strmatch(lines, 'LEXP*science*'), nsci)
   sciline = lines[isci[0]]
   IF nsci EQ 0 THEN message, 'No science exposures in this plan file'
   basename = repstr(coadd_file, 'coadd_2d-', '')
   coadd_sciline = repstr(sciline, basename, coadd_file)
   qtrim[isci] = 0
   ikeep = where(qtrim)
   new_lines = lines[ikeep]
   new_lines = [new_lines, coadd_sciline]
   hdr_line = WHERE(strmatch(new_lines, '*indir*'), nhdr)
   IF nhdr EQ 0 THEN message, 'No indir found in this plan file'
   new_lines[hdr_line] = "indir './' # Raw data diretory"
   newplanfile = repstr(planfile, '.par', '-coadd2d.par')
   stop
   nlines = n_elements(new_lines)
   openw, 1, newplanfile
   FOR ii = 0L, nlines-1L DO printf, 1, new_lines[ii]
   close, 1
ENDIF
   ;;long_reduce, coadd_plan

END
