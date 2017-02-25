;+
; NAME:
;   xgmos_mkclscripts
;
; PURPOSE:
;   Reading in a plan file, this code then creates cl scripts
;    for the XGMOS reductions 
;
; CALLING SEQUENCE:
;   xgmos_mkclscripts, planfil
;
; INPUTS:
;  planfil -- Plan file.  Can be created with gmos_plan
;
; OPTIONAL INPUTS:
;   fileexpr   - File names in the input directory; default to '*.fits*'
;   indir      - Input directory(s) for reading files;
;                default to current directory
;   planfile   - Output plan file; default to 'plan.par'.
;                This file is put in the same directory as the raw data files.
;
; OUTPUT:
;
; COMMENTS:
;   One plan file is made for each input directory.
;
;   The following flavors of images are listed:
;     bias
;     domeflat
;     iflat (internal flat)
;     twiflat
;     arc
;     science
;
; EXAMPLES:
; gmos_plan,'*.fits','/b/martell/data_arx/09072005/Raw/',planfile='plan-master.par'
; BUGS:
;
; PROCEDURES CALLED:
;
; INTERNAL SUPPORT ROUTINES:
;
;
; REVISION HISTORY:
;   15-Jun-2009  Written by JXP
;   13-Apr-2011  GW: change of default flatfield fitting order (now a string)
;-
;------------------------------------------------------------------------------
pro xgmos_mkclscripts, planfile, BIASFIL=biasfil, FLATROOT=flatroot, $
                       SKIP_BIAS=skip_bias

   if not keyword_set(BIASFIL) then biasfil = 'XGMOS_BIAS'
   if not keyword_set(FLATROOT) then flatroot = 'XGMOSFLAT_'

   ;; Read
   planstr = yanny_readone(planfile, hdr=planhdr, /anonymous)
   
   nplan = n_elements(planstr)

   close, /all
   openw, 10, 'xgmos_master.cl'
   cd, './', curr=pwd

   ;; Create Bias script
   bias = where(strmatch(planstr.flavor,'bias'),nbias)
   if nbias NE 0 then begin
       openw, 2, 'bias_files.lst'
       for kk=0L,nbias-1 do begin
           pos = strpos(planstr[bias[kk]].filename, '.fits')
           printf, 2, strmid(planstr[bias[kk]].filename,0,pos)
       endfor
       close, 2
       openw, 1, 'xgmos_bias.cl'
       printf, 1, 'xgmosBias("@bias_files.lst", "'+pwd+'/",'+$
               '"'+biasfil+'")'
       close, 1
       ;; Master
       printf, 10, 'cl < xgmos_bias.cl'
   endif else begin
      if not keyword_set(SKIP_BIAS) then begin
         print, 'xgmos_mkclscripts:  You need to get some biases and try again!'
         return
      endif else print, 'xgmos_mkclscripts:  Skipping bias creation!'
   endelse
   
   
   ;; Create instrument list 
   indx = where(strmatch(planstr.flavor,'science') OR $
                strmatch(planstr.flavor, 'std'), nid)
   
   ;; Grating list
   mask_list = planstr[indx].grating+ replicate('_',nid)+$
               strtrim(long(planstr[indx].wave),2) 
   mask_list = mask_list[uniq(mask_list, sort(mask_list))]
   nsetup = n_elements(mask_list)
   
   ;; cl script
   openw, 11, 'xgmos_flats.cl'
   openw, 12, 'xgmos_proc.cl'
   
   ;; Loop on Setups
   for ii=0L,nsetup-1 do begin
       
       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; Flats
       flats = where(strmatch(planstr.flavor,'domeflat') AND $
                     strmatch(planstr.grating+ replicate('_',nplan)+$
                              strtrim(long(planstr.wave),2), $
                              mask_list[ii]), nflat)
       if nflat EQ 0 then begin
           print, 'xgmos_mkclscripts: No flats for ', mask_list[ii]
           print, 'xgmos_mkclscripts: Try again..'
           close, /all
           return
       endif
       ;; Lists
       filnm= 'flats_'+mask_list[ii]+'.lst'
       printf, 11, 'xgmosFlat("@'+filnm+'", "'+pwd+'/",'+$
               '"'+flatroot+mask_list[ii]+'","'+biasfil+'","15,15,25")'
       openw, 2, filnm
       for kk=0L,nflat-1 do begin
           pos = strpos(planstr[flats[kk]].filename, '.fits')
           printf, 2, strmid(planstr[flats[kk]].filename,0,pos)
           ;; Delete gs files
           printf, 11, 'imdel("gs'+strmid(planstr[flats[kk]].filename,0,pos)+$
                   '",verif-)'
       endfor
       close, 2

       ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; Arc
       arcs = where(strmatch(planstr.flavor,'arc') AND $
                    strmatch(planstr.grating+ replicate('_',nplan)+$
                             strtrim(long(planstr.wave),2), $
                             mask_list[ii]), narc)
       if narc EQ 0 then begin
           print, 'xgmos_mkclscripts: No arcs for ', mask_list[ii]
           print, 'xgmos_mkclscripts: Try again..'
           close, /all
           return
       endif

       ;;;;;;;;;;;;;;;;;;;;;;;;;;
       ;; Process script
       allproc = where((strmatch(planstr.flavor,'domeflat') OR $
                        strmatch(planstr.flavor, 'science') OR $
                        strmatch(planstr.flavor, 'std')) AND $
                       strmid(planstr.filename,0,1) NE 'X' AND $
                       strmatch(planstr.grating+ replicate('_', nplan)+$
                                strtrim(long(planstr.wave), 2), $
                                mask_list[ii]), nflat)
       arcs = where((strmatch(planstr.flavor, 'filtarc') OR $
                     strmatch(planstr.flavor, 'arc')) AND $
                    strmatch(planstr.grating+ replicate('_', nplan)+$
                             strtrim(long(planstr.wave), 2), $
                             mask_list[ii]), narc)
       ;; Lists
       if narc EQ 0 then begin
          print, 'xgmos_mkclscripts: No arcs for ', mask_list[ii]
          print, 'xgmos_mkclscripts: Try again..'
          close, /all
          return
       endif
       filnm= 'arcs_'+mask_list[ii]+'.lst'
       openw, 2, filnm
       for kk=0L,narc-1 do begin
           pos = strpos(planstr[arcs[kk]].filename, '.fits')
           printf, 2, strmid(planstr[arcs[kk]].filename,0,pos)
        endfor
       close, 2
       ;; Script
       printf, 12, 'xgmosArcProc("@'+filnm+'", "'+pwd+'/",'+$
               '"'+biasfil+'",'+$
               '"'+flatroot+mask_list[ii]+'")'

;       allproc = [allproc,arcs] 
       allproc = [allproc]
       filnm= 'allproc_'+mask_list[ii]+'.lst'
       openw, 2, filnm
       for kk=0L,n_elements(allproc)-1 do begin
           pos = strpos(planstr[allproc[kk]].filename, '.fits')
           printf, 2, strmid(planstr[allproc[kk]].filename,0,pos)
       endfor
       close, 2
       ;; Script
       printf, 12, 'xgmosProc("@'+filnm+'", "'+pwd+'/",'+$
               '"'+biasfil+'",'+$
               '"'+flatroot+mask_list[ii]+'")'
   endfor

   ;; Master
   printf, 10, 'cl < xgmos_flats.cl'
   printf, 10, 'cl < xgmos_proc.cl'

   close, /all

   print, 'xgmos_mkclscripts: Move the scripts (and lists) where you wish'
   print, 'xgmos_mkclscripts: Launch IRAF from the proper spot'
   print, 'xgmos_mkclscripts: Load gemini gmos x_gmos'
   print, 'xgmos_mkclscripts: Then try --   cl < xgmos_master.cl'
   print, 'xgmos_mkclscripts: If successful, proceed to IDL'
   ;; 

return
end
;------------------------------------------------------------------------------
