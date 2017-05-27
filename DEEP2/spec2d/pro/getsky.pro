;+
; NAME:
; getsky
;
; PURPOSE:
; finds all the sky slits on a mask and computes important information
; (e.g., mask position).  Also, can tweak and recompute bspline
; information if desired.
;
; CALLING SEQUENCE:
; getsky, maskn, datadir, color, sky_slitn 
;
; INPUTS:
; 
; maskn    - number of mask on which to find sky slitlets
; datadir  - data directory in which spSlit files for this
;            mask can be found. 
; color    - 'R' or 'B'
;
; OPTIONAL INPUTS:
;
; nbuf      - number of pixels by which to buffer slit ends before bspline	
;
; KEYWORDS:
;
; bspline  - If set, the routine performs a variable skytweak and runs
;            a new bspline for each sky slit (THIS SHOULD CHANGE TO A
;            BKPT TWEAK AFTER V0_10 IS RUN).  New bsplines are stored
;            in the local directory in files called sky_sset.slitn.frame.sav
; silent   - set to run silently.
;
; OUTPUTS:
;
; sky_slitn - an array containing the slit numbers of the sky-only slits.
;
; OPTIONAL OUTPUTS:
;
; sky_slitx - an array containing the x-coordinates of the sky slits
; sky_slity - three guesses...
; sky_minmaxwave - 2-d array containing min and max wavelength for
;                  each sky-only spectrum
; bkspace   - if bspline keyword is set, this is the spacing of the
;             bkpts passed to bspline_iterfit
; skyspfile - file names for the sky slit spSlit files.
; nexp      - number of exposures to expect in each spSlit file (default:3).
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
; Calls find_skyslits.pro
;
; REVISION HISTORY:
; Written by BFG in summer 02.
; Made presentable on 2002Oct18.
;
;----------------------------------------------------------------------

pro getsky,  maskn, datadir, color, sky_slitn, $
             sky_slitx=sky_slitx,  sky_slity=sky_slity, $
             sky_minmaxwave=sky_minmaxwave, nbuf=nbuf, $
             bspline=bspline,  bkspace=bkspace, skyspfile=skyspfile,$ 
             nexp=nexp, silent=silent
  

;---find sky only slit numbers.
  sky_slitn = find_skyslits(datadir, slitx=sky_slitx, slity=sky_slity) 

   usegz=1
   nungzip=n_elements(findfile(datadir+'spS*.fits'))
   ngzip=n_elements(findfile(datadir+'spS*.fits.gz'))

   if nungzip gt 5 then usegz=0 else if ngzip eq 0 then message,'No spSlit files found in directory '+datadir

  n_slits =  n_elements(sky_slitn)

  if not keyword_set(nexp) then nexp=3
;  if not keyword_set(dlam) then dlam=0.2
  
;---construct spSlit filenames and read in spSlit files.
  filestem = datadir+'spSlit.'$
              +string(maskn, format='(i4)')+'.'

  for i=0, n_slits-1 do begin
 
    if sky_slitn[i] lt 10 then begin
      spfile = filestem+'00'+string(sky_slitn[i],format='(i1)')+color+'.fits'
    endif else if sky_slitn[i] lt 100 then begin
      spfile = filestem+'0'+string(sky_slitn[i],format='(i2)')+color+'.fits'
    endif else begin
      spfile = filestem+string(sky_slitn[i],format='(i3)')+color+'.fits'
    endelse

    if usegz eq 1 then spfile=spfile+'.gz'

    if i eq 0 then skyspfile = replicate(spfile, n_slits)

      skyspfile[i] = spfile
      if not keyword_set(silent) then $
        print,  'Attepting to read file '+spfile+', HDU '+$
                  string(1,  format='(i1)')
	
      skyslit = mrdfits(spfile, 1,/SILENT,status=status)
	      if i eq 0 then sky_minmaxwave = fltarr(n_slits, 2)
	
      if status ge 0 then begin
;---do variable sky tweak

	      if keyword_set(bspline) then begin
	          if not keyword_set(silent) then $
        	    print,  'Attepting to read file '+spfile+', HDU '+$
                	  string(1,  format='(i1)')
	         skysset =  mrdfits(spfile, 2,/SILENT)

	          tweakfit = deimos_skytweak_var(skyslit, skysset, color)
	          skyslit.lambdax = tweakfit
	      endif

	      wave = lambda_eval(skyslit)

  hasinfo=total(tag_names(skyslit) eq 'INFOMASK') gt 0

  if hasinfo then whgood = where((skyslit.infomask AND 1b) eq 0b $
                AND skyslit.mask eq 0b) $
        else whgood = where(skyslit.mask eq 0b)

;---find min and max wavelengths for sky slit spectra

	      sky_minmaxwave[i, *] = minmax(wave[whgood]) 
                                        ;^throws out vignetted regions.

;---redo bsplines and save if called
;***THIS SHOULD CHANGE TO A BKPT TWEAK AFTER V0_10 HAS BEEN RUN   
	    if keyword_set(bspline) then begin
    
	      for j=0, nexp-1 do begin
	        if not keyword_set(silent) then $  
        	  print,  'Attepting to read file '+spfile+', HDU '+$
	                  string(2*(j+1)-1,  format='(i1)')
	        skyslit =  mrdfits(spfile, 2*(j+1)-1)
	        if not keyword_set(silent) then $
	          print,  'Attepting to read file '+spfile+', HDU '+$
                  string(2*(j+1),  format='(i1)')
	        skysset =  mrdfits(spfile, 2*(j+1))
	        tweakfit = deimos_skytweak_var(skyslit, skysset, color)
        	skyslit.lambdax = tweakfit
	        wave = lambda_eval(skyslit)
	        skyflux =  skyslit.flux
        	flag_cr, skyslit,  newskyivar
	        skymask = newskyivar eq 0
        	skyflux = djs_maskinterp(skyflux,  skymask, iaxis=1)

;-------get rid of any other bad pix        
	        newskyivar =  newskyivar*(skyslit.mask eq 0)
;-------buffer ends of slits
        	if keyword_set(nbuf) then begin
	          skyy =  n_elements(newskyivar[1, *])
        	  newskyivar[*, 0:nbuf-1] = 0
	          newskyivar[*, skyy-nbuf:skyy-1] = 0
	        endif

        	if not keyword_set(bkspace) then bkspace = 0.2
	        isort =  bsort(wave)
        	sortwave = wave[isort]
	        sortflux = skyflux[isort]
	        sortivar = newskyivar[isort]
        	sky_sset=  bspline_iterfit(sortwave,  $
				   sortflux, invvar=sortivar, $
                                   upper=30,  lower=30,  maxiter=3, $
                                   bkspace=bkspace, /silent)
	        save, sky_sset, file='sky_sset.'+string(i, format='(i1)') $
			+'.'+ string(j, format='(i1)')+color+'.sav'
	      endfor    
	    endif
	endif else begin

; if file does not exist, make sure we do not use it!
              sky_minmaxwave[i, *] = [1E5,-1E5]
 		sky_slitx[i]=1E6
		sky_slity[i]=1E6
	endelse
  endfor

end
  
