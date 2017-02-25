;+ 
; NAME:
;  wfc3_g280_reduce_qso
;
; PURPOSE:
;   This code performs a full reduction on a given QSO for WFC3/G280
;   grism data.  It is tuned for quasars, but might work well on other
;   objects (especially point sources). Reduces beam A and beam C.
;
; CALLING SEQUENCE:
;   
;   wfc3_g280_reduce_qso, img_fil, spec_fil, QADIR=qadir, SKYDIR=, FITSDIR=, $
;                          FIND=, FRACTION=, NPIX=, $
;                          GUESS=, NAME=, EXBOX=, SRCH=, $
;                          JXP_KLUDGE=, AXE1=, OLD=, EXT_BOX=, $
;                          RANGE=, NIMG=, BOX_SIZE=, $
;                          NBEAMS=, BADPIX=, NOFIT=, BOXCAR=, $
;                          SIGMA=, NOTWEAK=, WAVERG=, STRETCH=, $
;                          IMGSZ=, SAVESKY=, USEUPPERSKY=, $
;                          USELOWERSKY=, INTQA=, NOSHIFT=
;
; INPUTS:
;   img_fil -- Acquisition image file (usually has a flt extension)
;   spec_fil -- HST processed spectral image (usually has a crj
;               extension)
;
; RETURNS:
;
; OUTPUTS:
;   FITS file with the structure of the reduction, alternatively it
;   can return FITS file in 'standard format' with the NOTWEAK
;   keyword. Intermediate structure files are saved with a _g280I.fits
;   extension, whereas final files are saved with a _g280F.fits extension
;
; OPTIONAL KEYWORDS:
;   QADIR= -- Directory for QA output [default: 'QA/']
;   SKYDIR= -- Directory for sky output [default: 'Sky/']
;   FITSDIR= -- Directory for fits output [default: 'Fits/']
;   FIND= -- runs the find_sources program to find sources
;   FRACTION= -- the fraction that sets the threshold level to be
;                considered a detection
;   NPIX= -- the minimum number of pixels in an object to be
;            considered a real feature (not a hot pixel, etc.)
;   GUESS= -- Guess of the QSO position [default: [2047L,1085L,2]]
;   NAME= -- Name use for the qso for filename generation
;   EXBOX= -- Number of pixels to extend analysis for object
;             centroiding in wfc3_g280_center_direct.pro [default: 10] 
;   SRCH=  -- Number of pixels to extend search for QSO in direct
;             image in wfc3_g280_center_direct.pro [default: 10]
;   AXE1= -- Use the old AXE1 data (DEPRECATED)
;   JXP_KLUDGE= -- Offset in x for generating the wavelength solution
;                  [default: 3.0; based on eyeball comparison to SDSS
;                  spectra] (DEPRECATED)
;   OLD= -- Use the old (circa 2010) code (DEPRECATED)
;   EXT_BOX= -- Extraction box [default:10L]
;   RANGE= -- Range to calculate the tracewidth [default: [2700d,4500]]
;   NIMG=   -- Number of images used in crj stack (Default=2)
;   BOX_SIZE= -- Total width in pixels of the boxcar if doing boxcar ext. [default: 9L]
;   NBEAMS= -- Number of beams to reduce: 1=BEAM_A, 2=BEAM_A+BEAM_C
;   BADPIX= -- Threshold above which pixels are flagged bad in direct image
;   NOFIT= -- Do not centroid the direct image
;   BOXCAR= -- Do boxcar extraction in stead of optimal
;   SIGMA= -- Overwrite the sigma used for the optimal extraction
;   NOTWEAK= -- Do not perform a manual tweak on the image - output
;               directly 'final' files
;   WAVERG= -- Wavelength range to plot in QA file
;   STRETCH= -- Stretch to use in 2D image in QA file
;   IMGSZ= -- Image size of the direct image in QA file
;   SAVESKY= -- Save the sky image
;   USEUPPERSKY= -- Use the upper sky band for sky subtraction
;   USELOWERSKY= -- Use the lower sky band for sky subtraction
;   INTQA= -- Generate intermediate QA files (for debugging)
;   NOSHIFT= -- Do not apply a y-pixel shift in fitting the trace 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfc3_g280_reduce_qso, img_fil, spec_fil, QADIR=qadir, SKYDIR=skydir, FITSDIR=fitsdir, $
;                         FIND=find, FRACTION=fraction, NPIX=npix, $
;                         GUESS=guess, NAME=name, EXBOX=exbox, SRCH=srch, $
;                         JXP_KLUDGE=jxp_kludge, AXE1=axe1, OLD=old, EXT_BOX=ext_box, $
;                         RANGE=range, NIMG=nimg, BOX_SIZE=box_size, $
;                         NBEAMS=nbeams, BADPIX=badpix, NOFIT=nofit, BOXCAR=boxcar, $
;                         SIGMA=sigma, NOTWEAK=notweak, WAVERG=waverg, STRETCH=stretch, $
;                         IMGSZ=imgsz, SAVESKY=savesky, USEUPPERSKY=useuppersky, $
;                         USELOWERSKY=uselowersky, INTQA=intqa, NOSHIFT=noshift
;
; PROCEDURES CALLED:
;  wfc3_g280_find_sources
;  wfc3_g280_comb_beams
;  wfc3_g280_mantweak
;  wfc3_g280_mkfits
;  wfc3_g280_optimal
;  wfc3_g280_qa_final
;  wfc3_g280_qa_position
;  wfc3_g280_radec
;  wfc3_g280_center_direct
;  wfc3_g280_trace_wave
;  wfc3_g280_qa_trace
;  wfc3_g280_skysub
;  wfc3_g280_boxcar
;  wfc3_g280_flux
;
; REVISION HISTORY:
;   23-Dec-2010 Written by JXP/JMO
;   10-Jun-2016 Major rewrite if the code to allow for multiple beams
;               to be reduced, save the data in a structure, plus much
;               more by MN
;------------------------------------------------------------------------------
pro wfc3_g280_reduce_qso, img_fil, spec_fil, QADIR=qadir, SKYDIR=skydir, FITSDIR=fitsdir, $
                          FIND=find, FRACTION=fraction, NPIX=npix, $
                          GUESS=guess, NAME=name, EXBOX=exbox, SRCH=srch, $
                          JXP_KLUDGE=jxp_kludge, AXE1=axe1, OLD=old, EXT_BOX=ext_box, $
                          RANGE=range, NIMG=nimg, BOX_SIZE=box_size, $
                          NBEAMS=nbeams, BADPIX=badpix, NOFIT=nofit, BOXCAR=boxcar, $
                          SIGMA=sigma, NOTWEAK=notweak, WAVERG=waverg, STRETCH=stretch, $
                          IMGSZ=imgsz, SAVESKY=savesky, USEUPPERSKY=useuppersky, $
                          USELOWERSKY=uselowersky, INTQA=intqa, NOSHIFT=noshift
  
  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
          'wfc3_g280_reduce_qso, img_fil, spec_fil, [fin_strct], QADIR=, NAME=, ' + $
          'XGUESS=, YGUESS=, SEARCH= [v1.0]'
    return
  endif 

  if not keyword_set(SKYDIR) then skydir = 'Sky/'
  if not keyword_set(QADIR) then qadir = 'QA/'
  if not keyword_set(FITSDIR) then fitsdir = 'Fits/'
  if not keyword_set(SRCH) then srch = 10L
  if not keyword_set(EXBOX) then exbox = 10L
  if not keyword_set(BADPIX) then badpix = 1d10
  if not keyword_set(EXT_BOX) then ext_box = 10L
  if not keyword_set(BOX_SIZE) then box_size = 9L
  if not keyword_set(RANGE) then range = [2700d,4500] 
  if not keyword_set(NIMG) then NIMG=2 ;; Default for our program.
  if not keyword_set(NBEAMS) then $
     if keyword_set(OLD) then NBEAMS=1 else NBEAMS=2
  if keyword_set(OLD) and not keyword_set(guess) then guess=[2047L,1085L,2]
  
  ;; First define the structure
  wfc3_g280_strct=replicate({wfc3_g280_strct},1000)

  ;; Add the image/spectrum files
  wfc3_g280_strct(*).img_fil=img_fil
  wfc3_g280_strct(*).spec_fil=spec_fil

  ;; Find the sources
  if keyword_set(find) then begin
     wfc3_g280_find_sources, wfc3_g280_strct, fraction=fraction, npix=npix
  endif else begin
     if not keyword_set(guess) then begin
        print, 'wfc3_g280_reduce_qso: Need to supply some initial guesses.'
        return
     endif
     wfc3_g280_strct=wfc3_g280_strct(0:n_elements(guess)/3-1)
     wfc3_g280_strct(*).xguess=reform(guess(0,*))
     wfc3_g280_strct(*).yguess=reform(guess(1,*))
     wfc3_g280_strct(*).chip=reform(guess(2,*))
  endelse

  ;; Centroid the direct image
  if keyword_set(nofit) then begin
     wfc3_g280_strct(*).x0=wfc3_g280_strct(*).xguess
     wfc3_g280_strct(*).y0=wfc3_g280_strct(*).yguess
  endif else $
     wfc3_g280_center_direct, wfc3_g280_strct, EXTENDBOX=EXBOX, SRCH=srch, $
                              BADPIX=badpix
  
  ;; Assign a name and ra/dec of the object
  wfc3_g280_radec, wfc3_g280_strct
  ;; Overwrite the default names by the user supplied names
  if keyword_set(name) then $
     wfc3_g280_strct.name=name
  
  ;; QA for position in direct image
  if keyword_set(intqa) then $
     wfc3_g280_qa_position, wfc3_g280_strct, SRCH=srch, QADIR=qadir
  
  ;;;;;;
  ;; Loop through the sources 
  for ii=0L,n_elements(wfc3_g280_strct)-1 do begin
     
     ;; remove the objects whose trace falls outside the
     ;; spectral image completely
     print, wfc3_g280_strct[ii].name
     if wfc3_g280_strct(ii).y0 ge 1950. then begin
        print, 'wfc3_g280_reduce_qso: Trace falls outside the chip. Not'+$
               ' reducing this object'
        continue
     endif

     for beam_idx=0,nbeams-1 do begin ;;loop over the beams (currently only A and C)
        
        if beam_idx eq 0 then begin
           
           ;; Get information from the headers + read spec_image
           head=headfits(wfc3_g280_strct(ii).spec_fil,ext=0)
           wfc3_g280_strct.exptime=sxpar(head, 'EXPTIME')
           specim=xmrdfits(wfc3_g280_strct(ii).spec_fil, $
                           7-3*wfc3_g280_strct(ii).chip)
        endif
        
        ;; Trace     
        wfc3_g280_trace_wave, wfc3_g280_strct, ii, specim, AXE1=axe1, $
                              JXP_KLUDGE=jxp_kludge, $
                              CALIB=calib, OLD=old, EXT_BOX=ext_box, $
                              RANGE=range, BEAM=beam_idx, NOSHIFT=noshift

        ;; QA of Trace
        if keyword_set(intqa) then $
           wfc3_g280_qa_trace, wfc3_g280_strct, ii, specim, BEAM=beam_idx, $
                               QADIR=qadir

        ;; Sky Subtraction
        specim_sky = wfc3_g280_skysub(wfc3_g280_strct, ii, specim, $
                                      BEAM=beam_idx, SKYDIR=skydir,$
                                      savesky=savesky, useuppersky=useuppersky, $
                                      uselowersky=uselowersky)
        
        ;; Variance image
        var = (specim>1.) + NIMG*calib.readno^2

        ;; Extract spectrum
        if keyword_set(old) or keyword_set(boxcar) then $
           wfc3_g280_boxcar, wfc3_g280_strct, ii, specim_sky, var, $
                             BOX_SIZE=box_size, BEAM=beam_idx $
        else $
           wfc3_g280_optimal, wfc3_g280_strct, ii, specim_sky, var, $
                              BEAM=beam_idx, sigma=sigma
        
        ;; Flux
        wfc3_g280_flux, wfc3_g280_strct, ii, calib, OLD=old, $
                        BEAM=beam_idx        
      
     endfor

     ;; combine the beams
     wfc3_g280_comb_beams, wfc3_g280_strct, ii
     
     ;; Final QA
     wfc3_g280_qa_final, wfc3_g280_strct, ii, IMGSZ=imgsz, $
                        QADIR=qadir, WAVERG=waverg, STRETCH=stretch, $
                        SIGMA=sigma, BOXCAR=boxcar, BOX_SIZE=box_size

     ;; save the structure
     fitsname=fitsdir+wfc3_g280_strct(ii).name+'_g280I.fits'
     mwrfits, wfc3_g280_strct(ii), fitsname, /create

     if keyword_set(notweak) then $
        wfc3_g280_mantweak, wfc3_g280_strct(ii).name+'_g280I.fits', '0'
     
  endfor
  
end
