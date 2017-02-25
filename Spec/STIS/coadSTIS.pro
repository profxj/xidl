;+
; Name:
; coadSTIS
;   Version 1.1
;
; PURPOSE:
;    Automated procedure for coadding STIS 1D or echelle spectra
;
; CALLING SEQUENCE:
;
;   coadSTIS, infile, echelle=, nfits=, tol=, /print, view=
;
; INPUTS:
;   infile  - list of STIS spectra to coadd (e.g. *x1d.fits)
;
; RETURNS:
;
; OPTIONAL OUTPUTS:
;   [filename].fits - multi-extension FITS of 1D coadded spectrum
;                     (flux, error, data quality flags)
;   [filename]_F.fits/_E.fits - FITS of coadded flux/error
;   [filename]_O.fits - multi-extension FITS of coadded orders
;                       (flux, error, data qualtiy flags
;   Note: [filename] represents root of infile name
;
; OPTIONAL KEYWORDS:
;   echelle=  - coadd echelle spectral orders
;               = 2 - coadded orders and into 1D spectrum 
;   tol=      - S/N tolerance for weighting spectra (default, 0.5)
;   view=     - display coadded spectrum with error and compared to
;               reference spectrum (highest S/N) with x_splot
;               = 1 - for echelle spectra, display each coadded order
;                     as well as 1D spectrum
;               = 2 - for echelle spectra, display only 1D spectrum
;
; OPTIONAL OUTPUTS:
;   nfits=    - number of FITS to create
;               = 0 - coadded spectrum not saved
;               = 1 - Default; multi-extension FITS (flux, error, data
;                     quality flag); [filename].fits
;               = 2 - 2 single-extension FITS of flux and
;                     error (= -1 where data qualtiy flags non-zero);
;                     [filename]_F.fits/_E.fits 
;               = 3 - combination of 1 and 2
;               Note: [filename]_O.fits (flux, error and data quality
;               flags for coadded orders) created when echelle keyword
;               set (except for nfits=0)
;   /print    - print wavelength relation, flux scale and S/N for
;               spectra to [filename].wave/.flux/.snr, respectively
;   /generic - set to coadd already coadded spectra
;   /zqso - infile includes column of quasar redshifts to include in
;           generic coaddingXS
;
; COMMENTS:
;   - Output FITS contain modified header of reference spectrum
;     (highest S/N); includes history of coadding process
;
; EXAMPLES:
;   coadSTIS, FJ2155_0922_E140M, echelle=2, nfits=3, /print, view=2
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_specrebin
;  x_combspec
;  x_splot
;
; REVISION HISTORY:
;   23-August-2004:  Written by Kathy Cooksey
;   27-September-2004: Modified use of DQ by making them a mask
;                      initially so that some non-zero DQ would be 
;                      considered acceptable (16, 32, 1024). 
;                      More cleanly truncated data written to FITS.
;   07-September-2005: Under specified assumptions, now handles
;                      echelle spectra with differing orders more 
;                      robustly
;   12-September-2005: Changed default FITS names to reflect being
;                      unnormalized
;   17-October-2005: fixed bug with nfits=2, 1D spectra, flux scaling (biggie),
;                    and trimming coadded spectrum (even bigger)
;   18-October-2005: trim_array includes trimming -9.99; made default
;                    number -9.99 (not -999.99)
;   19-October-2005: calculate global_pix needed; set flux=0. in
;                    trouble spots
;   25-October-2005: attempted to add generic coaddition ability; 
;                    changed FILENAME in keyword header
;   22-August-2007: Increase trim_array nbins by 1, add compile measureSN
;   27-September-2007: Minor bug fixes
;   04-October-2007: read in multiple extensions (exposures); no FLOOR
;-
;------------------------------------------------------------------------------

function measureSN,spec,err=err,inflg=inflg,plot=plot,wave=wave,$
                   cont=cont,dev=dev,region=region,silent=silent
if not keyword_set(inflg) then inflg = 0
case inflg of
   4: begin
      flux = spec
      if not keyword_set(err) then stop,'measureSN: must set err'
      error = err
      if not keyword_set(wave) then stop,'measureSN: must set wave'
      svspec = {wave:wave,flux:spec,error:err} ; save to prevent mistakes
      npix = n_elements(wave)
   end
   else: flux = x_readspec(spec,inflg=inflg,head=hd,sig=error,wav=wave,$
                           fil_sig=err,npix=npix)
endcase 

;;For those freak cases where some element is NaN, eliminate it in
;;both flux and error
bd = where(finite(flux,/nan) or finite(error,/nan),nbd,complement=gd)
if nbd ne 0 then begin
    flux = flux[gd]
    error = error[gd]
    wave = wave[gd]
    if not keyword_set(silent) then $
      print,'measureSN: elements of flux or error NaN, trim ',spec
endif 

if not keyword_set(cont) then begin
;;Histogram data and determine peak
    div = 5.

;;To prevent trying to histogram large range of flux and/or error
;;limit the number of histogram bins by setting sztol
    sztol = 100*n_elements(flux)    

    if median(flux) le 0. then fluxbin = 5.e-15 $
    else fluxbin = median(flux)/div
    while (max(flux,iimx,min=mn,subscript_min=iimn)-mn)/fluxbin gt sztol $
      do begin
        if abs(flux[iimx]-fluxbin) gt abs(mn-fluxbin) then begin
            ;;eliminate max flux as larger discrepency
            gd = where(flux ne flux[iimx],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endif else begin
            ;;eliminate min flux as larger discrepency
            gd = where(flux ne flux[iimn],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endelse 
        if not keyword_set(silent) then $
          print,'measureSN: flux range too large; trim flux and error ',spec
    endwhile 

    if median(error) le 0 then errorbin = 5.e-16 $
    else errorbin = median(error)/div
    while (max(error,iimx,min=mn,subscript_min=iimn)-mn)/errorbin gt sztol $
      do begin
        if abs(error[iimx]-errorbin) gt abs(mn-errorbin) then begin
            ;;eliminate max error as larger discrepency
            gd = where(error ne error[iimx],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endif else begin
            ;;eliminate min error as larger discrepency
            gd = where(error ne error[iimn],ngd)
            flux = flux[gd]
            error = error[gd]
            wave = wave[gd]
        endelse 
        if not keyword_set(silent) then $
          print,'measureSN: error range too large; trim flux and error ',spec
    endwhile 

;;Limit to positive flux and corresponding error
;;Flux <= 0. is non-physical
    gd = where(flux gt 0.,ngd)  ;ge 0. causes funniness
    if ngd gt 0 then begin
       ;;center the histograms
       nbins = ceil((max(flux[gd],min=mn)-mn)/fluxbin)
       mn = median(flux[gd]) - 0.5*nbins*fluxbin
        fluxhist = histogram(flux[gd],binsize=fluxbin,locations=fluxloc,min=mn)
        if n_elements(fluxloc) le 3 then begin
            ;; Arrays must be big enough to fit
            fluxhist = [0,0,0,fluxhist,0,0,0]
            fluxloc = (lindgen(n_elements(fluxhist))-3)*fluxbin + fluxloc[0]
            if not keyword_set(silent) then $
              print,'measureSN: buffered flux histogram'
        endif 

       nbins = ceil((max(error[gd],min=mn)-mn)/errorbin)
       mn = median(error[gd]) - 0.5*nbins*errorbin
        errorhist = histogram(error[gd],binsize=errorbin,$
                              locations=errorloc,min=mn)
        if n_elements(errorloc) le 3 then begin
            ;; Arrays must be big enough to fit
            errorhist = [0,0,0,errorhist,0,0,0]
            errorloc = (lindgen(n_elements(errorhist))-3)*errorbin + errorloc[0]
            if not keyword_set(silent) then $
              print,'measureSN: buffered error histogram'
        endif 
        ;; Median values and deviation will be from fitted Gaussian
        ;; (b/c of measure_errors weighting, mdflux/error exactly the
        ;; same) 
        !quiet = 1 ;suppress messages
        gfit = gaussfit(fluxloc,fluxhist,coeff,nterms=3,$
                        measure_errors=sqrt(fluxhist))
        mdflux = coeff[1]
        devflux = coeff[2]
        gfit = gaussfit(errorloc,errorhist,coeff,nterms=3,$
                        measure_errors=sqrt(errorhist))
        mderror = coeff[1]
        deverror = coeff[2]
        snr = mdflux/mderror
        devsnr = snr*sqrt((devflux/mdflux)^2+(deverror/mderror)^2)
    endif else begin
        if not keyword_set(silent) then $
          print,'measureSN: flux <= 0. for ',spec
        mdflux = -999.99
        mderror = -999.99
        snr = -999.99
        devflux = -999.99
        deverror = -999.99
        devsnr = -999.99
    endelse
;stop
    if keyword_set(plot) then $
      x_splot,fluxloc,fluxhist,xtwo=errorloc,ytwo=errorhist,/block
    delvarx,fluxhist,errorhist


endif else begin
;;;;;;;;;;;;;;;
;; Continuum fit and eliminate absorption features
;;;;;;;;;;;;;;;
    safe = {wave:wave,flux:flux,error:error}

    if keyword_set(plot) then $
      autocontfit,wave,flux,error,buffer=2,binsize=7,region=region,$
      /refin,/silent,/view $
    else autocontfit,wave,flux,error,buffer=2,binsize=7,region=region,$
      /refin,/silent

    mdflux = median(flux[region],/even)
    mderror = median(error[region],/even)
    snr = mdflux/mderror

    devflux = stddev(flux[region])
    deverror = stddev(error[region])
    devsnr = snr*sqrt((devflux/mdflux)^2+(deverror/mderror)^2)

endelse                         ;/cont 

!QUIET = 0                      ;continue suppressing errors

;; Restore arrays
if inflg eq 4 then begin
   wave = svspec.wave
   spec = svspec.flux
   err = svspec.error
endif 

dev = [devflux,deverror,devsnr]
return,[mdflux,mderror,snr]

end


PRO trim_array,a,lead_zero,nbins
;Determines leading/trailing zeroes or -9.99 (defaults) in array a  
;Returns location of first non-zero element and resulting number of
;bins; array of indices used to actually truncate a can be created
;by: lead_zero + INDGEN(nbins)
size = N_ELEMENTS(a)
done = 0
nn = 0
lead_zero = -1
trail_zero = size
WHILE NOT done DO BEGIN
    IF (a[nn] NE 0D) and (a[nn] ne -9.99) $
      AND (lead_zero EQ -1) THEN lead_zero = nn
    IF (a[size-nn-1] NE 0D) and (a[size-nn-1] ne -9.99) $
      AND (trail_zero EQ size) THEN $
      trail_zero = size-nn-1
    IF (lead_zero NE -1) AND (trail_zero NE size) THEN done = 1 
    IF nn EQ size-1 THEN done = 1 ELSE nn = nn + 1
ENDWHILE
IF (lead_zero EQ -1) THEN lead_zero = 0
nbins = trail_zero - lead_zero + 1
END



FUNCTION BINARY, number
;http://www.mpimet.mpg.de/en/misc/software/idl/html/src/martin_schultz/binary.pro
;+
; NAME:
;   binary (function)
; PURPOSE:
;   This function returns the binary representation
;   of a number. Numbers are converted to LONG integers
;   if necessary.
; CATEGORY:
;   General Programming
; EXAMPLE:
;   Binary representation of 11B:
;     IDL> print, binary(11B)                                           
;     0 0 0 0 1 0 1 1
;   If data extraction is used instead of conversion ->
;   Binary representation of pi (little endian IEEE representation):
;     IDL> print, format='(z9.8,5x,4(1x,8a1))', long(!pi,0), binary(!pi)
;     40490fdb      01000000 01001001 00001111 11011011
; AUTHOR:
;      Kevin Ivory                           Tel: +49 5556 979 434
;      Max-Planck-Institut fuer Aeronomie    Fax: +49 5556 979 240
;      Max-Planck-Str. 2                     mailto:Kevin.Ivory@linmpi.mpg.de
;      D-37191 Katlenburg-Lindau, GERMANY
;
;-
  ;KC On_Error, 1
  ;KC s = SIZE(number)
  ;KC type = s[s[0] + 1]
  ;KC IF type EQ 0 THEN Message, 'Number parameter must be defined.'
bit = [0,1]                     ;KC ['0','1']
  ;KC IF type EQ 1 OR type EQ 2 THEN BEGIN
  ;KC  bitvalue = 2^INDGEN(8*type)
  ;KC ENDIF ELSE BEGIN
  ;KC  Print, 'Converting "number" to LONG...'
number = LONG(number)           ; data conversion
;   If you want the binary representation of the floating point value,
;   use extraction instead of conversion:
;   number = LONG(number, 0)    ; data extraction
bitvalue = 2L^LINDGEN(16)       ;KC (32)
  ;KC ENDELSE

  RETURN, REVERSE(bit((number AND bitvalue) EQ bitvalue))
END


FUNCTION makemask,dq
MAGIC_DQ_MASK = binary(64462)   ;DQ = {0,16,32,1024} acceptable
num = N_ELEMENTS(dq)
mask = REPLICATE(1B,num)
FOR nn=0L,num-1 DO BEGIN
    dq_bin = binary(dq[nn])
    test = WHERE((dq_bin*MAGIC_DQ_MASK) NE 0,n_nonzero)
    IF n_nonzero EQ 0 THEN mask[nn] = 0B
ENDFOR
RETURN,mask
END



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

PRO coadSTIS,infile,echelle=echelle,nfits=nfits,tol=tol,$
             print=print,view=view,generic=generic,zqso=zqso

;;root of files to be created
OUTFILE = STRMID(infile,0,STRPOS(infile,'.',/REVERSE_SEARCH)) 

caldat,systime(/julian),mo,da,yr
date=strtrim(da,2)+'.'+strtrim(mo,2)+'.'+strtrim(yr,2)

history = 'Modifications made with coadSTIS.pro '+date 
IF ARG_PRESENT(nfits) THEN nfits = 1 ;default

;;;;;;;;;;;;;;;;
; Skip to /generic option
;;;;;;;;;;;;;;;;
if keyword_set(generic) then goto,gnrc

;;;;;;;;;;;;;;;
; Read list of spectra to coadd
;;;;;;;;;;;;;;;
N_SPEC = FILE_LINES(infile)
readcol,infile,files,format='A',/silent
PRINT,'Coadding from spectra listed in ',infile

;;minimum S/N used to scale data
IF KEYWORD_SET(tol) THEN TOLERANCE = tol ELSE TOLERANCE = 0.5 
;BUFFER = 50                     ;extra pixels for rebinning
;;minimum number of pixels and orders accepted to define quantity
;;(e.g. S/N) 
N_PIX_MIN = 100                
N_ORDER_MIN = 20                ;might need to be adjusted for grating
;;Notes to be written to fits file

history = [history,'Program paramters:']
history = [history,'S/N tolerance = '+STRING(TOLERANCE)]
history = [history,'Minimum number of pixels used to define S/N and '+$
           'flux scaling: '+STRING(N_PIX_MIN)]
IF KEYWORD_SET(echelle) THEN $
  history = [history,'Minimum number of orders used to define S/N and '+$
             'flux scaling: '+STRING(N_ORDER_MIN)]
history = [history,'Flux and error measured in units of '+$
           'erg/s/cm**2/Angstrom']

IF N_SPEC LE 1 THEN BEGIN
    ;;In case supplSTIS.pro makes *.list with only one spectrum listed
    sxhread,files[0],hd
    num = fxpar(hd,'NEXTEND')
    IF num EQ 1 THEN BEGIN
        PRINT,'EXITING: None or one spectra listed in ',infile
        GOTO,QUIT
    ENDIF
ENDIF


;;;;;;;;;;;;;;;
; Instantiate quasar structure from reduced STIS data
;;;;;;;;;;;;;;;
quasar = PTRARR(100*N_SPEC,/NOZERO)
nwfiles = strarr(100*N_SPEC)
history = [history,'Coadd spectrum: '+files]
;;Data element contains: wavelength, flux, error and data quality
;;flags (good EQ 0)
TOT_EXP = 0
TOT_SPEC = 0
FOR ii=0,N_SPEC-1 DO BEGIN
    sxhread,files[ii],hd
    num = fxpar(hd,'NEXTEND')
    TOT_EXP = TOT_EXP + fxpar(hd,'EXPTIME')
    FOR jj=1,num DO BEGIN       ; accessing extensions > 0
        dat = xmrdfits(files[ii],jj,/silent)
        quasar[TOT_SPEC+jj-1] = PTR_NEW(dat,/no_copy)
    ENDFOR                      ;loop sub-array
    nwfiles[TOT_SPEC:TOT_SPEC+num-1] = files[ii]
    TOT_SPEC = TOT_SPEC+num     ;increment
ENDFOR
;; Trim
N_SPEC = TOT_SPEC
quasar = quasar[0:tot_spec-1]
files = nwfiles[0:tot_spec-1]


IF KEYWORD_SET(echelle) THEN BEGIN
;;;;;;;;;;;;;;;
; COADD ECHELLE SPECTRA
;;;;;;;;;;;;;;;
    N_SPORDER = N_ELEMENTS((*quasar[0]).sporder)
    ORDER_O = (*quasar[0])[0].sporder
    problem_spec = -1
    FOR ii=1,N_SPEC-1 DO BEGIN
        tmp = N_ELEMENTS((*quasar[ii]).sporder)
        IF tmp LT N_SPORDER THEN BEGIN
            PRINT,'WARNING: Number of orders in ',files[ii],$
              ' not same in all listed spectra'
            stop
                                ;GOTO,QUIT
            problem_spec = ii
            N_SPORDER = N_ELEMENTS((*quasar[ii]).sporder)
            ORDER_O = (*quasar[ii])[0].sporder
        ENDIF
    ENDFOR 
    ;;The following is a template for fixing the problem where all
    ;;spectra do not have the same number of orders
    ;;Assume the beginning orders (longer wavelengths) are trouble ones
    IF problem_spec NE -1 THEN BEGIN
        FOR ii=0,N_SPEC-1 DO BEGIN
            ;;Must verify that no more spectra are problematic
            IF n_elements((*quasar[ii]).sporder) ne n_sporder THEN BEGIN
                ;;Shift higher orders (shorter wavelengths) to
                ;;beginning to match problematic spectra; must access
                ;;appropriately in future (since all not same size)
                diff = n_sporder - n_elements((*quasar[ii]).sporder) 
                (*quasar[ii]).wavelength = $
                  SHIFT((*quasar[ii]).wavelength,0,diff)
                (*quasar[ii]).flux = SHIFT((*quasar[ii]).flux,0,diff)
                (*quasar[ii]).error = SHIFT((*quasar[ii]).error,0,diff)
                (*quasar[ii]).dq = SHIFT((*quasar[ii]).dq,0,diff)
            ENDIF
        ENDFOR
        history = [history,files[problem_spec]+' has less orders; '+$
                   'other spectra shifted and trimmed to match']
        n_pix = (*quasar[problem_spec])[0].nelem
    ENDIF else N_PIX = (*quasar[0])[0].nelem
    history = [history,STRTRIM(STRING(N_SPORDER),2)+$
               ' orders beginning with '+STRTRIM(STRING(ORDER_O),2)]
    

;;;;;;;;;;;;;;;
; Manipulate DQ so that it's 0's (good) and 1's (bad)
;;;;;;;;;;;;;;;
    FOR ii=0,N_SPEC-1 DO BEGIN
        FOR jj=0,N_SPORDER-1 DO BEGIN
            mask = makemask((*quasar[ii])[jj].dq)
            (*quasar[ii])[jj].dq = mask
        ENDFOR
    ENDFOR
    history = [history,'Acceptable STIS DQ include any combination of '+$
               '0, 16, 32, and 1024']

;;;;;;;;;;;;;;;
; Determine reference spectrum (highest S/N)
;;;;;;;;;;;;;;;
    ;;Median S/N for each order of each spectra and median S/N for
    ;;entire spectrum (stored in [*,N_SPORDER])
    sn_ratio = REPLICATE(-9.99,N_SPEC,N_SPORDER+1)
    ;;Median lambda relation for each order of each spectra and median
    ;;dlambda/lambda  (i.e. (lambda_i+1 - lambda_i)/lambda_i+1) for
    ;;                                    entire spectrum 
    lambda_rel = REPLICATE(-9.99,N_SPEC,N_SPORDER+1)
    FOR jj=0,N_SPORDER-1 DO BEGIN
        FOR ii=0,N_SPEC-1 DO BEGIN
            good = WHERE(((*quasar[ii])[jj].dq EQ 0) AND $
                         (*quasar[ii])[jj].error NE 0D,n_good)
            IF n_good GE N_PIX_MIN THEN sn_ratio[ii,jj] = $
              MEDIAN((*quasar[ii])[jj].flux[good]/$
                     (*quasar[ii])[jj].error[good],/EVEN)
            lambda_rel[ii,jj] = $
              MEDIAN(((SHIFT((*quasar[ii])[jj].wavelength,-1))[0:N_PIX-2]-$
                      (*quasar[ii])[jj].wavelength[0:N_PIX-2])/$
                     (SHIFT((*quasar[ii])[jj].wavelength,-1))[0:N_PIX-2],/EVEN)
                                ;or (*quasar[ii])[jj].wavelength[0:N_PIX-2],/EVEN)
        ENDFOR
    ENDFOR 
    ;;Compare spectra to determine absolute wavelength range and
    ;;spectrum with highest S/N
    REF_SPEC_INDX = 0
    LAMBDA_MIN = MIN((*quasar[REF_SPEC_INDX]).wavelength,MAX=LAMBDA_MAX)
    FOR ii=0,N_SPEC-1 DO BEGIN
        snr_good = WHERE(sn_ratio[ii,0:N_SPORDER-1] NE -9.99,n_snr_good)
        IF n_snr_good GE N_ORDER_MIN THEN BEGIN
            sn_ratio[ii,N_SPORDER] = MEDIAN(sn_ratio[ii,snr_good],/EVEN)
            lambda_rel[ii,N_SPORDER] = MEDIAN(lambda_rel[ii,snr_good],/EVEN)
        ENDIF 
        IF sn_ratio[ii,N_SPORDER] GT sn_ratio[REF_SPEC_INDX,N_SPORDER] THEN $
          REF_SPEC_INDX = ii
        test_min = MIN((*quasar[ii]).wavelength,MAX=test_max)
        LAMBDA_MIN = test_min < LAMBDA_MIN
        LAMBDA_MAX = test_max > LAMBDA_MAX
    ENDFOR 

    ;;global wavelength array needs to only encompass shortest spectra
    if problem_spec ne -1 then $
      lambda_min = min((*quasar[problem_spec]).wavelength,max=lambda_max)
    history = [history,'Reference spectrum: '+files[REF_SPEC_INDX]+$
               ' (highest S/N)']
    history = [history,'Header from reference spectrum ']


;;;;;;;;;;;;;;;
; Generate new wavelength array for re-binning
;;;;;;;;;;;;;;;
;    GLOBAL_PIX = N_SPORDER*N_PIX + BUFFER
    dlambda_lambda = MEDIAN(lambda_rel[*,N_SPORDER],/EVEN)
    CDELT1 = ALOG10(1+dlambda_lambda) 
    GLOBAL_PIX = ceil(alog10(lambda_max/lambda_min)/cdelt1)+1
    global_wave = 10^(ALOG10(LAMBDA_MIN)+DINDGEN(GLOBAL_PIX)*CDELT1)


;;;;;;;;;;;;;;;
; Rebin spectra 
;;;;;;;;;;;;;;;
    rebin_spec = {wavelength:PTRARR(N_SPEC,N_SPORDER,/NOZERO),$
                  flux:PTRARR(N_SPEC,N_SPORDER,/NOZERO),$
                  error:PTRARR(N_SPEC,N_SPORDER,/NOZERO),$
                  dq:PTRARR(N_SPEC,N_SPORDER,/NOZERO)}
    local_pix = 0               ;need max number of bins in any order
    FOR jj=0,N_SPORDER-1 DO BEGIN
        ;;corresponding order of spectra to be on same wavelength scale 
        local_min = MIN((*quasar[REF_SPEC_INDX])[jj].wavelength,$
                        MAX=local_max)
        wave_good = WHERE((global_wave GT FLOOR(local_min)) AND $
                          (global_wave LT CEIL(local_max)),n_wave_good)
        local_pix = n_wave_good > local_pix
        FOR ii=0,N_SPEC-1 DO BEGIN
            var = ((*quasar[ii])[jj].error)^2
            var_bad = WHERE((*quasar[ii])[jj].dq NE 0,n_var_bad)
            IF n_var_bad NE 0 THEN var[var_bad] = 0. 
            x_specrebin,(*quasar[ii])[jj].wavelength,$
              (*quasar[ii])[jj].flux,$
              global_wave[wave_good],new_flx,$
              VAR=var,NWVAR=new_var,/SILENT
            new_dq = REPLICATE(0B,N_ELEMENTS(new_var))
            new_dq_bad = WHERE(new_var LE 0,n_new_dq_bad)
            IF n_new_dq_bad NE 0 THEN BEGIN
                ;;Neighbors of bad pixels bad
                new_dq[(new_dq_bad-1)[1:n_new_dq_bad-1]] = 25B
                new_dq[(new_dq_bad+1)[0:n_new_dq_bad-2]] = 25B
                ;;DQ if new_var EQ 0 or -1
                new_dq[new_dq_bad] = 17B
                ;;Change new_var = -1 to 0
                new_var[new_dq_bad] = 0.
            ENDIF 
            rebin_spec.wavelength[ii,jj] = PTR_NEW(global_wave[wave_good])
            rebin_spec.flux[ii,jj] = PTR_NEW(new_flx,/NO_COPY)
            rebin_spec.error[ii,jj] = PTR_NEW(SQRT(new_var),/NO_COPY)
            rebin_spec.dq[ii,jj] = PTR_NEW(new_dq,/NO_COPY)
        ENDFOR
        ;;For comparisons/debugging, will make giant array of original
        ;;reference spectrum 
        IF jj EQ 0 THEN BEGIN
            ref_spec_wave = (*quasar[REF_SPEC_INDX])[jj].wavelength
            ref_spec_flux = (*quasar[REF_SPEC_INDX])[jj].flux
            ref_spec_err = (*quasar[REF_SPEC_INDX])[jj].error
            ref_spec_dq = (*quasar[REF_SPEC_INDX])[jj].dq
            rebin_spec_wave = *rebin_spec.wavelength[REF_SPEC_INDX,jj]
            rebin_spec_flux = *rebin_spec.flux[REF_SPEC_INDX,jj]
            rebin_spec_err = *rebin_spec.error[REF_SPEC_INDX,jj]
            rebin_spec_dq = *rebin_spec.dq[REF_SPEC_INDX,jj]
        ENDIF ELSE BEGIN
            ref_spec_wave = [ref_spec_wave,$
                             (*quasar[REF_SPEC_INDX])[jj].wavelength]
            ref_spec_flux = [ref_spec_flux,$
                             (*quasar[REF_SPEC_INDX])[jj].flux]
            ref_spec_err = [ref_spec_err,$
                            (*quasar[REF_SPEC_INDX])[jj].error]
            ref_spec_dq = [ref_spec_dq,(*quasar[REF_SPEC_INDX])[jj].dq]
            rebin_spec_wave = [rebin_spec_wave,$
                               *rebin_spec.wavelength[REF_SPEC_INDX,jj]]
            rebin_spec_flux = [rebin_spec_flux,$
                               *rebin_spec.flux[REF_SPEC_INDX,jj]]
            rebin_spec_err = [rebin_spec_err,$
                              *rebin_spec.error[REF_SPEC_INDX,jj]]
            rebin_spec_dq = [rebin_spec_dq,$
                             *rebin_spec.dq[REF_SPEC_INDX,jj]]
        END
    ENDFOR
    order = SORT(ref_spec_wave)
    ref_spec_wave = ref_spec_wave[order]
    ref_spec_flux = ref_spec_flux[order]
    ref_spec_err = ref_spec_err[order]
    ref_spec_dq = ref_spec_dq[order]
    ref_spec_good = WHERE(ref_spec_err NE 0D,n_ref_spec_good)
    ref_spec_snr = REPLICATE(-1,N_ELEMENTS(ref_spec_wave))
    IF n_ref_spec_good NE 0 THEN ref_spec_snr[ref_spec_good] = $
      ref_spec_flux[ref_spec_good]/ref_spec_err[ref_spec_good]
    order = SORT(rebin_spec_wave)
    rebin_spec_wave = rebin_spec_wave[order]
    rebin_spec_flux = rebin_spec_flux[order]
    rebin_spec_err = rebin_spec_err[order]
    rebin_spec_dq = rebin_spec_dq[order]
    rebin_spec_good = WHERE(rebin_spec_err NE 0D,n_rebin_spec_good)
    rebin_spec_snr = REPLICATE(-1,N_ELEMENTS(rebin_spec_wave))
    IF n_rebin_spec_good NE 0 THEN rebin_spec_snr[rebin_spec_good] = $
      rebin_spec_flux[rebin_spec_good]/rebin_spec_err[rebin_spec_good]

    history = [history,'All orders rebinned to same '+$
               'wavelength scale with x_specrebin.pro']
    history = [history,'DQ = 17 if rebinned spectrum had variance <= 0']
    history = [history,'DQ = 25 for edges of variance <= 0']


;;;;;;;;;;;;;;;
; Determine flux scaling
;;;;;;;;;;;;;;;
    flux_scale = REPLICATE(-9.99,N_SPEC,N_SPORDER+1)
    FOR jj=0,N_SPORDER-1 DO BEGIN
        FOR ii=0,N_SPEC-1 DO BEGIN
            ;;scale on: (a) good data (ii-th and reference spetrum
            ;;with DQ = 0 and  ii-th spectrum new_var > 0); (b)
            ;;prevent dividing by zero; and (c) only scale on areas
            ;;with S/N GE TOLERANCE 
            good = WHERE((*rebin_spec.dq[ii,jj] EQ 0) AND $
                         ((*rebin_spec.flux[REF_SPEC_INDX,jj]) NE 0D) $
                         AND ((*rebin_spec.error[ii,jj]) NE 0D) AND $
                         ((*rebin_spec.flux[ii,jj])/$ 
                          (*rebin_spec.error[ii,jj]) GE TOLERANCE),n_good)
            IF n_good GE N_PIX_MIN THEN $
              flux_scale[ii,jj] = $
              MEDIAN((*rebin_spec.flux[ref_spec_indx,jj])[good]/$
                     (*rebin_spec.flux[ii,jj])[good],/EVEN) 
        ENDFOR                 
    ENDFOR                    
    FOR ii=0,N_SPEC-1 DO BEGIN
        ;;scale for spectra stored in last bin
        good = WHERE(flux_scale[ii,0:N_SPORDER-1] NE -9.99,n_good)
        IF n_good GE N_ORDER_MIN THEN $
          flux_scale[ii,N_SPORDER] = MEDIAN(flux_scale[ii,good],/EVEN)
    ENDFOR 


;;;;;;;;;;;;;;;
; Print wavelength relation, flux ratio/scale, and S/N
;;;;;;;;;;;;;;;
    IF KEYWORD_SET(print) THEN BEGIN
        OPENW,15,OUTFILE+'.wave'
        OPENW,16,OUTFILE+'.flux'
        OPENW,17,OUTFILE+'.snr'
        FOR jj=0,N_SPORDER DO BEGIN
            FOR ii=0,N_SPEC-1 DO BEGIN
                CASE ii OF
                    N_SPEC-1: BEGIN
                        PRINTF,15,FORMAT='(D11.6,TR3)',lambda_rel[ii,jj]
                        PRINTF,16,FORMAT='(D11.6,TR3)',flux_scale[ii,jj]
                        PRINTF,17,FORMAT='(D11.6,TR3)',sn_ratio[ii,jj]
                    END
                    0: BEGIN
                        PRINTF,15,FORMAT='(I2,TR3,D11.6,TR3,$)',jj+1,$
                          lambda_rel[ii,jj]
                        PRINTF,16,FORMAT='(I2,TR3,D11.6,TR3,$)',jj+1,$
                          flux_scale[ii,jj]
                        PRINTF,17,FORMAT='(I2,TR3,D11.6,TR3,$)',jj+1,$
                          sn_ratio[ii,jj]
                    END
                    ELSE: BEGIN
                        PRINTF,15,FORMAT='(D11.6,TR3,$)',lambda_rel[ii,jj]
                        PRINTF,16,FORMAT='(D11.6,TR3,$)',flux_scale[ii,jj]
                        PRINTF,17,FORMAT='(D11.6,TR3,$)',sn_ratio[ii,jj]
                    END
                ENDCASE
            ENDFOR
        ENDFOR
        CLOSE,15
        history = [history,'NOTE: '+OUTFILE+'.wave created']
        CLOSE,16
        history = [history,'NOTE: '+OUTFILE+'.flux created']
        CLOSE,17
        history = [history,'NOTE: '+OUTFILE+'.snr created']
    ENDIF 
    

;;;;;;;;;;;;;;;
; Coadd orders
;;;;;;;;;;;;;;;
    ;;this structure is oversized and will be reassigned later
    tmp_coad_order = {wavelength:REPLICATE(0.d,N_SPORDER,local_pix),$
                      flux:REPLICATE(0.d,N_SPORDER,local_pix),$
                      error:REPLICATE(-1.,N_SPORDER,local_pix),$
                      dq:REPLICATE(21B,N_SPORDER,local_pix)}
    ;;Need to scale all orders by same flux
    scale_echelle = flux_scale[*,N_SPORDER]
    history = [history,'Flux scale for each spectra:']
    history = [history,files+STRING(scale_echelle)]
    history = [history,'Each order weighted by its median S/N']

    new_local_pix = 0  ;need max number of bins with flux in any order
    new_lambda_min = LAMBDA_MAX ;looking for minimum for CRVAL1
    FOR jj=0,N_SPORDER-1 DO BEGIN
        n_pix_order = N_ELEMENTS(*rebin_spec.wavelength[REF_SPEC_INDX,jj])
        flux_order = DBLARR(n_pix_order,N_SPEC,/NOZERO)
        var_order = DBLARR(n_pix_order,N_SPEC,/NOZERO)
        ;;Need to weight each orders with its S/N
        snr_order = sn_ratio[*,jj]
        FOR ii=0,N_SPEC-1 DO BEGIN
            flux_order[*,ii] = *rebin_spec.flux[ii,jj]
            var_order[*,ii] = (*rebin_spec.error[ii,jj])^2
            err_bad = WHERE(*rebin_spec.dq[ii,jj] NE 0,n_err_bad)
            IF n_err_bad NE 0 THEN var_order[err_bad,ii] = 0.
        ENDFOR
        x_combspec,flux_order,var_order,flux_comb,var_comb,$
          SCALE=scale_echelle,SNR=snr_order
        dq_comb = REPLICATE(0B,N_ELEMENTS(var_comb))
        var_bad = WHERE(var_comb LE 0,n_var_bad)
        IF n_var_bad NE 0 THEN BEGIN
            ;;Change var_comb = -1 to 0
            var_comb[var_bad] = 0.
            dq_comb[var_bad] = 19B
        ENDIF

        ;;Truncate leading/trailing zeroes (because new_spec wavelength
        ;;> original wavelength) Remaining zeroes changed to -9.99 later
        ;;Can't do this 'spec_good' stuff because will assign features to
        ;;wrong wavelengths when 'spec_good' conditions not met in spectrum
;        order_good = WHERE((flux_comb NE -9.99) AND (var_comb NE 0.) AND $
;                           (dq_comb EQ 0),n_order_good)
;        IF n_order_good NE 0 THEN BEGIN
;            trim_array,flux_comb[order_good],lead_zero,n_flux_good
;            flux_good = lead_zero + INDGEN(n_flux_good)
;            flux_good = order_good[flux_good]
;        ENDIF ELSE BEGIN
        trim_array,flux_comb,lead_zero,n_flux_good
        flux_good = lead_zero + INDGEN(n_flux_good)
;        ENDELSE
        ;;orders will be different sizes due to truncation 
        new_local_pix = n_flux_good > new_local_pix
        flux_comb = flux_comb[flux_good]
        wave_comb = (*rebin_spec.wavelength[REF_SPEC_INDX,jj])[flux_good]
        error_comb = SQRT(var_comb[flux_good])
        dq_comb = dq_comb[flux_good]

        new_lambda_min = MIN(wave_comb) < new_lambda_min
        tmp_pix = N_ELEMENTS(wave_comb)
        tmp_coad_order.wavelength[jj,0:tmp_pix-1] = wave_comb
        tmp_coad_order.flux[jj,0:tmp_pix-1] = flux_comb
        tmp_coad_order.error[jj,0:tmp_pix-1] = error_comb
        tmp_coad_order.dq[jj,0:tmp_pix-1] = dq_comb
        ;;For comparisons/debugging, will make giant array of coadded
        ;;echelle spectrum 
        IF jj EQ 0 THEN BEGIN
            coad_order_wave = wave_comb
            coad_order_flux = flux_comb
            coad_order_err = error_comb
            coad_order_dq = dq_comb
        ENDIF ELSE BEGIN
            coad_order_wave = [coad_order_wave,wave_comb]
            coad_order_flux = [coad_order_flux,flux_comb]
            coad_order_err = [coad_order_err,error_comb]
            coad_order_dq = [coad_order_dq,dq_comb]
        ENDELSE

        IF KEYWORD_SET(view) THEN BEGIN
            IF view EQ 1 THEN BEGIN
                PRINT,'Coadd order: ',ORDER_O+jj
                PRINT,'Plotting combined spectrum (black) and error (red)'
                order_good = WHERE(tmp_coad_order.dq EQ 0,n_order_good)
                IF n_order_good NE 0 THEN BEGIN
                    x_splot,tmp_coad_order.wavelength[jj,order_good],$
                      tmp_coad_order.flux[jj,order_good],$
                      ytwo=tmp_coad_order.error[jj,order_good],/block
                    PRINT,'Plotting combined spectrum (black) and original ',$
                      'reference spectrum (red)'
                    ref_good = $
                      WHERE((*quasar[REF_SPEC_INDX])[jj].dq EQ 0,n_ref_good)
                    IF n_ref_good NE 0 THEN $
                      x_splot,tmp_coad_order.wavelength[jj,order_good],$
                      tmp_coad_order.flux[jj,order_good],$
                      xtwo=(*quasar[REF_SPEC_INDX])[jj].wavelength[ref_good],$
                      ytwo=(*quasar[REF_SPEC_INDX])[jj].flux[ref_good],$
                      /block $
                    ELSE PRINT,'No plot because no good data in reference'
                ENDIF ELSE PRINT,'No plot because no good data in coadded'
                IF jj EQ 0 THEN history = $
                  [history,'Coadded orders verified visually']
            ENDIF 
        ENDIF 
    ENDFOR 
    history = [history,'DQ = 19 if coadded variance <= 0']
    ;;Truncate tmp_coad_order because eliminated leading and trailing
    ;;zeroes in flux 
    coad_order = {wavelength:DBLARR(N_SPORDER,new_local_pix),$
                  flux:DBLARR(N_SPORDER,new_local_pix),$
                  error:DBLARR(N_SPORDER,new_local_pix),$
                  dq:DBLARR(N_SPORDER,new_local_pix)}
    STRUCT_ASSIGN,tmp_coad_order,coad_order
    history = [history,'Truncated FLUX to eliminate leading/trailing zeroes']
    history = [history,'Regions with no data set to 0.d with DQ = 21'+$
               ' in coadded orders']
    order = SORT(coad_order_wave)
    coad_order_wave = coad_order_wave[order]
    coad_order_flux = coad_order_flux[order]
    coad_order_err = coad_order_err[order]
    coad_order_dq = coad_order_dq[order]
    coad_order_good = WHERE(coad_order_err NE 0D,n_coad_order_good)
    coad_order_snr = REPLICATE(-1,N_ELEMENTS(coad_order_wave))
    IF n_coad_order_good NE 0 THEN coad_order_snr[coad_order_good] = $
      coad_order_flux[coad_order_good]/coad_order_err[coad_order_good]
    history = [history,'Orders coadded']

;;;;;;;;;;;;;;;
; Write orders to FITS (template: esi_echcoaddfin.pro)
;;;;;;;;;;;;;;;
    IF nfits NE 0 THEN BEGIN
        outfile_echelle = OUTFILE+'_O.fits'
        tmp = mrdfits(files[REF_SPEC_INDX],0,head,/silent)
        sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
          files[REF_SPEC_INDX]
        sxaddpar,head,'COMMENT',' '
        sxaddpar,head,'COMMENT',' *** Added Keywords *** '
        sxaddpar,head,'COMMENT',' '
        sxaddpar, head, 'CRVAL1', ALOG10(new_lambda_min)
        sxaddpar, head, 'CDELT1', CDELT1
        sxaddpar, head, 'CRPIX1', 1   
        sxaddpar, head, 'CTYPE1', 'LINEAR' 
        sxaddpar, head, 'DC-FLAG', 1  
        sxaddpar, head, 'BITPIX', -32 
        sxaddpar, head, 'TOT_EXP', TOT_EXP,$
          'total exposure time (seconds)--calculated'
                                ;sxaddpar, head, 'NAXIS', 2  
                                ;sxaddpar, head, 'NAXIS1', N_SPORDER
                                ;sxaddpar, head, 'NAXIS2', new_local_pix
        sxaddhist,history,head
        mwrfits, coad_order.flux, outfile_echelle, head, /create, /silent
        mwrfits, coad_order.error, outfile_echelle, /silent
        mwrfits, coad_order.wavelength, outfile_echelle, /silent
        mwrfits, coad_order.dq, outfile_echelle, /silent
        PRINT,'NOTE: ',outfile_echelle,' created'
    ENDIF 


;;;;;;;;;;;;;;;
; Coadd spectrum orders into 1D spectrum
;;;;;;;;;;;;;;;
    IF NOT STRCMP(echelle,1) THEN BEGIN
        ;;Template set by Gabe's mike_1dspec.pro; basically taking a
        ;;weighted mean to get coadded flux 
        n_pix_echelle = N_ELEMENTS(coad_order.wavelength)
        tot_flux = REPLICATE(0.d,n_pix_echelle)
        tot_wave = 10^(ALOG10(new_lambda_min)+DINDGEN(n_pix_echelle)*CDELT1)
        tot_err = REPLICATE(0.d,n_pix_echelle)
        weight = DBLARR(n_pix_echelle)
        tot_dq = REPLICATE(23B,n_pix_echelle)
        FOR jj=0,N_SPORDER-1 DO BEGIN
            indx = WHERE(ABS(tot_wave - coad_order.wavelength[jj,0]) LT $
                         1.5E-4,n_indx)
            IF n_indx EQ 0 THEN BEGIN
                PRINT,'EXITING: Echelle wavelength does not coincide with ',$
                  'global wavelength, order index:',jj
                GOTO,QUIT  
            ENDIF
            indx = indx[0] 
            var_good = WHERE(coad_order.dq[jj,*] EQ 0,n_var_good)
            IF n_var_good NE 0 THEN BEGIN
                tmp_var = (coad_order.error[jj,var_good])^2
                tmp_flux = tot_flux[indx+var_good]+$
                  coad_order.flux[jj,var_good]/tmp_var
                FOR nn=0L,n_var_good-1 DO BEGIN
                    tot_flux[indx+var_good[nn]] = tmp_flux[nn]
                    weight[indx+var_good[nn]] = weight[indx+var_good[nn]] + $
                      1./tmp_var[nn]
                    tot_dq[indx+var_good[nn]] = coad_order.dq[jj,var_good[nn]]
                ENDFOR 
            ENDIF 
        ENDFOR 
        weight_good = WHERE(weight NE 0D,n_weight_good)
        IF n_weight_good EQ 0 THEN BEGIN
            PRINT,'EXITING: Failed to coadd orders into 1D spectrum'
            GOTO,QUIT  
        ENDIF            
        tot_flux[weight_good] = tot_flux[weight_good]/weight[weight_good]
        tot_err[weight_good] = SQRT(1./weight[weight_good])
        tot_dq[weight_good] = 0B
        history = [history,'Flux = 0.d, DQ = 23 in bad regions'+$
                   ' of coadded 1D spectrum']
        ;;Truncate leading/trailing zeroes (because coadded spectrum has less
        ;;pixels than n_pix_echelle
        ;;Can't do this 'spec_good' stuff because will assign features to
        ;;wrong wavelengths when 'spec_good' conditions not met in spectrum
        trim_array,tot_flux,lead_zero,n_flux_good
        flux_good = lead_zero + LINDGEN(n_flux_good)
        coad_echelle = {wavelength:tot_wave[flux_good],$
                        flux:tot_flux[flux_good],error:tot_err[flux_good],$
                        dq:tot_dq[flux_good]}
        IF KEYWORD_SET(view) THEN BEGIN
            PRINT,'Plotting coadded 1D echelle spectra (black) and error (red)'
            good = WHERE(coad_echelle.dq EQ 0,n_good)
            IF n_good NE 0 THEN BEGIN
                x_splot,coad_echelle.wavelength[good],coad_echelle.flux[good],$
                  ytwo=coad_echelle.error[good],/block
                PRINT,'Plotting coadded 1D echelle reference spectrum (black)',$
                  ' and coadded spectra (red)'
                ref_good = WHERE(ref_spec_dq EQ 0,n_ref_good)
                IF n_ref_good NE 0 THEN $
                  x_splot,ref_spec_wave[ref_good],ref_spec_flux[ref_good],$
                  xtwo=coad_echelle.wavelength[good],ytwo=coad_echelle.flux[good],$
                  /block $
                ELSE PRINT,'No plot because no good data in reference'
            ENDIF ELSE PRINT,'No plot because no good data in coadded'
            history = [history,'Coadded 1D spectrum from echelle spectra '+$
                       'verified visually']
        ENDIF
        ;;For comparisons/debugging, will make S/N array for coadded
        ;;echelle spectrum 
        coad_echelle_snr = REPLICATE(-1,n_flux_good)
        coad_echelle_good = WHERE(coad_echelle.error NE 0D,n_coad_echelle_good)
        IF n_coad_echelle_good NE 0 THEN coad_echelle_snr[coad_echelle_good]= $
          coad_echelle.flux[coad_echelle_good]/$
          coad_echelle.error[coad_echelle_good]
        history = [history,'Orders coadded into 1D spectrum']

        
;;;;;;;;;;;;;;;
; Write 1D to FITS (template: esi_echcoaddfin.pro)
;;;;;;;;;;;;;;;
        IF nfits NE 0 THEN BEGIN
            ;;Need information from zeroeth extension
            tmp = mrdfits(files[REF_SPEC_INDX],0,head,/silent)
            sxaddpar,head,'COMMENT',' '
            sxaddpar,head,'COMMENT',' *** Added Keywords *** '
            sxaddpar,head,'COMMENT',' '
            sxaddpar, head, 'CRVAL1', ALOG10(MIN(coad_echelle.wavelength))
            sxaddpar, head, 'CDELT1', CDELT1
            sxaddpar, head, 'CRPIX1', 1   
            sxaddpar, head, 'CTYPE1', 'LINEAR' 
            sxaddpar, head, 'DC-FLAG', 1  
            sxaddpar, head, 'BITPIX', -32 
            sxaddpar, head, 'TOT_EXP', TOT_EXP,$
              'total exposure time (seconds)--calculated'
            sxaddhist,history,head
            CASE nfits OF 
                2: BEGIN
                    ;;flux FITS
                    outfile_echelle = OUTFILE+'_F.fits'
                    sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
                      files[REF_SPEC_INDX]
                    mwrfits,coad_echelle.flux,outfile_echelle,head,$
                      /create,/silent
                    PRINT,'NOTE: ',outfile_echelle,' created'
                    ;;error FITS
                    outfile_echelle = OUTFILE+'_E.fits'
                    sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
                      files[REF_SPEC_INDX]
                    tmp_err = coad_echelle.error 
                    err_bad = WHERE(coad_echelle.dq NE 0,n_err_bad)
                    IF n_err_bad NE 0L THEN $
                      coad_echelle.error[err_bad] = -1.
                    mwrfits,coad_echelle.error,outfile_echelle,head,$
                      /create,/silent
                    fxhmodify,outfile_echelle,'HISTORY',$
                      'Error = -1 where DQ non-zero'
                    PRINT,'NOTE: ',outfile_echelle,' created'
                    coad_echelle.error = tmp_err
                END
                3: BEGIN
                    ;;flux, error, DQ FITS
                    outfile_echelle = OUTFILE+'.fits'
                    sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
                      files[REF_SPEC_INDX]
                    mwrfits,coad_echelle.flux,outfile_echelle,head,$
                      /create,/silent
                    mwrfits,coad_echelle.error,outfile_echelle,/silent
                    mwrfits,coad_echelle.dq,outfile_echelle,/silent
                    PRINT,'NOTE: ',outfile_echelle,' created'
                    ;;flux FITS
                    outfile_echelle = OUTFILE+'_F.fits'
                    sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
                      files[REF_SPEC_INDX]
                    mwrfits,coad_echelle.flux,outfile_echelle,head,$
                      /create,/silent
                    PRINT,'NOTE: ',outfile_echelle,' created'
                    ;;error FITS
                    outfile_echelle = OUTFILE+'_E.fits'
                    sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
                      files[REF_SPEC_INDX]
                    tmp_err = coad_echelle.error 
                    err_bad = WHERE(coad_echelle.dq NE 0,n_err_bad)
                    IF n_err_bad NE 0L THEN $
                      coad_echelle.error[err_bad] = -1.
                    mwrfits,coad_echelle.error,outfile_echelle,head,$
                      /create,/silent
                    fxhmodify,outfile_echelle,'HISTORY',$
                      'Error = -1 where DQ non-zero'
                    PRINT,'NOTE: ',outfile_echelle,' created'
                    coad_echelle.error = tmp_err
                END
                ELSE: BEGIN 
                    outfile_echelle = OUTFILE+'.fits'
                    sxaddpar,head,'FILENAME',outfile_echelle,' previous: '+$
                      files[REF_SPEC_INDX]
                    mwrfits,coad_echelle.flux,outfile_echelle,head,$
                      /create,/silent
                    mwrfits,coad_echelle.error,outfile_echelle,/silent
                    mwrfits,coad_echelle.wavelength,outfile_echelle,/silent
                    mwrfits,coad_echelle.dq,outfile_echelle,/silent
                    PRINT,'NOTE: ',outfile_echelle,' created'
                END
            ENDCASE 
        ENDIF                   ;end nfits option 
    ENDIF                       ;end coadd orders to 1D option
ENDIF ELSE BEGIN                ;end echelle spectra option



;;;;;;;;;;;;;;;
; COADD 1D SPECTRA
;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;
; Manipulate DQ so that it's 0's (good) and 1's (bad)
;;;;;;;;;;;;;;;
    FOR ii=0,N_SPEC-1 DO BEGIN
        mask = makemask((*quasar[ii]).dq)
        (*quasar[ii]).dq = mask
    ENDFOR
    history = [history,'Acceptable STIS DQ include any combination of '+$
               '0, 16, 32, and 1024']


;;;;;;;;;;;;;;;
; Check wavelength coverages overlap
;;;;;;;;;;;;;;;
    problem_spec = -1
    wv_lim = [min((*quasar[0]).wavelength,max=mx),mx]
    FOR ii=1,N_SPEC-1 DO BEGIN
        tmp = [min((*quasar[ii]).wavelength,max=mx),mx]
        ;;If ii-th range completely outside of first spectrum
        if tmp[1] lt wv_lim[0] or tmp[0] gt wv_lim[1] then begin
            PRINT,'WARNING: Wavelength coverage non-overlapping in ',$
              files[ii],', exiting'
;            GOTO,QUIT
            STOP
            if problem_spec[0] eq -1 then problem_spec = [ii] $
            else problem_spec = [problem_spec,ii]
        endif 
    ENDFOR 

    ;;Assuming only two pointings in set so only two different
    ;;wavelength ranges
    if size(problem_spec,/dimensions) ne 0 then begin
        print,'NOTE: Limiting input spectra to those that overlap'
        tmp = lindgen(n_spec)
        tmp[problem_spec] = -1
        gd = where(tmp ne -1,n_spec)
        files = files[gd]
        quasar = quasar[gd]
        mask = mask[gd]
        history = [history,'Not including (wavelength out of bounds): '+$
                   files[problem_spec]]
    endif 


;;;;;;;;;;;;;;;
; Determine reference spectrum (highest S/N)
;;;;;;;;;;;;;;;
    N_PIX = (*quasar[0]).nelem
    ;;Median S/N for each pixel of each spectra and median S/N for
    ;;entire spectrum (stored in [*,N_PIX])
    sn_ratio = REPLICATE(-9.99,N_SPEC,N_PIX+1)
    ;;Median lambda relation for each order of each spectra and median
    ;;dlambda (i.e. (lambda_i+1 - lambda_i)) for entire spectrum
    lambda_rel = REPLICATE(-9.99,N_SPEC,N_PIX+1)
    FOR ii=0,N_SPEC-1 DO BEGIN
        lambda_rel[ii,0:N_PIX-2] = $
          ((SHIFT((*quasar[ii]).wavelength,-1))[0:N_PIX-2]-$
           (*quasar[ii]).wavelength[0:N_PIX-2]) ;/$ ;if log relation
                                ;(*quasar[ii]).wavelength[0:N_PIX-2]
        lambda_rel[ii,N_PIX-1] = lambda_rel[ii,N_PIX-2]
        good = WHERE((*quasar[ii]).dq EQ 0 AND $
                     (*quasar[ii]).error NE 0D,n_good)
        IF n_good NE 0 THEN sn_ratio[ii,good] = $
          (*quasar[ii]).flux[good]/(*quasar[ii]).error[good]
    ENDFOR 
    ;;Compare spectra to determine absolute wavelength range and
    ;;spectrum with highest S/N
    REF_SPEC_INDX = 0
    FOR ii=0,N_SPEC-1 DO BEGIN
        snr_good = WHERE(sn_ratio[ii,0:N_PIX-1] NE -9.99,n_snr_good)
        IF n_snr_good GE N_PIX_MIN THEN BEGIN
            sn_ratio[ii,N_PIX] = MEDIAN(sn_ratio[ii,snr_good],/EVEN)
            lambda_rel[ii,N_PIX] = MEDIAN(lambda_rel[ii,snr_good],/EVEN)
        ENDIF 
        IF sn_ratio[ii,N_PIX] GT sn_ratio[REF_SPEC_INDX,N_PIX] THEN $
          REF_SPEC_INDX = ii
    ENDFOR 
    LAMBDA_MIN = MIN((*quasar[REF_SPEC_INDX]).wavelength,MAX=LAMBDA_MAX)
    history = [history,'Reference spectrum: '+files[REF_SPEC_INDX]+$
               ' (highest S/N)']
    ;;For comparisons/debugging, will make S/N array of original
    ;;reference spectrum 
    ref_spec_snr = REPLICATE(-1,N_PIX)
    ref_spec_good = WHERE((*quasar[REF_SPEC_INDX]).error NE 0D,n_ref_spec_good)
    IF n_ref_spec_good NE 0 THEN ref_spec_snr[ref_spec_good] = $
      (*quasar[REF_SPEC_INDX]).flux[ref_spec_good]/$
      (*quasar[REF_SPEC_INDX]).error[ref_spec_good]


;;;;;;;;;;;;;;;
; Generate new wavelength array for re-binning
;;;;;;;;;;;;;;;
    dlambda = MEDIAN(lambda_rel[*,N_PIX],/EVEN)
    CDELT1 = dlambda
    GLOBAL_PIX = ceil((lambda_max-lambda_min)/cdelt1)
    global_wave = LAMBDA_MIN + DINDGEN(GLOBAL_PIX)*CDELT1


;;;;;;;;;;;;;;;
; Rebin spectra 
;;;;;;;;;;;;;;;
    rebin_spec = {wavelength:PTRARR(N_SPEC,/NOZERO),$
                  flux:PTRARR(N_SPEC,/NOZERO),$
                  error:PTRARR(N_SPEC,/NOZERO),dq:PTRARR(N_SPEC,/NOZERO)}
    FOR ii=0,N_SPEC-1 DO BEGIN
        var = ((*quasar[ii]).error)^2
        var_bad = WHERE((*quasar[ii]).dq NE 0,n_var_bad)
        IF n_var_bad NE 0 THEN var[var_bad] = 0.
        x_specrebin,(*quasar[ii]).wavelength,(*quasar[ii]).flux,$
          global_wave,new_flx,VAR=var,NWVAR=new_var,/SILENT
        new_dq = REPLICATE(0B,N_ELEMENTS(new_var))
        new_dq_bad = WHERE(new_var LE 0,n_new_dq_bad)
        CASE n_new_dq_bad OF
            0:                  ;continue
            1: BEGIN
                
                CASE new_dq_bad[0] OF
                    0: new_dq[new_dq_bad+1] = 25B
                    N_ELEMENTS(new_var): new_dq[new_dq_bad-1] = 25B
                    ELSE:
                ENDCASE 
                ;;DQ if new_var EQ 0 or -1
                new_dq[new_dq_bad] = 17B
                ;;Change new_var = -1 to 0
                new_var[new_dq_bad] = 0.
            END 
            ELSE: BEGIN
                ;;Neighbors of bad pixels bad
                new_dq[(new_dq_bad-1)[1:n_new_dq_bad-1]] = 25B
                new_dq[(new_dq_bad+1)[0:n_new_dq_bad-2]] = 25B
                ;;DQ if new_var EQ 0 or -1
                new_dq[new_dq_bad] = 17B
                ;;Change new_var = -1 to 0
                new_var[new_dq_bad] = 0.
            END 
        ENDCASE 
        rebin_spec.wavelength[ii] = PTR_NEW(global_wave)
        rebin_spec.flux[ii] = PTR_NEW(new_flx,/NO_COPY)
        rebin_spec.error[ii] = PTR_NEW(SQRT(new_var),/NO_COPY)
        rebin_spec.dq[ii] = PTR_NEW(new_dq,/NO_COPY)
    ENDFOR 
    history = [history,'Spectra rebinned to same wavelength scale']
    history = [history,'DQ = 17 if rebinned spectrum had variance <= 0']
    history = [history,'DQ = 25 for edges of variance <= 0']

;;;;;;;;;;;;;;;
; Determine flux scaling, based on S/N
;;;;;;;;;;;;;;;
    flux_scale = REPLICATE(-9.99,N_SPEC,GLOBAL_PIX+1)
    FOR ii=0,N_SPEC-1 DO BEGIN
        ;;scale on: (a) good data (ii-th and reference spectrum with
        ;;DQ = 0 and ii-th spectrum new_var > 0); (b) prevent dividing
        ;;by zero; and (c) only scale on areas with S/N GE TOLERANCE
        good = WHERE((*rebin_spec.dq[ii] EQ 0) AND $
                     ((*rebin_spec.flux[REF_SPEC_INDX]) NE 0D) $
                     AND ((*rebin_spec.error[ii]) NE 0D) AND $
                     ((*rebin_spec.flux[ii])/$
                      (*rebin_spec.error[ii]) GE TOLERANCE),n_good)
        IF n_good NE 0 THEN flux_scale[ii,good] = $
          (*rebin_spec.flux[ref_spec_indx])[good]/$
          (*rebin_spec.flux[ii])[good]
    ENDFOR
    FOR ii=0,N_SPEC-1 DO BEGIN
        good = WHERE(flux_scale[ii,0:GLOBAL_PIX-1] NE -9.99,n_good)
        IF n_good GE N_PIX_MIN THEN $
          flux_scale[ii,GLOBAL_PIX] = MEDIAN(flux_scale[ii,good],/EVEN)
    ENDFOR 
    history = [history,'Flux scale and S/N for each spectra:']
    history = [history,files+STRING(flux_scale[*,GLOBAL_PIX])+$
               STRING(sn_ratio[*,N_PIX])]


;;;;;;;;;;;;;;;
; Print wavelength relation, flux ratio/scale, and S/N
;;;;;;;;;;;;;;;
    IF KEYWORD_SET(print) THEN BEGIN
        OPENW,15,OUTFILE+'.wave' 
        OPENW,17,OUTFILE+'.snr'
        FOR kk=0,N_PIX DO BEGIN
            FOR ii=0,N_SPEC-1 DO BEGIN
                CASE ii OF
                    N_SPEC-1: BEGIN
                        PRINTF,15,FORMAT='(D11.6)',lambda_rel[ii,kk]
                        PRINTF,17,FORMAT='(D11.6)',sn_ratio[ii,kk]
                    END
                    0: BEGIN
                        PRINTF,15,FORMAT='(I4,TR3,D11.6,TR3,$)',kk+1,$
                          lambda_rel[ii,kk]
                        PRINTF,17,FORMAT='(I4,TR3,D11.6,TR3,$)',kk+1,$
                          sn_ratio[ii,kk]
                    END
                    ELSE: BEGIN
                        PRINTF,15,FORMAT='(D11.6,TR3,$)',lambda_rel[ii,kk]
                        PRINTF,17,FORMAT='(D11.6,TR3,$)',sn_ratio[ii,kk]
                    END
                ENDCASE
            ENDFOR
        ENDFOR
        CLOSE,15
        history = [history,'NOTE: '+OUTFILE+'.wave created']
        CLOSE,17
        history = [history,'NOTE: '+OUTFILE+'.snr created']

        OPENW,16,OUTFILE+'.flux'
        FOR kk=0,GLOBAL_PIX DO BEGIN
            FOR ii=0,N_SPEC-1 DO BEGIN
                CASE ii OF
                    N_SPEC-1: PRINTF,16,FORMAT='(D11.6)',flux_scale[ii,kk]
                    0: PRINTF,16,FORMAT='(I4,TR3,D11.6,TR3,$)',kk+1,$
                      flux_scale[ii,kk]
                    ELSE: PRINTF,16,FORMAT='(D11.6,TR3,$)',flux_scale[ii,kk]
                ENDCASE
            ENDFOR
        ENDFOR
        CLOSE,16
        history = [history,'NOTE: '+OUTFILE+'.flux created']
    ENDIF 
    

;;;;;;;;;;;;;;;
; Coadd spectra
;;;;;;;;;;;;;;;
    flux_1d = DBLARR(GLOBAL_PIX,N_SPEC,/NOZERO)
    var_1d = DBLARR(GLOBAL_PIX,N_SPEC,/NOZERO) 
    scale_1d = flux_scale[*,GLOBAL_PIX]
                                ;snr_1d = replicate(-1.,N_SPEC)
    snr_1d = sn_ratio[*,N_PIX]
    FOR ii=0,N_SPEC-1 DO BEGIN
        flux_1d[*,ii] = (*rebin_spec.flux[ii])
        var_1d[*,ii] = (*rebin_spec.error[ii])^2
        err_bad = WHERE(*rebin_spec.dq[ii] NE 0,n_err_bad)
        IF n_err_bad NE 0 THEN var_1d[err_bad,ii] = 0.
    ENDFOR 
    x_combspec,flux_1d,var_1d,flux_comb,var_comb,SCALE=scale_1d,SNR=snr_1d
    dq_comb = REPLICATE(0B,GLOBAL_PIX)
    var_bad = WHERE(var_comb LE 0D,n_var_bad)
    IF n_var_bad NE 0 THEN BEGIN
        ;;change var_comb = -1 to 0
        var_comb[var_bad] = 0.
        ;;DQ if var_comb EQ 0 or -1
        dq_comb[var_bad] = 19B
    ENDIF 

    ;;Want to truncate leading/trailing zeroes (because new_spec wavelength
    ;;GT original wavelength)
    ;;Can't do this 'spec_good' stuff because will assign features to
    ;;wrong wavelengths when 'spec_good' conditions not met in spectrum
    trim_array,flux_comb,lead_zero,n_flux_good 
    flux_good = lead_zero + INDGEN(n_flux_good)
    flux_comb = flux_comb[flux_good]
    wave_comb = global_wave[flux_good]
    error_comb = SQRT(var_comb[flux_good])
    dq_comb = dq_comb[flux_good]
    coad_1d = {wavelength:wave_comb,flux:flux_comb,error:error_comb,$
               dq:dq_comb}
    IF KEYWORD_SET(view) THEN BEGIN
        PRINT,'Plotting combined spectrum and error'
        good = WHERE(coad_1d.dq EQ 0,n_good)
        IF n_good NE 0 THEN BEGIN
            x_splot,coad_1d.wavelength[good],coad_1d.flux[good],$
              ytwo=coad_1d.error[good],/block
            PRINT,'Plotting rebinned reference spectrum (black) and ',$
              'combined spectrum (red)'
            ref_good = WHERE(*rebin_spec.dq[REF_SPEC_INDX] EQ 0,n_ref_good)
            IF n_ref_good NE 0 THEN $
              x_splot,(*rebin_spec.wavelength[REF_SPEC_INDX])[ref_good],$
              (*rebin_spec.flux[REF_SPEC_INDX])[ref_good],$
              xtwo=coad_1d.wavelength[good],ytwo=coad_1d.flux[good],/block $
            ELSE PRINT,'No plot because no good data in reference'
        ENDIF ELSE PRINT,'No plot because no good data in coadded'
        history = [history,'Coadded 1D spectrum verified visually']
    ENDIF 
    ;;For comparisons/debugging, will make S/N array for coadded 1D
    ;;spectrum 
    coad_1d_snr = REPLICATE(-1,n_flux_good)
    coad_1d_good = WHERE(coad_1d.error NE 0D,n_coad_1d_good)
    IF n_coad_1d_good NE 0 THEN coad_1d_snr[coad_1d_good] = $
      coad_1d.flux[coad_1d_good]/coad_1d.error[coad_1d_good]
    history = [history,'Truncated FLUX to eliminate leading/trailing zeroes']
    history = [history,'Spectra coadded']
    history = [history,'DQ = 19 if coadded variance <= 0']
    

;;;;;;;;;;;;;;;
; Write spectrum to FITS (template: esi_echcoaddfin.pro)
;;;;;;;;;;;;;;;
    IF nfits NE 0 THEN BEGIN
        ;;Need information from zeroeth extension
        tmp = mrdfits(files[REF_SPEC_INDX],0,head,/silent)    
        sxaddpar,head,'COMMENT',' '
        sxaddpar,head,'COMMENT',' *** Added Keywords *** '
        sxaddpar,head,'COMMENT',' '
        sxaddpar,head, 'TOT_EXP', TOT_EXP,$
          'total exposure time (seconds)--calculated'
        sxaddpar,head,'CRVAL1', MIN(coad_1d.wavelength)
        sxaddpar,head,'CDELT1', CDELT1,$
          'Median delta(lambda) of original spectra' 
        sxaddpar,head,'CRPIX1', 1
        sxaddpar,head,'CTYPE1', 'LINEAR' 
        sxaddpar,head,'DC-FLAG', 0 
        sxaddpar,head,'BITPIX', -32 
        sxaddhist,history,head
        CASE nfits OF
            2: BEGIN
                ;;flux FITS
                outfile_1d = OUTFILE+'_F.fits'
                sxaddpar,head,'FILENAME',outfile_1d,' previous: '+$
                  files[REF_SPEC_INDX]
                mwrfits, coad_1d.flux, outfile_1d, head, /create, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
                ;;error FITS
                outfile_1d = OUTFILE+'_E.fits'
                sxaddpar,head,'FILENAME',outfile_1d,' previous: '+$
                  files[REF_SPEC_INDX]
                tmp_err = coad_1d.error 
                err_bad = WHERE(coad_1d.dq NE 0,n_err_bad)
                IF n_err_bad NE 0L THEN $
                  coad_1d.error[err_bad] = -1.
                mwrfits, coad_1d.error, outfile_1d, head, /create, /silent
                fxhmodify,outfile_1d,'HISTORY',$
                  'Error = -1 where DQ non-zero'
                PRINT,'NOTE: ',outfile_1d,' created'
                coad_1d.error = tmp_err
            END
            3: BEGIN
                ;;flux, error, DQ FITS
                outfile_1d = OUTFILE+'.fits'
                sxaddpar,head,'FILENAME',outfile_1d,' previous: '+$
                  files[REF_SPEC_INDX]
                mwrfits, coad_1d.flux, outfile_1d, head, /create, /silent
                mwrfits, coad_1d.error, outfile_1d, /silent
                mwrfits, coad_1d.dq, outfile_1d, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
                ;;flux FITS
                outfile_1d = OUTFILE+'_F.fits'
                sxaddpar,head,'FILENAME',outfile_1d,' previous: '+$
                  files[REF_SPEC_INDX]
                mwrfits, coad_1d.flux, outfile_1d, head, /create, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
                ;;error FITS
                outfile_1d = OUTFILE+'_E.fits'
                sxaddpar,head,'FILENAME',outfile_1d,' previous: '+$
                  files[REF_SPEC_INDX]
                tmp_err = coad_1d.error 
                err_bad = WHERE(coad_1d.dq NE 0,n_err_bad)
                IF n_err_bad NE 0L THEN $
                  coad_1d.error[err_bad] = -1.
                mwrfits, coad_1d.error, outfile_1d, head, /create, /silent
                fxhmodify,outfile_1d,'HISTORY',$
                  'Error = -1 where DQ non-zero'
                PRINT,'NOTE: ',outfile_1d,' created'
                coad_1d.error = tmp_err
            END
            ELSE: BEGIN
                ;;flux, error, DQ FITS
                outfile_1d = OUTFILE+'.fits'
                sxaddpar,head,'FILENAME',outfile_1d,' previous: '+$
                  files[REF_SPEC_INDX]
                mwrfits, coad_1d.flux, outfile_1d, head, /create, /silent
                mwrfits, coad_1d.error, outfile_1d, /silent
                mwrfits, coad_1d.wavelength, outfile_1d, /silent
                mwrfits, coad_1d.dq, outfile_1d, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
            END
        ENDCASE
    ENDIF                       ;end nfits option
ENDELSE                         ;end 1D spectra option



gnrc:                           ;/generic set
;;;;;;;;;;;;;;;;;;;;;;;
;; COADD GENERIC SPECTRA
;;;;;;;;;;;;;;;;;;;;;;;
if keyword_set(generic) then begin
    if size(infile,/dimension) eq 0 then begin
        ;; Read in either one or two columns
        if not keyword_set(zqso) then $
          readcol,infile,files,format='(a)'  $
        else readcol,infile,files,zqso,format='(a,f)'
    endif else begin 
        ;; infile is already string array and presumes zqso is too
        files = infile
    endelse 

    n_spec = n_elements(files)

    ;;Conform input spectra to previous style
    data = {wavelength:ptrarr(n_spec),$
            flux:ptrarr(n_spec),$
            error:ptrarr(n_spec),$
            mask:ptrarr(n_spec) $
           }
    
    cdelt1 = fltarr(n_spec)
    crval1 = fltarr(n_spec)
    lambda_min = 1e6
    lambda_max = 0.
    tot_exp = 0.
    snr = fltarr(n_spec)
    scale = fltarr(n_spec)
    for ii=0,n_spec-1 do begin
        ;;Assume the error files are thus named
        tmp = files[ii]
        strput,tmp,'e.fits',strpos(tmp,'f.fits')
        
        fx = x_readspec(files[ii],0,head=hdr,wav=wv,sig=err,fil_sig=tmp)
        msk = replicate(0,n_elements(wv))
        if keyword_set(zqso) then wv = wv/(1+zqso[ii])
        gd = where(fx ne 0.,ngd)
        if ngd ne 0 then msk[gd] = 1 $
        else stop,'coadSTIS: no good data in'+files[ii]
        
        dum = measureSN(files[ii],err=tmp,/silent)
        snr[ii] = dum[2]
        scale[ii] = dum[0]
        
        cdelt1[ii] = fxpar(hdr,'CDELT1') 
        crval1[ii] = fxpar(hdr,'CRVAL1') 
        tot_exp = tot_exp+fxpar(hdr,'TOT_EXP')
        lambda_min = min(wv,max=tmp) < lambda_min
        lambda_max = tmp > lambda_max

        data.wavelength[ii] = ptr_new(wv,/no_copy)
        data.flux[ii] = ptr_new(fx,/no_copy)
        data.error[ii] = ptr_new(err,/no_copy)
        data.mask[ii] = ptr_new(msk,/no_copy)

    endfor 

    ;; Create new error array
    global_pix = ceil(alog10(lambda_max/lambda_min)/median(cdelt1))+1
    cdelt1 = median(cdelt1)
    crval1 = alog10(lambda_min)
    global_wave = 10^(crval1 + dindgen(global_pix)*cdelt1)
    
    ;;Will hold rebinned spectra in format read for x_combspec
    flux_comb = fltarr(global_pix,n_spec)
    error_comb = fltarr(global_pix,n_spec)
    mask_comb = fltarr(global_pix,n_spec)
    var_comb = fltarr(global_pix,n_spec)
    
    for ii=0,n_spec-1 do begin
        var = (*data.error[ii])^2
        var_bad = where((*data.mask[ii]) eq 0,n_var_bad)
        if n_var_bad ne 0 then var[var_bad] = 0.
        x_specrebin,(*data.wavelength[ii]),(*data.flux[ii]),$
          global_wave,new_flx,var=var,nwvar=new_var,/silent

        new_msk = replicate(1,global_pix)
        new_msk_bad = where(new_var le 0,n_new_msk_bad)
        if n_new_msk_bad ne 0 then begin
            ;;Neighbors of bad pixels bad
            new_msk[(new_msk_bad-1)[1:n_new_msk_bad-1]] = 0
            new_msk[(new_msk_bad+1)[0:n_new_msk_bad-2]] = 0
            ;;MSK if new_var EQ 0 or -1
            new_msk[new_msk_bad] = 0
            ;;Change new_var = -1 to 0
            new_var[new_msk_bad] = 0.
        endif
        flux_comb[*,ii] = new_flx
        error_comb[*,ii] = SQRT(new_var)
        mask_comb[*,ii] = new_msk
        var_comb[*,ii] = new_var
                                ;stop
    ENDFOR

    tmp = max(snr,ref_spec_indx)
    history = [history,'Reference spectrum: '+files[REF_SPEC_INDX]+$
               ' (highest S/N)']
    history = [history,'Header from reference spectrum ']

    history = [history,'All orders rebinned to same '+$
               'wavelength scale with x_specrebin.pro']

    ;;Scale to highest S/N spectrum
    scale = scale[ref_spec_indx]/scale
    x_combspec,flux_comb,var_comb,new_flx,new_var,scale=scale,snr=snr

    gd = where(new_var ne -1,ngd,complement=bd,ncomplement=nbd)
    new_msk = replicate(0,global_pix)
    if ngd ne 0 then new_msk[gd] = 1
    new_err = replicate(-1.,global_pix)
    new_err[gd] = sqrt(new_var[gd])

    history = [history,'Flux scale for each spectra:']
    history = [history,files+STRING(scale)]
    history = [history,'Each spectrum weighted by its median S/N']


;;;;;;;;;;;;;;;
; Write spectrum to FITS (template: esi_echcoaddfin.pro)
;;;;;;;;;;;;;;;
    IF nfits NE 0 THEN BEGIN
        ;;Need information from zeroeth extension
        tmp = mrdfits(files[REF_SPEC_INDX],0,head,/silent)    
        sxaddpar,head,'COMMENT',' '
        sxaddpar,head,'COMMENT',' *** Added Keywords *** '
        sxaddpar,head,'COMMENT',' '
        sxaddpar,head, 'TOT_EXP', TOT_EXP,$
          'total exposure time (seconds)--calculated'
        sxaddpar,head,'CRVAL1', crval1
        sxaddpar,head,'CDELT1', CDELT1,$
          'Median CDELT1 of original spectra' 
        sxaddpar,head,'CRPIX1', 1
        sxaddpar,head,'CTYPE1', 'LINEAR' 
        sxaddpar,head,'DC-FLAG', 1
        sxaddpar,head,'BITPIX', -32 
        sxaddhist,history,head
        CASE nfits OF
            2: BEGIN
                ;;flux FITS
                outfile_1d = OUTFILE+'_F.fits'
                mwrfits, new_flx, outfile_1d, head, /create, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
                ;;error FITS
                outfile_1d = OUTFILE+'_E.fits'
                tmp_err = new_err 
                err_bad = WHERE(new_msk NE 1,n_err_bad)
                IF n_err_bad NE 0L THEN $
                  new_err[err_bad] = -1.
                mwrfits, new_err, outfile_1d, head, /create, /silent
                fxhmodify,outfile_1d,'HISTORY',$
                  'Error = -1 where mask eq 0'
                PRINT,'NOTE: ',outfile_1d,' created'
                new_err = tmp_err
            END
            3: BEGIN
                ;;flux, error, DQ FITS
                outfile_1d = OUTFILE+'.fits'
                mwrfits, new_flx, outfile_1d, head, /create, /silent
                mwrfits, new_err, outfile_1d, /silent
                mwrfits, new_msk, outfile_1d, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
                ;;flux FITS
                outfile_1d = OUTFILE+'_F.fits'
                mwrfits, new_flx, outfile_1d, head, /create, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
                ;;error FITS
                outfile_1d = OUTFILE+'_E.fits'
                tmp_err = new_err 
                err_bad = WHERE(new_msk NE 1,n_err_bad)
                IF n_err_bad NE 0L THEN $
                  new_err[err_bad] = -1.
                mwrfits, new_err, outfile_1d, head, /create, /silent
                fxhmodify,outfile_1d,'HISTORY',$
                  'Error = -1 where mask ne 1'
                PRINT,'NOTE: ',outfile_1d,' created'
                new_err = tmp_err
            END
            ELSE: BEGIN
                ;;flux, error, DQ FITS
                outfile_1d = OUTFILE+'.fits'
                mwrfits, new_flx, outfile_1d, head, /create, /silent
                mwrfits, new_err, outfile_1d, /silent
                mwrfits, global_wave, outfile_1d, /silent
                mwrfits, new_msk, outfile_1d, /silent
                PRINT,'NOTE: ',outfile_1d,' created'
            END
        ENDCASE
    ENDIF                       ;end nfits option

    stop
endif                           ;end generic option



PRINT,'Finished coadding spectra listed in ',infile
QUIT:
PTR_FREE,quasar
PTR_FREE,rebin_spec.wavelength
PTR_FREE,rebin_spec.flux
PTR_FREE,rebin_spec.error
PTR_FREE,rebin_spec.dq
END
