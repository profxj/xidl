;+ 
; NAME:
; lris_zfind
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   lris_zfind, fspec_fil
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lris_zfind
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;    4-Jun-2008 Adopted wfccd_zfind, KLC
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro lris_zfind, specfil, chip, ZMIN=zmin, ZMAX=zmax, PLOT=plot, $
                NPOLY=npoly,inflg=inflg,zfil=zfil,fitstar=fitstar,$
                eigengal=eigengal, eigenstar=eigenstar, eigendir=eigendir, $
                zans=zans, _extra=extra
;
if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
      'lris_zfind, specfil, chip, ZMIN=, ZMAX=, /PLOT, NPOLY= [v1.0]'
    return
endif 
chip = strupcase(strtrim(chip,2))

;  Optional Keywords

if not keyword_set( npoly ) then npoly = 3L
if not keyword_set( zmin ) then zmin = -0.03
if not keyword_set( zmax ) then zmax = 0.5
if not keyword_set( SKYLIN ) then begin
    nskylin = 6L
    skylin = fltarr(2,10)
    skylin[*,0] = [5569.5, 5593.0]
    skylin[*,1] = [5885.3, 5905.0]
    skylin[*,2] = [6290.3, 6317.0]
    skylin[*,3] = [6351.4, 6380.5]
    skylin[*,4] = [6855.5, 6906.] ; A band
    skylin[*,5] = [7586.5, 7693.] ; B band
endif else begin
    sz_skylin = size(skylin, /dimensions)
    nskylin = sz_skylin[1]
endelse
;if not arg_present(inflg) then inflg = 2 ;LONGSLIT format

;; Select eigenspec
if not keyword_set(eigendir) then $
   eigendir = getenv('LONGSLIT_DIR')+'/pro/LRIS/'
case chip of 
    'BLUE': begin
        eigengal = 'lrisEigenGal-52223b400D560.fits'
        eigenstar = 'lrisEigenStar-52374b400D560.fits'
    end
    'RED': begin
        eigengal = 'lrisEigenGal-52223r600-10000.fits'
        eigenstar = 'lrisEigenStar-52374r600-10000.fits'
    end 
    else: begin
       if not keyword_set(eigengal) and not keyword_set(eigenstar) then $
          stop,'lris_zfind: eigengal and/or eigenstar not set'
    end
endcase 

;  Loop
nspec = n_elements(specfil)
for qq=0L,nspec-1 do begin
    splog, 'Examining '+strtrim(specfil[qq],2)

    ;; Read in 
    spec_fx = x_readspec(specfil[qq],inflg=inflg,wav=spec_wv,$
                         sig=spec_sig,npix=npix,_extra=extra) ; fil_sig=

    if chip eq 'BLUE' then begin
        ;; Make sure to go from blue to red
        srt = sort(spec_wv)
        spec_wv = spec_wv[srt]
        spec_fx = spec_fx[srt]
        spec_sig = spec_sig[srt]
    endif 
    spec_var = spec_sig^2


    ;; Trim
    gdvar = where(spec_sig GT 0., ngdvar,complement=bd,ncomplement=nbd)
    if nbd ne 0 then spec_var[bd] = -1.

    if(ngdvar gt 0) then begin

        mnpix = min(gdvar, MAX=mxpix)

        ;; Final arrays
        npix = mxpix-mnpix+1
        fin_wv = spec_wv[mnpix:mxpix]
        fin_fx = double(spec_fx[mnpix:mxpix])
        fin_var = spec_var[mnpix:mxpix]
        
        fin_ivar = dblarr(npix)
        
        posvar = where(fin_var GT 0.)
        fin_ivar[posvar] = 1./fin_var[posvar]
        
        ;; Kill Sky Lines
        for i=0L,nskylin-1 do begin
            skypix = where(fin_wv GE skylin[0,i] AND $
                           fin_wv LE skylin[1,i], nskypix)
            if nskypix NE 0 then fin_ivar[skypix] = 0.d
        endfor

        ;; Create the header
        mkhdr, head, 1, [5000L]
        sxaddpar, head, 'COEFF0', double(alog10(min(fin_wv)))
        dwv = fin_wv - shift(fin_wv,1)
        sxaddpar, head, 'COEFF1', alog10(1+median(dwv/fin_wv))
        
        ;; DO GALAXY REDSHIFTS 

        ;; Launch zfind
        npoly=3
        zans_gal = x_zfind(fin_fx, fin_ivar, hdr=head, $
                           eigenfile=eigengal,eigendir=eigendir,$
                           zmin=zmin, zmax=zmax, doplot=plot, npoly=npoly, $
                           /silent)
        splog, 'found z= '+string(zans_gal.z)
        
        ;; Check on zans_gal
        if zmax-zans_gal.z LT 0.02 OR zans_gal.z_err LT 0. then begin
            splog, 'zans_gal ~ zmax!!  Retrying..'
            zans_gal = x_zfind(fin_fx, fin_ivar, hdr=head, $
                               eigenfile=eigengal,eigendir=eigendir,$
                               zmin=zmin, zmax=zmax+0.1, DOPLOT=plot, $
                               npoly=(npoly+2L)<6L, /silent)
            splog, 'found z= '+string(zans_gal.z)
        endif

        splog, zans_gal.z, zans_gal.z_err, zans_gal.rchi2

        ;; Check again
        if (zmax-zans_gal.z)+0.1 LT 0.02 OR zans_gal.z_err LT 0. then begin
            splog, 'Still trouble in zfind!'
        endif


                                ; PLOT
        if keyword_set( PLOT ) then begin
            splog,'Chi^2 = ',zans_gal.rchi2
            synth = x_synthspec(zans_gal, loglam=alog10(fin_wv),$
                                eigendir=eigendir)
            x_splot, fin_wv, fin_fx, YTWO=synth, /block
        endif


        zans_gal.class='GALAXY'

        if keyword_set(fitstar) then begin
           ;; DO STELLAR REDSHIFTS (ADDED HERE BY MRB BASED ON
           ;; idlspec2d/pro/spec1d/spreduce1d.pro )
           
           ;; Select the stars eigen-file here to detemine how many
           ;; templates are in it 
           shdr = headfits(djs_filepath(eigenstar, root_dir=eigendir))
           nstar = sxpar(shdr, 'NAXIS2') > 1
           npoly = 3
           zmin = -0.0033       ; -1000 km/sec
           zmax = 0.0033        ; +1000 km/sec
           pspace = 1
           nfind = 1

            
            zans_star=0
            for istar=0, nstar-1 do begin
                subclass = strtrim( sxpar(shdr, 'NAME'+$
                                          strtrim(string(istar),2)), 2)
                plottitle = subclass + '-Star Redshift'
                
                splog, 'Compute STAR (' + subclass + ') redshifts:', $
                  ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
                t0 = systime(1)
                tmp_zans_star = x_zfind(fin_fx, fin_ivar, hdr=head, $
                                        eigenfile=eigenstar, columns=istar, $
                                        eigendir = eigendir, $
                                        npoly=npoly, zmin=zmin, zmax=zmax, $
                                        pspace=1, $
                                        nfind=nfind, width=5*pspace, $
                                        plottitle=plottitle, doplot=plot, $
                                        debug=debug, /silent)
                splog, 'found z= '+string(tmp_zans_star.z)
                splog, 'CPU time to compute STAR redshifts = ', systime(1)-t0
                

                ;; PLOT
;      if keyword_set( PLOT ) then begin
;          splog,'Chi^2 = ',tmp_zans_star.rchi2
;          synth = x_synthspec(tmp_zans_star, loglam=alog10(fin_wv),$
;            eigendir=eigendir)
;          x_splot, fin_wv, fin_fx, YTWO=synth, /block
;      endif

                tmp_zans_star.class = 'STAR'
                tmp_zans_star.subclass = subclass
                
                if(n_tags(zans_star) gt 0) then $
                  zans_star=[zans_star, tmp_zans_star] $
                else $
                  zans_star=tmp_zans_star

            endfor

            ;;  replace soln if zans_star is better
            isort = sort(zans_star.rchi2 + (zans_star.dof EQ 0)* $
                         max(zans_star.rchi2))
            if(zans_star[isort[0]].rchi2 lt zans_gal.rchi2) then $
              zans_gal=zans_star[isort[0]]
        endif                   ;/fitstar


        ;; Save
        ;; ??? mrb 02-19-04 deal with change in structure
        tmp_zans = zans_gal
        struct_assign, zans_gal, tmp_zans, /nozero
        inegone=where(tmp_zans.tcolumn eq -1L, nnegone)
        if(nnegone gt 0) then begin
            imaxnegone=max(inegone)
            if(imaxnegone lt n_elements(tmp_zans.tcolumn)-1L) then $
              tmp_zans.tcolumn[imaxnegone+1: $
                               n_elements(tmp_zans.tcolumn)-1L]=-1L
        endif 
    endif                       ;if good pixels to fit

    splog, 'Final redshift = ', tmp_zans.z,'+/-',tmp_zans.z_err,$
           tmp_zans.rchi2

    ;; Save to pass back with zans
    if qq eq 0 then zans = replicate(tmp_zans,nspec) $ ; will overwrite
    else zans[qq] = tmp_zans
endfor                          ;loop nspec

print,''
if keyword_set(zfil) then begin
    close,1
    openw,1,zfil
    printf,1,'Spectrum','z','zerr','rchi2',format='(a35,1x,a9,1x,a9,1x,a5)'
    writecol,zfil,specfil,zans.z,zans.z_err,$
      zans.rchi2,fmt='(a35,1x,f9.6,1x,f9.6,1x,f5.2)',filnum=1
    close,1
    splog, 'created ',zfil
endif

return
end
