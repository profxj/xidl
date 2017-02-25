;+ 
; NAME:
; wfccd_zfind
;    Version 1.0
;
; PURPOSE:
;   Converts SDSS eigenfunctions into WFCCD
;
; CALLING SEQUENCE:
;   
;   wfccd_zfind, fspec_fil
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   XSIZE      - Size of gui in screen x-pixels (default = 1000)
;   YSIZE      - Size of gui in screen y-pixels (default = 600)
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_zfind, wfccd, maskid, expsr
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   19-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro wfccd_zfind, fspec_fil, obj_nm, ZMIN=zmin, ZMAX=zmax, PLOT=plot, $
                 NPOLY=npoly, WVMNX=wvmnx

;
if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'wfccd_zfind, fspec_fil, [obj_nm], ZMIN=, ZMAX=, /PLOT, NPOLY= [v1.0]'
    return
endif 

;  Optional Keywords

if not keyword_set( npoly ) then npoly = 3L
if not keyword_set( zmin ) then zmin = -0.03
if not keyword_set( zmax ) then zmax = 0.5
if not keyword_set( WVMNX ) then wvmnx = [3800., 9000.]
if not keyword_set( SCL_VAR ) then scl_var = 1.0
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

;  Read in

wfccd_wrfspec, wffspec, fspec_fil, /read

;  Spec

if not keyword_set( OBJ_NM ) then gdsp = where(wffspec.flg_anly NE 0, ngd) $
else begin
    gdsp = lonarr(1)
    gdsp[0] = x_getobjnm(wffspec,obj_nm)
    ngd = 1L
endelse


;  Loop

nbad = 0L
obj_bad = strarr(1000)

for q=0L,ngd-1 do begin
    print, 'wfccd_zfind: Examining '+strtrim(wffspec[gdsp[q]].slit_id,2)+$
      strtrim(wffspec[gdsp[q]].obj_id,2)+';  Number '+strtrim(q,2)+$
      ' of '+strtrim(ngd-1,2)

    npix = wffspec[gdsp[q]].npix  

    spec_wv = wffspec[gdsp[q]].wave[0:npix-1] 
    spec_fx = wffspec[gdsp[q]].fx[0:npix-1]

                                ; Trim
    mnwv = min( abs(spec_wv-wvmnx[0]), imnwv)
    mxwv = min( abs(spec_wv-wvmnx[1]), imxwv)
    gdvar = where(wffspec[gdsp[q]].var[0:npix-1] GT 0., ngdvar)

    if(ngdvar gt 0) then begin

        mnvar = min(gdvar, MAX=mxvar)
        
        mnpix = imnwv > mnvar
        mxpix = imxwv < mxvar
        
        ;; Trim off last 3 pix (troubles with rebinning)
        mxpix = mxpix - 3L
        
        ;; Final arrays
        npix = mxpix-mnpix+1
        fin_wv = spec_wv[mnpix:mxpix]
        fin_fx = double(spec_fx[mnpix:mxpix])
        fin_var = wffspec[gdsp[q]].var[mnpix:mxpix]
        
        fin_ivar = dblarr(npix)
        
        posvar = where(fin_var GT 0.)
        fin_ivar[posvar] = 1./(fin_var[posvar]*scl_var)
        
        ;; Kill Sky Lines
        for i=0L,nskylin-1 do begin
            skypix = where(fin_wv GE skylin[0,i] AND $
                           fin_wv LE skylin[1,i], nskypix)
            if nskypix NE 0 then fin_ivar[skypix] = 0.d
        endfor
        
;      if keyword_set( PLOT ) then begin
;          clr = getcolor(/load)
;          !p.color = clr.black
;          !p.background = clr.white
;          posivar = where(fin_ivar GT 0.)
;          fin_sig = fltarr(npix)
;          fin_sig[posivar] = 1./sqrt(fin_ivar[posivar])
;          x_splot, fin_wv, fin_fx, YTWO=fin_sig, /block
;      endif
        
        ;; Create the header
        mkhdr, head, 1, [5000L]
        sxaddpar, head, 'COEFF0', double(alog10(fin_wv[0]))
        sxaddpar, head, 'COEFF1', 100.d/2.9979e5/alog(10.d)
        
; DO GALAXY REDSHIFTS 

        ;; Launch zfind
        npoly=3
;        zmin = -0.0033   ; zmin = -1000 km/s
;        zmax = 0.5
        zans = x_zfind(fin_fx, fin_ivar, hdr=head, $
                       eigenfile='wfEigenGal-52223.fits', $
                       eigendir=getenv('XIDL_DIR')+'Spec/WFCCD/Analysis/',$
                       zmin=zmin, zmax=zmax, doplot=plot, npoly=npoly, /silent)
        splog, 'found z= '+string(zans.z)
        
                                ; Check on zans
        if zmax-zans.z LT 0.02 OR zans.z_err LT 0. then begin
            print, 'wfccd_zfind: zans ~ zmax!!  Retrying..'
            zans = x_zfind(fin_fx, fin_ivar, hdr=head, $
                           eigenfile='wfEigenGal-52223.fits', $
                           eigendir=getenv('XIDL_DIR')+'Spec/WFCCD/Analysis/',$
                           zmin=zmin, zmax=zmax+0.1, DOPLOT=plot, $
                           npoly=(npoly+2L)<6L, /silent)
            splog, 'found z= '+string(zans.z)
        endif

        print, 'wfccd_zfind: ', zans.z, zans.z_err, zans.rchi2

                                ; Check again
        if (zmax-zans.z)+0.1 LT 0.02 OR zans.z_err LT 0. then begin
            print, 'wfccd_zfind: Still trouble in zfind!'
            nbad = nbad + 1
            obj_bad(nbad) = strtrim(wffspec[gdsp[q]].slit_id,2)+$
              strtrim(wffspec[gdsp[q]].obj_id,2)
        endif


       ; PLOT
        if keyword_set( PLOT ) then begin
            print,zans.rchi2
            synth = x_synthspec(zans, loglam=alog10(fin_wv))
            x_splot, fin_wv, fin_fx, YTWO=synth, /block
        endif



        zans.class='GALAXY'

; DO STELLAR REDSHIFTS (ADDED HERE BY MRB BASED ON
; idlspec2d/pro/spec1d/spreduce1d.pro )

; Select the stars eigen-file here to detemine how many templates are in it
        eigendir = concat_dir(getenv('XIDL_DIR'), '/Spec/WFCCD/Analysis/')
        allfiles = findfile(djs_filepath('wfEigenStar-52374.fits', $
                                         root_dir=eigendir), count=ct)
        if (ct EQ 0) then $
          message, 'Unable to find EIGENFILE matching '+eigenfile
        eigenfile = fileandpath(allfiles[ (reverse(sort(allfiles)))[0] ])
        shdr = headfits(djs_filepath(eigenfile, root_dir=eigendir))
        nstar = sxpar(shdr, 'NAXIS2') > 1
        npoly = 3
        zmin = -0.0033          ; -1000 km/sec
        zmax = 0.0033           ; +1000 km/sec
        pspace = 1
        nfind = 1
        
        zans_star=0
        for istar=0, nstar-1 do begin
            subclass = strtrim( sxpar(shdr, 'NAME'+strtrim(string(istar),2)), 2)
            plottitle = subclass + '-Star Redshift'
            
            splog, 'Compute STAR (' + subclass + ') redshifts:', $
              ' ZMIN=', zmin, ' ZMAX=', zmax, ' PSPACE=', pspace
            t0 = systime(1)
            tmp_zans_star = x_zfind(fin_fx, fin_ivar, hdr=head, $
                                    eigenfile=eigenfile, columns=istar, $
                                    eigendir = eigendir, $
                                  npoly=npoly, zmin=zmin, zmax=zmax, pspace=1, $
                                  nfind=nfind, width=5*pspace, $
                                  plottitle=plottitle, doplot=plot, $
                                  debug=debug, /silent)
            splog, 'found z= '+string(tmp_zans_star.z)
            splog, 'CPU time to compute STAR redshifts = ', systime(1)-t0
            

; PLOT
;      if keyword_set( PLOT ) then begin
;          print,zans.rchi2
;          synth = x_synthspec(tmp_zans_star, loglam=alog10(fin_wv))
;          x_splot, fin_wv, fin_fx, YTWO=synth, /block
;      endif

            tmp_zans_star.class = 'STAR'
            tmp_zans_star.subclass = subclass
            
            if(n_tags(zans_star) gt 0) then $
              zans_star=[zans_star, tmp_zans_star] $
            else $
              zans_star=tmp_zans_star

        endfor

;;     replace soln if zans_star is better
        isort = sort(zans_star.rchi2 + (zans_star.dof EQ 0)* $
                     max(zans_star.rchi2))
        if(zans_star[isort[0]].rchi2 lt zans.rchi2) then $
          zans=zans_star[isort[0]]

        ;; Save
        ;; ??? mrb 02-19-04 deal with change in structure
        tmp_zans=wffspec[gdsp[q]].zans
        struct_assign, zans, tmp_zans, /nozero
        inegone=where(tmp_zans.tcolumn eq -1L, nnegone)
        if(nnegone gt 0) then begin
            imaxnegone=max(inegone)
            if(imaxnegone lt n_elements(tmp_zans.tcolumn)-1L) then $
              tmp_zans.tcolumn[imaxnegone+1: $
                               n_elements(tmp_zans.tcolumn)-1L]=-1L
        endif
        wffspec[gdsp[q]].zans = tmp_zans
    endif

    print, 'wfccd_zfind: Final redshift = ', tmp_zans.z
endfor


;; bad fits
if nbad GT 0L then begin
    print, 'wfccd_zfind: Bad fits!  Check output!'
    print, obj_bad[1:nbad]
endif

;; Write back
wfccd_wrfspec, wffspec, fspec_fil

return
end
