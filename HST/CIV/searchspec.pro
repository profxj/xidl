;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; searchspec.pro
; Author: Kathy Cooksey     Date: 22 Sep 2005
; Project: O VI absorber survey and HST Archive Metal-Line
;          System Survey with Jason Prochaska
; Description: Searches given spectrum for plausible features,
;              based primarily on S/N and Doppler paramter
; Inputs:
;  wave -- wavelength array
;  flux -- flux array
;  errr -- error array
;  nsigma -- number of sigma used to screen possible lines
;  dopb -- sigma (width) for Gaussian distribution used to
;          search for lines (see Notes)
; Optional:
;  ndopb -- number of dopb used to determine region included in
;           EW measurement (default: 8)
;  name -- name of quasar for records
;  fitgauss -- attempt to fit gaussian or not (not recommended)
;  header -- variable name to hold returned header
;  ostrct -- output structure
;  chisqtol -- 2 element array constraining acceptable reduced
;              chisq for Gaussian fitting (default: [0.,1.2]
;  saturate -- cutoff fraction of saturated pixels to fit Gaussian
;  instr -- instrument ID number (default: instr = inflg) for linlst file
;           128 STIS
;            64 LiF 2A
;            32 LiF 1B
;            16 SiC 1A
;             8 LiF 1A
;             4 LiF 2B
;             2 SiC 2A
;             1 SiC 1B
;  dw0 -- shift for wavelength
;  w0 -- reference wavelength for stretching and shifting wavelength
;  strct_fil -- name of file to save search information 
;              (default: searchspec.fits)
;  table -- prints information for each possible feature 
;           and search parameters in header
;           (Centroid, bounds, dopb, ew and error, nsigma, chi^2,
;           weight type)
;  linlst -- create explicitly formatted line list
;  group -- option for how to group hits as features:
;            1: centroid falling within wavelength limits
;            2: wavelength limits falling within limits (most
;               liberal, minimize features)
;            3: consecutive (most conservative, maximize features)
;  /nochng -- however wavelength limits set with group option, don't
;             change them
;  /shift -- shift wavelength as specified by dw0 and w0
;  /silent -- don't print messages
;  /view -- show spectra with potential features marked
;  /comp -- measure EW and sigEW 3 different ways
;  /debug -- stop at useful places
; Output: 
;  strct_fil: structure of autocorrelation with
;         detections in ext. 1, all in ext. 2, (comp in ext. 3)
;  ostrct: or can be structure
;  header: header of structure
; Optional output:
;  table: table of information
; Notes:
;   Gaussian = 1/(sigmag*sqrt(2*pi))*exp(-0.5*((x-mu)/sigmag)^2)
;   Weight or Doppler fit = 1/(b*sqrt(pi))*exp(-((v-mu)/b)^2)
;   b = sqrt(2)*sigmag_v = sqrt(2)*c*sigmag_wv/wv_mu
;   x-mu = +/-nsigmag*sigmag = +/-ndopb*b
; Examples: 
;   searchspec,wv,fx,err,3,20,/view
; History:
;   19 Apr 05     created by KLC
;    3 May 05     change sigg to Doppler parameter b and 
;                 Gaussian to Doppler profile
;    4 May 05     added range keyword; write search parameters
;                 to output
;   17 May 05     Revamp with correct 3*sigma ranges and new
;                 synthetic test case; modified output structure
;   17 May 05     complete revmap of format
;   18 May 05     finished debugging new format; added /plot
;   22 Sep 05     searchspec -> searchspec2 to handle 1 spec
;   28 Sep 05     all EW previously calculated wrong; no longer
;                 divide by total(g*dwv) but total(g); changed 
;                 to be in velocity space  
;   30 Sep 05     searchspec2 -> searchspec3 with new method
;    7 Oct 05     mostly completed overhaul and works well
;    3 Nov 05     accept arrays, not name of spectrum, inflg,
;                 and fil_sig; add name (replace spec)
;    4 Nov 05     check that independent variables for Gaussian
;                 fit is more than parameters
;   15 Nov 05     added calls to search_wrtab, search_compew,
;                 plotauto; cleaned up code;
;                 add additional check for new wrest when fitting
;                 2nd Gaussian
;   13 Dec 05     saved "duplicates" because search_overlap not
;                 perfect
;   17 Jan 06     velcoity = c*(wave/refwv - 1.), etc.
;   20 Jan 06     velocity finally correct
;    7 Apr 06     set bad chi^2 = -9.99, set saturation cut-off
;   25 Jan 07     modify weighting and handling of flux
;   14 May 07     modify nsigma cut per pixel; dopb-->FWHM
;   19 Feb 08     restrict -|error| <= flux <= 1+error and
;                 exclude gaps
;    3 Nov 10     now use civ_search_overlap
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro searchspec, wave, flux, error, nsigma, b_fwhm, name=name,$
                fitgauss=fitgauss,header=header, ostrct=ostrct,$
                ndopb=ndopb, chisqtol=chisqtol,saturate=saturate,instr=instr,$
                dw0=dw0, w0=w0, strct_fil=strct_fil, table=table, linlst=linlst, $
                group=group,nochng=nochng,shift=shift,$
                each=each,silent=silent,view=view,debug=debug,comp=comp
if not keyword_set(strct_fil) and not keyword_set(ostrct) then $
   strct_fil = 'searchspec.fits'
if not keyword_set(ndopb) then ndopb = 3
if not keyword_set(chisqtol) then chisqtol = [0.,2.5]
if not keyword_set(saturate) then saturate = 0.15 ;15%
if not keyword_set(group) then group = 1
if not keyword_set(instr) then instr = 0

dopb = b_fwhm/2.354             ;Gaussian sigma 

;;Constants
sol = 2.998e5                   ; km/s
rt2 = sqrt(2.)                  
rtpi = sqrt(!pi)
;; Match the elements of fuselinstrct
detection = {ion:' ', $
             wrest:-1.d, $
             wv_lim:dblarr(2), $
             flg:0, $
             instr:instr, $
             ew:0.d, $
             sigew:0.d,  $
             Ncolm: 0., $
             sigNcolm: 0., $
             zabs:0.d, $
             zsig: 0.d, $
             doppler:0.d, $
             doperr:0.d, $
             amp:0., $
             amperr:0., $
             chisq:0.}
;;each detection flg = 0 for non-hit and 1 for hit
;;auto detection flg = 2 for boxcar and 4 for Gaussian-integrated EW

;;Shift and stretch wavelength
if keyword_set(shift) then wave = wave + dw0*w0/wave
dwv=wave-shift(wave,1)
dwv[0]=dwv[1]

;; Trim outliers in data
sverror = error
svflux = flux
bd = where(error le 0.,nbd,complement=gd,ncomplement=ngd)
if nbd ne 0 then begin
    error[bd] = 0.
    flux[bd] = -1.
endif
if ngd eq 0 then begin
    print,'searchspec: no data with sig > 0'
    ndet = 0
    auto = -1            ;empty
    each = -1
    nrej = 0
    print,'searchspec: No lines detected with S/N > ',nsigma
    if keyword_set(strct_fil) then begin 
       bdstrct_fil = strmid(strct_fil,0,strpos(strct_fil,'.',/reverse_search))+'_nogddata.fits'
       goto,no_gd_data          ;skip to write structure
    endif 
    ostrct = -1
endif
bd = where(flux[gd] lt -abs(error[gd]),nbd)
if nbd ne 0 then flux[gd[bd]] = -abs(error[gd[bd]]) ;floor
bd = where(flux[gd] gt 1+error[gd],nbd)
if nbd ne 0 then flux[gd[bd]] = 1+error[gd[bd]] ;ceil

npix = n_elements(flux)

dv = dwv/wave*sol
;velocity = dblarr(npix)
;for ii=0l,npix-1 do velocity[ii] = total(dv[0:ii])
velocity = total(dv,/cumulative)


;;Instantiate
dum = min(abs(velocity-velocity[0]-ndopb*dopb),imu) ;first centroid
dum = min(abs(velocity-velocity[npix-1]+ndopb*dopb),ilast) ;last cent
maxdet = ilast-imu+1            ;Max usable pixels (inclusive)

;;Most time-intensive part, if can reuse previous version, do so
if not keyword_set(each) then begin
    each = replicate(detection,maxdet) ;Values for individual pixels
    
    each.wrest = wave[imu:ilast]
    each.doppler = dopb
    vello = velocity-ndopb*dopb
    velhi = velocity+ndopb*dopb

    ;;Assign min values for every possible wavelength
    for ii=0L, maxdet-1 do begin
        dum = min(velocity-vello[imu],ilo,/absolute)
        dum = min(velocity-velhi[imu],ihi,/absolute)
        rng = lindgen(ihi-ilo+1)+ilo
        gd = where(error[rng] gt 0.,ngd) ;exclude gaps
        if ngd eq 0 then begin
            each[ii].wv_lim = [wave[ilo],wave[ihi]]
            each[ii].ew = 0.
            each[ii].sigew = 0.
            each[ii].flg = 0    ;just to be explicit
        endif else begin
            rng = rng[gd]

            ;;Record the EW and sigEW of every point
            sig = wave[imu]*dopb/(rt2*sol) ;A
            wght = 1./(rt2*rtpi*sig)*$
              exp(-0.5*((wave[rng]-wave[imu])/sig)^2)
            ;;Normalize weight 
            nwght = total(wght) ;*dwv[rng])/(2*sig*rtpi)

            ;;Modify 1-flux when flux < -|error| and flux > 1+|error|
            nfx = where(flux[rng] lt -abs(error[rng]),nnfx)
            pfx = where(flux[rng] gt 1.+abs(error[rng]),npfx)
            modfx = 1-flux[rng]
            if nnfx ne 0 then modfx[nfx] = 1+2*abs(error[rng[nfx]])
            if npfx ne 0 then modfx[pfx] = -2*abs(error[rng[pfx]])
            
            each[ii].wv_lim = [wave[ilo],wave[ihi]]
            each[ii].ew = total(modfx*dwv[rng]*wght)/nwght
            each[ii].sigew = sqrt(total((error[rng]*dwv[rng])^2*wght)/$
                                  nwght)
        endelse

        ;;Check S/N criterion with Gaussian-convolved EW/sigEW
        if each[ii].ew/each[ii].sigew ge 1. then $ ;will trim later
          each[ii].flg = 1      ;met detection threshold
        
        imu++
    endfor
endif                           ;each structure not supplied


;;Group positive detections into probable absorption features
srt = sort(each.wv_lim[0])
each = each[srt]
hit = where(each.flg eq 1,nhit) ;;used to be before the sort and unclear why and how that possibly worked 
auto = replicate(detection,nhit) ;max possible


if keyword_set(debug) then stop,'searchspec: Completed evaluating each pixel'


;;Churn through 'hits' and combine into features
clumps = replicate(-1L,nhit,2)

ii = 0L
while ii lt nhit-2 do begin
    done = 0
    clump = 0                   ;false
    rng = hit[ii]

    while not done do begin
        case group of         
            1: begin            ;default
                ;;If next hit rest wavelength <= upper wavelength limit of
                ;;current hit and visa versa, then consider still part
                ;;of one feature 
                if each[hit[ii+1]].wrest le each[hit[ii]].wv_lim[1] and $
                  each[hit[ii+1]].wv_lim[0] le each[hit[ii]].wrest then $
                  clump = 1 else clump = 0
            end
            2: begin
                ;;If next hit rest lower wavelength <= upper
                ;;wavelength limit of current hit, then consider still
                ;;part of one feature 
                if each[hit[ii+1]].wv_lim[0] le each[hit[ii]].wv_lim[1] then $
                  clump = 1 else clump = 0
            end
            3: begin            
                ;;If next hit is consecutive to current hit
                if hit[ii]+1 eq hit[ii+1] then clump = 1 $ ;true
                else clump = 0  ;false
            end
            else: begin
                stop,'searchspec: Not valid option group = ',group
            end
        endcase 

        ;;Clump 
        if clump  then rng = [rng,hit[ii+1]] $
        else begin
            clumps[ii,*] = [rng[0],rng[n_elements(rng)-1]]
            done = 1
        endelse 

        ii++
        if ii eq nhit-2 then begin
            done = 1
            ;;Process last clump
            if clump then clumps[ii-1,*] = [rng[0],rng[n_elements(rng)-1]]
        endif 
    endwhile                    ;not done
endwhile                        ;loop hits


;;;;;;;;;;;;;;;
;; Fit Features
;;;;;;;;;;;;;;;
gd = where(clumps[*,0] ne -1,nclumps)
clumps = clumps[gd,*]
for ii=0,nclumps-1 do begin
    ;;Troubleshooting
                                ;tst = where(each[rng].wrest ge 1081.0 and $
                                ;            each[rng].wrest le 1082.0,ntst)
    ntst = 0
    if ntst ne 0 then begin
        testing = 1
        print,min(each[clump[ii,0]:clump[ii,1]].wrest,max=mx1),mx1,$
          min(each[clump[ii,0]:clump[ii,1]].wv_lim,max=mx2),mx2
        stop,'searchspec testing: determined clump'     
    endif else testing = 0
    ;;End troubleshooting

    ;;Evaluate one feature (even single satisfactory points
    ;;will be subject to chi^2, nsigma and dopb (re)check

    dum = min(abs(wave-each[clumps[ii,0]].wv_lim[0]),ilo)
    dum = min(abs(wave-each[clumps[ii,1]].wv_lim[1]),ihi)
    auto[ii].wv_lim = [wave[ilo],wave[ihi]] ;set range #1
    

    ;;;;;;;;;;;;;;;
    ;;Fit Gaussian 
    ;;;;;;;;;;;;;;;
    if keyword_set(fitgauss) then begin
        ;;Initial estimates for Gaussian fit
        a1 = total((1.-(flux[ilo:ihi]<1.))*velocity[ilo:ihi])$
          /total(1.-(flux[ilo:ihi]<1.)) ;flux-weighted centroid
        dum = min(abs(velocity-a1),a0)
        if 1.-flux[a0] gt 0. then a0 = 1.-flux[a0] $
        else begin
            fxmx = max(flux[ilo:ihi],min=fxmn)
            a0 = fxmx-fxmn      ;Amplitude
        endelse 
        a2 = (velocity[ihi]-velocity[ilo])/(rt2*ndopb) ;width in sigma
        gfit = gaussfit(velocity[ilo:ihi],1.-flux[ilo:ihi],$
                        coeff,chisq=chisq,sigma=coefferr,$
                        measure_errors=error[ilo:ihi],yerror=gerr,$
                        nterms=3,estimates=[a0,a1,a2])

        ;;Prevent unphysical results or fit corruption due to
        ;;nearby feature (when gaussfit fails to converge) or
        ;;feature saturated and shouldn't be fit by Gaussian
        sat = where(flux[ilo:ihi] LT error[ilo:ihi]/5. OR $
                    flux[ilo:ihi] LT 0.05, nsat)

        neg = where(coeff le 0.,nneg) 
        if nneg eq 0 and coeff[1] le velocity[ihi] and $
          coeff[1] ge velocity[ilo] and $
          nsat lt saturate*(ihi-ilo) then begin
            linterp,velocity,wave,coeff[1],auto[ii].wrest
            auto[ii].doppler = coeff[2]*rt2
            auto[ii].doperr = coefferr[2]*rt2
            auto[ii].amp = coeff[0]
            auto[ii].amperr = coefferr[0]
            auto[ii].chisq = chisq
            ;;?should use fitted pars?
            if testing eq 1 then $
              stop,'Accept fit #1:',auto[ii].wrest,auto[ii].wv_lim
        endif else begin
            ;;Set reasonable values and bypass much of following
            auto[ii].wrest = total((1.-(flux[ilo:ihi]<1.))*wave[ilo:ihi])$
              /total(1.-(flux[ilo:ihi]<1.)) ;flux weighted centroid
            auto[ii].doppler = a2*rt2 ;best guess
            auto[ii].amp = a0
            auto[ii].chisq = -9.99 ;just bad
            if keyword_set(debug) then stop
            if not keyword_set(silent) then $
              print,'searchspec: Unreasonable fit, forcing wrest = ',$
              string(auto[ii].wrest,format='(f9.4)')
            if testing eq 1 then $
              stop,'Reject fit #1:',auto[ii].wrest,auto[ii].wv_lim
        endelse 



        ;;Chi^2 check #1
        if (auto[ii].chisq ge chisqtol[0] and $
            auto[ii].chisq le chisqtol[1]) then begin
            ;;Re-evaluate wavelength limits
            ;;Checked chi^2 tolerance in order to ensure that this
            ;;is just an improvement on an already good fit
            ;;(preserves potential for blends or assymmetric features)
            dum = min(abs(velocity-coeff[1]+ndopb*coeff[2]*rt2),ilo)
            dum = min(abs(velocity-coeff[1]-ndopb*coeff[2]*rt2),ihi)
            
            if ihi-ilo le 3 then begin
                ;;Revert to previous fit; too few independent params
                dum = min(abs(wave-each[clumps[ii,0]].wv_lim[0]),ilo)
                dum = min(abs(wave-each[clumps[ii,1]].wv_lim[1]),ihi)
                
                gfit = gaussfit(velocity[ilo:ihi],1.-flux[ilo:ihi],$
                                coeff,chisq=chisq,sigma=coefferr,$
                                measure_errors=error[ilo:ihi],$
                                yerror=gerr,$
                                nterms=3,estimates=[a0,a1,a2])
                if keyword_set(debug) then $
                  stop,'searchspec debug: ihi eq ilo'
                if testing eq 1 then stop,'Revert to fit #1 (too few):',$
                  auto[ii].wrest,auto[ii].wv_lim
            endif else begin
                ;;Can fit Gaussian
                gfit = gaussfit(velocity[ilo:ihi],1.-flux[ilo:ihi],$
                                coeff,chisq=chisq,sigma=coefferr,$
                                measure_errors=error[ilo:ihi],$
                                yerror=gerr,nterms=3,$
                                estimates=coeff)
                
                                ;?Put in a check here that if new
                                ;chisq more out of chisq tol, then abort?
                
                ;;Prevent unphysical results or fit corruption due to
                ;;nearby feature (when gaussfit fails to converge)
                ;;and contain new wrest to old wv_lim
                neg = where(coeff le 0.,nneg) 
                linterp,velocity,wave,coeff[1],dum
                if nneg eq 0 and coeff[1] le velocity[ihi] and $
                  coeff[1] ge velocity[ilo] and $
                  dum le auto[ii].wv_lim[1] and $
                  dum ge auto[ii].wv_lim[0] then begin
                    ;;Proceed
                    if not keyword_set(nochng) then begin
                        ;;Allow changing wavelength limits
                        ;;(enables overlap) (set range #2)
                        auto[ii].wv_lim = [wave[ilo],wave[ihi]] 
                    endif else begin
                        ;;Prevent changing of range but enable new
                        ;;fit pars 
                        dum2 = min(abs(wave-each[clumps[ii,0]].wv_lim[0]),ilo)
                        dum2 = min(abs(wave-each[clumps[ii,1]].wv_lim[1]),$
                                   ihi)
                    endelse 
                    auto[ii].wrest = dum
                    auto[ii].doppler = coeff[2]*rt2
                    auto[ii].doperr = coefferr[2]*rt2
                    auto[ii].amp = coeff[0]
                    auto[ii].amperr = coefferr[0]
                    auto[ii].chisq = chisq     
                    if testing eq 1 then stop,'Accept fit #2:',$
                      auto[ii].wrest,auto[ii].wv_lim
                endif else begin
                    ;;Revert to previous fit
                    dum = min(abs(wave-each[clumps[ii,0]].wv_lim[0]),ilo)
                    dum = min(abs(wave-each[clumps[ii,1]].wv_lim[1]),ihi)
                    
                    gfit = gaussfit(velocity[ilo:ihi],1.-flux[ilo:ihi],$
                                    coeff,chisq=chisq,sigma=coefferr,$
                                    measure_errors=error[ilo:ihi],$
                                    yerror=gerr,$
                                    nterms=3,estimates=[a0,a1,a2])
                    if keyword_set(debug) then $
                      stop,'searchspec debug: revert to previous fit'
                    if testing eq 1 then stop,$
                      'Revert to fit #1 (bad fit):',$
                      auto[ii].wrest,auto[ii].wv_lim
                endelse 
            endelse             ;ilo ne ihi-3
            if testing eq 1 then stop,'End chi^2 check #1:',$
              auto[ii].wrest,auto[ii].wv_lim
        endif                   ;chi^2 check #1

        ;;Check chi^2 #2 (potentially new boundaries)
        if (auto[ii].chisq ge chisqtol[0] and $
            auto[ii].chisq le chisqtol[1]) then $
          auto[ii].flg = 4 else $ ;Gaussian-weight
          auto[ii].flg = 2      ;boxcar-weight
    endif else begin            ;not keyword_set(fitgauss)
        ;;Assume values for auto structure w/o fitting
        a1 = total((1.-(flux[ilo:ihi]<1.))*velocity[ilo:ihi])$
          /total(1.-(flux[ilo:ihi]<1.)) ;flux-weighted centroid
        dum = min(abs(velocity-a1),a0)
        if 1.-flux[a0] gt 0. then a0 = 1.-flux[a0] $
        else begin
            fxmx = max(flux[ilo:ihi],min=fxmn)
            a0 = fxmx-fxmn      ;Amplitude
        endelse 
        a2 = (velocity[ihi]-velocity[ilo])/(rt2*2.*ndopb) ;width in sigma
        
        auto[ii].wrest = total((1.-(flux[ilo:ihi]<1.))*wave[ilo:ihi])$
          /total(1.-(flux[ilo:ihi]<1.)) ;flux weighted centroid
        auto[ii].doppler = a2*rt2 ;best guess
        auto[ii].amp = a0
        auto[ii].chisq = -9.99  ;just bad
        if keyword_set(debug) then stop
        if not keyword_set(silent) then $
          print,'searchspec: Unreasonable fit, forcing wrest = ',$
          string(auto[ii].wrest,format='(f9.4)')
        if testing eq 1 then $
          stop,'Reject fit #1:',auto[ii].wrest,auto[ii].wv_lim
        
        auto[ii].flg = 2 
    endelse 

    ;;Finally measure EW and sigEW
    case auto[ii].flg of 
        2: begin                ;boxcar extraction
            rng = lindgen(ihi-ilo+1)+ilo
            gd = where(error[rng] gt 0.,ngd) ;exclude gaps
            if ngd eq 0 then begin
                stop,'searchspec: should not have bad data here'
                auto[ii].ew = 0.
                auto[ii].sigew = 0.
                auto[ii].flg = 0
            endif else begin
                rng = rng[gd]
                auto[ii].ew = total((1-flux[rng])*dwv[rng])
                auto[ii].sigew = sqrt(total((error[rng]*dwv[rng])^2))
            endelse
        end
        4: begin                ;use Gaussian fit
            sig = auto[ii].wrest*auto[ii].doppler/(sol*rt2)
            sigerr = auto[ii].wrest/(sol*rt2)*auto[ii].doperr
            linterp,wave,velocity,auto[ii].wrest,vrest

            tlo = rt2*(velocity[ilo]-vrest)/auto[ii].doppler
            thi = rt2*(velocity[ihi]-vrest)/auto[ii].doppler
            gdiff = gaussint(thi) - gaussint(tlo)
            
            auto[ii].ew = rt2*rtpi*auto[ii].amp*sig*gdiff ;mA
            auto[ii].sigew = rt2*rtpi*gdiff*$
              sqrt((auto[ii].amperr[0]*sig)^2 + $
                   (auto[ii].amp*sigerr)^2) ;mA
        end
        else: stop,'searchspec: unexpected result, stopping'
    endcase 
    
    done = 1                    ;free to evaluate next feature
    if keyword_set(debug) then stop
    if testing eq 1 then stop,'Moving on:',$
      auto[ii].wrest,auto[ii].wv_lim
endfor                          ;clumps


;;Tie up loose ends
auto.ion = string(strtrim(auto.wrest,2),format='(a9)')
each.ion = string(strtrim(each.wrest,2),format='(a9)')
auto.ew = auto.ew*10^3          ;mA
auto.sigew = auto.sigew*10^3    ;mA
each.ew = each.ew*10^3          ;mA
each.sigew = each.sigew*10^3    ;mA

;;Doppler b cut can't be trusted b/c not everything really fit
gd = where(auto.flg gt 0 and auto.ew/auto.sigew ge nsigma,ndet,$
           complement=bd,ncomplement=nbd) ;trim
print,'searchspec: Eliminated ',nbd,' with S/N < ',nsigma

if ndet eq 0 then begin 
    auto = -1            ;empty
    nrej = 0
    print,'searchspec: No lines detected with S/N > ',nsigma
    if keyword_set(strct_fil) then begin 
       bdstrct_fil = strmid(strct_fil,0,strpos(strct_fil,'.',/reverse_search))+'_nofeatures.fits'
       goto,no_gd_data          ;just print structure
    endif 
    ostrct = -1
endif else begin
    auto = auto[gd]
    rslt = civ_search_overlap(auto,elim=elim,/divid) ;pair down duplicates
    nrej = n_elements(elim)
    if elim[0] eq -1 then nrej = 0
    if nrej ne 0 then reject = auto[elim]
    print,'searchspec: Eliminated ',n_elements(elim),' duplicates (saved)'
    auto = auto[rslt]
    ndet = n_elements(auto)
    print,'searchspec: Detected total ',ndet,' lines'
endelse

if keyword_set(comp) then begin        
    ;;Comparing different EW measurements
    search_compew,wave,flux,error,auto,comp,velocity=velocity

    ;;Print comparison
    if not keyword_set(silent) then begin
        print,''
        print,'Compare EW Measures:'
        printcol,comp.wrest,comp.ewbox,comp.sigewbox,comp.ewgwght,$
          comp.sigewgwght,comp.ewg,comp.sigewg,comp.flg,$
          format='(f9.4,3(2x,f7.1,1x,f7.1),2x,i1)'
    endif 
endif                           ;/comp

;;Write information to data file
if keyword_set(table) then begin
    tabhd = ''
    if keyword_set(name) then tabhd = [tabhd,name]
    tabhd = [tabhd,'No. sigma = '+string(nsigma)]
    tabhd = [tabhd,'Doppler b = '+string(dopb)]
    tabhd = [tabhd,'Chi^2 range: '+string(chisqtol)]
    tabhd = [tabhd,'No. Doppler b to determine width = '+string(ndopb)]
    if keyword_set(shift) then $
      tabhd = [tabhd,'Wavelength shifted: '+string(w0)+' '+string(dw0)] 

    close,/all
    openw,13,table
    search_wrtab,auto,filnum=13,head=tabhd

    if nrej ne 0 then begin
        tabhd = ['','~~~~~~~~','Rejected','~~~~~~~~']
        search_wrtab,reject,filnum=13,head=tabhd
    endif 
    close,13

    if not keyword_set(silent) then print,'searchspec: created ',table
endif 


if keyword_set(linlst) then begin
    if not keyword_set(instr) then instr = 0
    caldat,systime(/julian),mo,da,yr
    date=strtrim(da,2)+'.'+strtrim(mo,2)+'.'+strtrim(yr,2)

    head=strarr(2)
    if keyword_set(name) then $
      head[0] = string(name,nsigma,dopb,format='(a0,2x,f5.1,2x,f5.1)') $
    else head[0] = string(nsigma,dopb,format='(f5.1,2x,f5.1)') 
    head[1] = date

    if keyword_set(silent) then makelinlst,linlst,auto,head=head,/silent $
    else makelinlst,linlst,auto,head=head
endif


no_gd_data:                     ;no data with error > 0. or no auto features

;;Write multi-extension FITS file
header = strarr(1)
if keyword_set(bdstrct_fil) then $
   sxaddpar,header,'FILENAME',bdstrct_fil,'name of file' $
else if keyword_set(strct_fil) then $
   sxaddpar,header,'FILENAME',strct_fil,'name of file'
if keyword_set(name) then sxaddpar,header,'TARGNAME',name,'name of spectrum' $
else sxaddpar,header,'TARGNAME',' ','name of spectrum'
sxaddpar,header,'NSIGMA',nsigma,'Detection threshold parameter'
sxaddpar,header,'DOPB',dopb,'Doppler paramter (km/s)'
sxaddpar,header,'XSQLO',chisqtol[0],'Min reduced Chi^2 for good fit'
sxaddpar,header,'XSQHI',chisqtol[1],'Max reduced Chi^2 for good fit'
sxaddpar,header,'NDOPB',ndopb,'No. Doppler b to determine width'
sxaddpar,header,'GROUPING',group,'1: cent, 2: range, 3: consec'
if keyword_set(shift) then begin
   sxaddpar,header,'W0',w0,'Reference wavelength for shift'
   sxaddpar,header,'DWO',dw0,'Wavelength shift'
endif else begin
   sxaddpar,header,'W0',0.,'Reference wavelength for shift'
   sxaddpar,header,'DWO',0.,'Wavelength shift'
endelse 
   
if keyword_set(strct_fil) then begin
   if keyword_set(bdstrct_fil) then begin
      mwrfits,auto,bdstrct_fil,header,/create,/silent
      mwrfits,each,bdstrct_fil,/silent
      if nrej ne 0 then mwrfits,reject,bdstrct_fil,/silent
      if keyword_set(comp) then mwrfits,comp,bdstrct_fil,/silent
      if not keyword_set(silent) then print,'searchspec: created ',bdstrct_fil
   endif else begin
      mwrfits,auto,strct_fil,header,/create,/silent
      mwrfits,each,strct_fil,/silent
      if nrej ne 0 then mwrfits,reject,strct_fil,/silent
      if keyword_set(comp) then mwrfits,comp,strct_fil,/silent
      if not keyword_set(silent) then print,'searchspec: created ',strct_fil
   endelse
endif 
ostrct = auto         ; return structure


;;Display spectrum with extraction windows in red and wavelength
;;limits for each feature as blue points
if keyword_set(view) then begin
    print,'searchspec: plotting automatically detect lines'
    search_plot,wave,flux,auto
    if nrej ne 0 then begin
        print,'searchspec: plotting rejected duplicates'
        search_plot,wave,flux,reject
    endif
endif

if keyword_set(debug) then stop

;; Undo trimming
error = sverror
flux = svflux

end 
