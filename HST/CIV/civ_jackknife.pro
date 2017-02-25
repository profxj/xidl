;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; civ_jackknife.pro               
; Author: Kathy Cooksey                      Date: 4 Mar 2009
; Project: CIV absorption at z < 1 with Xavier Prochaska
; Description: Jackknife resampling to measure sample bias,
;              variance, and covariance matrix
; Input: 
; Optional:
; Output: 
; Example:
; History:
;   4 Mar 2009 -- by KLC; borrowed from
;                 SDDS/DLA/DR3/Analysis/pro/jack_corrm.pro
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
@civ_calcewn                    ; resolve civ_calcewn_ndblt()
pro civ_jackknife,option,matrix_only=matrix_only,clobber=clobber,$
                  view=view,by_los=by_los
  plot = 1
  if not keyword_set(option) then option = 0

  if keyword_set(by_los) then $
     cd,getenv('MLSS_DIR')+'/analysis/inputs/jackknife_los/' $
  else cd,getenv('MLSS_DIR')+'/analysis/inputs/jackknife/'
  
  case option of 
     0: begin
        ;; Full sample
        zlim=[0.,10]
        zstr = ''
        sens_fil =  getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal.fits'
     end
     1: begin
        ;; Lower-z bin
        zlim=[0.,0.6]
        zstr = 'zlt0.6'
        sens_fil =  getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal_zlt0.6.fits'
     end
     2: begin
        ;; Higher-z bin
        zlim=[0.6,10.]
        zstr = 'zgt0.6'
        sens_fil =  getenv('MLSS_DIR')+'/analysis/complete/cmplt_cumx_nogal_zgt0.6.fits'
     end
  endcase 

  ;; Loop rating = 6 and rating = 5,6
  for ii=0,1 do begin
     if ii eq 0 then begin
        civ_group,civcand,rating=[6],flg_sys=384,zlim=zlim,nciv=ntot
        root = 'rtgeq6'
     endif else begin
        civ_group,civcand,rating=[5,6],flg_sys=384,zlim=zlim,nciv=ntot
        root = 'rtg56'
     endelse 

     if keyword_set(by_los) then begin
        ;; Exclude lines of sight and then test
        civcand.qso = strtrim(civcand.qso,2)
        los = uniq(civcand.qso,sort(civcand.qso))
        ntot = n_elements(los)
     endif 
     
     ;; Output names
     fn_fil = strarr(ntot)
     few_fil = strarr(ntot)
     fnpars = dblarr(2,ntot)
     fewpars = dblarr(2,ntot)

     ;; Loops
     for jj=0,ntot-1 do begin
        num = strtrim(jj,2)

        if keyword_set(by_los) then begin
           rng = where(civcand.qso ne civcand[los[jj]].qso)
        endif else begin
           rng = lindgen(ntot-1)
           if jj lt ntot-1 then $
              rng[jj:ntot-2] = rng[jj:ntot-2] + 1
        endelse 

        ;; Ranges 
        n = civ_calcewn_ndblt(civcand[rng],signcolm=sn,flg_colm=fn,/silent)
        nintlim = [min(n,imn)-sn[imn], 10.^14.3, 10.^16]
        ewintlim = [min(civcand[rng].ew[0],imn,max=mx,subscript_max=imx)-$
                    civcand[rng[imn]].sigew[0],mx+civcand[rng[imx]].sigew[0],1.e4]
        ewintlim[2] = ewintlim[1]

        fn_fil[jj] = 'fn_'+root+'_'+zstr+'fit'+num+'.fits'
        few_fil[jj] = 'few_'+root+'_'+zstr+'fit'+num+'.fits'

        if not keyword_set(matrix_only) then begin
           ;; f(N) 
           if keyword_set(clobber) then $
              civ_maxsngpow,civcand[rng],sens_fil,/ncolm,ostrct_fil=fn_fil[jj],$
                            plot=plot,intlim=nintlim $
           else begin
              test = file_search(fn_fil[jj],count=ntest)
              if ntest eq 0 then $
                 civ_maxsngpow,civcand[rng],sens_fil,/ncolm,ostrct_fil=fn_fil[jj],$
                               plot=plot,intlim=nintlim
           endelse 
           
           ;; f(EW) 
           if keyword_set(clobber) then $
              civ_maxsngpow,civcand[rng],sens_fil,/ew,ostrct_fil=few_fil[jj],$
                            plot=plot,intlim=ewintlim $
           else begin
              test = file_search(few_fil[jj],count=ntest)
              if ntest eq 0 then $
                 civ_maxsngpow,civcand[rng],sens_fil,/ew,ostrct_fil=few_fil[jj],$
                               plot=plot,intlim=ewintlim
           endelse 

        endif                   ; matrix_only = 0

        ;; Store results
        fitstrct = xmrdfits(fn_fil[jj],1,/silent)
        fnpars[*,jj] = [fitstrct.coeff,fitstrct.alpha]

        fitstrct = xmrdfits(few_fil[jj],1,/silent)
        fewpars[*,jj] = [fitstrct.coeff,fitstrct.alpha]
     endfor                     ; loop while removing one

     ;; Correlation matrix
     wgt = double(ntot-1)/double(ntot)    ; b/c regular jackknife
     fnmeanpars = [mean(fnpars[0,*]),mean(fnpars[1,*])]
     fnsigpars = sqrt( wgt * $
                       [total((fnpars[0,*]-fnmeanpars[0])^2),$
                        total((fnpars[1,*]-fnmeanpars[1])^2)] )
     fncovar = dblarr(2,2)

     fewmeanpars = [mean(fewpars[0,*]),mean(fewpars[1,*])]
     fewsigpars = sqrt( wgt * $
                        [total((fewpars[0,*]-fewmeanpars[0])^2),$
                         total((fewpars[1,*]-fewmeanpars[1])^2)] )
     fewcovar = dblarr(2,2)

     for jj=0,1 do begin
        for kk=0,1 do begin
           fncovar[jj,kk] = wgt * $
                            total( (fnpars[jj,*] - fnmeanpars[jj]) * $
                                   (fnpars[kk,*] - fnmeanpars[kk]) ) / $
                            (fnsigpars[jj] * fnsigpars[kk])
           fewcovar[jj,kk] = wgt * $
                             total( (fewpars[jj,*] - fewmeanpars[jj]) * $
                                    (fewpars[kk,*] - fewmeanpars[kk]) ) / $
                             (fewsigpars[jj] * fewsigpars[kk])
        endfor                  ; loop kk
     endfor                     ; loop jj

     if keyword_set(view) then begin
        if view eq 1 then nbins = 10 else nbins = view
        for ii=0,1 do begin
           print,'civ_jackknife: column density ',ii
           hist = histogram(fnpars[ii,*],loc=loc,nbins=nbins)
           x_splot,[loc[0],loc,loc[nbins-1]],[0,hist,0],psym1=10,/block

           print,'civ_jackknife: equivalent width ',ii
           hist = histogram(fewpars[ii,*],loc=loc,nbins=nbins)
           x_splot,[loc[0],loc,loc[nbins-1]],[0,hist,0],psym1=10,/block
        endfor 
     endif                      ; /view

     ;; Output
     ofil = root+'_'+zstr+'corrm.fits'

     strct = {pars:fnpars, meanpars:fnmeanpars, $
              sigpars:fnsigpars, covar:fncovar}
     mwrfits,strct,'fn_'+ofil,/create,/silent
     print,'civ_jackknife: created fn_'+ofil

     strct = {pars:fewpars, meanpars:fewmeanpars, $
              sigpars:fewsigpars, covar:fewcovar}
     mwrfits,strct,'few_'+ofil,/create,/silent
     print,'civ_jackknife: created few_'+ofil

  endfor                        ; loop rating=6 and rating=5,6

end
