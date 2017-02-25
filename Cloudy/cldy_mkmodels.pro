;+
; NAME:
; cldy_mkmodels
;  V1.1
;
; PURPOSE:
;    Run a suite of cloudy models
; CALLING SEQUENCE:
;
;  cldy_mkmodels, infil, fluxfil, z, OUTDIR=,INDIR=,CLOUDYDIR=,NHI=,METALS=,U=,HDEN=
;                                    /CLOBBER
;
;  INPUT:
;    infil:  
;       input file with:
;           <Variable name>  <desired minimum value> <desired maximum value> <step size>
;
;
;           variables which can be input  are:  NHI, nH, U, metals
;           NOTE:  MUST use these specific variable tags
;
;     fluxfil: cldy_cuba output file defining flux
;     
;     z:  redshift 
;
;     NHI, METALS, U, HDEN:  specify a single value for these variables to be used 
;                            in all models run (HDEN is for NH)
;
;     OUTDIR:  output file directory
;     INDIR:  input file directory
;     TITLE:  cloudy title (default: model)
;     FILENAME: input/output file name lead (default: cloudy)
;     /CLOBBER: clobber previous runs with same FILENAME
;     /NODOUB: Do not use double optical depths
;     /OLD_VER:  Assumes version < v6 [not recommended!]
;     FNU:  Do not vary U.  Use this variable to set log f(nu)=4*pi*Jnu at 912 [cgs]
;
;  NHI 16.0 18.0 0.25
;  metals -1.0 0.0 1.0
;  U -6.0 0.0 0.25
;  nH -2.0 0.0 2.0
;
;
; REVISION HISTORY:
;   10-Sep-2005 Written by GEP
;-
;------------------------------------------------------------------------------

pro cldy_mkmodels, infil, fluxfil, z, OUTDIR=outdir, $
                   INDIR=indir, CLOUDYDIR=cloudydir, CIE=cie, $
                   NHI=nhi, METALS=metals, U=u, HDEN=hden, $
                   TITLE=title, FILNAME=filename, SHAPECMD=shapecmd, $
                   CLOBBER=clobber, NODOUB=nodoub, OLD_VER=old_ver, $
                   NPROC=nproc, FNU=fnu

  if (N_params() LT 3) then begin
      print, 'Syntax - cldy_mkmodels, infil, fluxfil, z, /OLD_VER, NPROC='
      return
  endif

  if not keyword_set(fluxfil) then begin
      if not (keyword_set(CIE) or keyword_set(SHAPECMD)) then stop
  endif
  if not keyword_set(CLOUDYDIR) then cloudydir=getenv('CLDYDIR')
  if not keyword_set(OUTDIR) then outdir = ''
  if not keyword_set(INDIR) then indir = ''
  if not keyword_set(NPROC) then nproc = 1L

  if not keyword_set(infil) then begin
    print, 'infil not set!'
    stop
  endif

  if not keyword_set(FILENAME) then filename = 'cloudy'


  readcol, infil, varp, minp, maxp, stepp, FORMAT='A,D,D,D'

  nvars = n_elements(minp)

  if not keyword_set(nhi) then nhi = 0
  if not keyword_set(nh) then nh = 0
  if not keyword_set(metals) then metals = 0
  if not keyword_set(u) then u = 0
  if not keyword_set(heat) then heat = 0

  if keyword_set(HDEN) then nh = hden

  for i=0L, nvars-1 do begin
      g = long((maxp[i]-minp[i]+1d-7)/stepp[i])+1
      case varp[i] of
          'NH': nhtot = findgen(g)*stepp[i]+minp[i]  ;; Total NH column
          'NHI': nhi = findgen(g)*stepp[i]+minp[i] ;; NHI column
          'nH': nh = findgen(g)*stepp[i]+minp[i]
          'metals': metals = findgen(g)*stepp[i]+minp[i]
          'U': u = findgen(g)*stepp[i]+minp[i]
          'HEAT': heat = findgen(g)*stepp[i]+minp[i]
          else: begin
              print, 'cldy_mkmodels:  not a valid variable!: '+varp[i]
              stop
          end
      endcase
  endfor

  ;; Stopping (JXP)
  if not keyword_set(nhi) then begin
     if not keyword_set(NHTOT) then stop else begin
        flg_NH = 1
        nhi = nhtot
     endelse
  endif else flg_NH = 0
  if n_elements(nh) LT 1 then stop
  if n_elements(metals) LT 1 then stop
  nu = 0
  if not keyword_set(u) then begin
     if keyword_set(FNU) then begin
        u = [-99.]
        nu = 1
     endif 
  endif else nu = n_elements(u)
  if nu EQ 0 then stop

  nnhi = n_elements(nhi)  ;; Overloaded with NH (set by flg_NH)
  nmetals = n_elements(metals)
  nnh = n_elements(nh)
  nheat = n_elements(heat)

  ntot = nnhi * nmetals * nu * nnh * nheat
  count = 0
  tot_count=0

  for i=0,nnhi-1 do begin
      for j=0, nmetals-1 do begin
          for k=0, nu-1 do begin
              for l=0, nnh-1 do begin
                 for m=0, nheat-1 do begin
                    ;; Counters (for NPROC)
                    tot_count = tot_count + 1
                    
                    if not keyword_set(SILENT) then begin
                       print, 'ijkl = ', i,j,k,l, nnhi,nmetals,nu,nnh, $
                              FORMAT='(a7,1x,8i5)'
                       if nheat GT 1 then print, 'm = ', m, nheat
                    endif
                    infil=filename+'_'+strtrim(string(i),2)+'_'+ $
                          strtrim(string(j),2)+'_'+ $
                          strtrim(string(k),2)+'_'+strtrim(string(l),2)+'_'+$
                          strtrim(string(m),2)+'.in'
                    
                    ;; Check for output file
                    outfil=filename+'_'+strtrim(string(i),2)+'_' $
                        +strtrim(string(j),2)+'_' $
                           +strtrim(string(k),2)+'_'+strtrim(string(l),2)+'_'$
                           +strtrim(string(m),2)+'.out'
                    if (x_chkfil(outfil,/silent) OR $
                        x_chkfil(outdir+outfil,/silent)) AND $
                       not keyword_set(CLOBBER) $
                    then continue
                    
                    ;; Heating?
                    if nheat GT 1 then iheat = heat[m] + nh[l]
                    
                    ;; Input file
                    if not ((x_chkfil(infil,/silent) OR $
                             x_chkfil(indir+infil,/silent)) AND $
                            not keyword_set(CLOBBER)) then begin
                       params = [z,nhi[i],metals[j],u[k],nh[l]]
                       cldy_wrin, params, fluxfil, OUTFIL=infil, TITLE=title, $
                                  NODOUB=nodoub, SHAPECMD=shapecmd, CIE=cie, OLD_VER=old_ver, $
                                  FNU=fnu, flg_NH=flg_NH, HEAT=iheat
;                      if keyword_set(v0602a) then $
;                        cldy_wrin, params, fluxfil, OUTFIL=infil, TITLE=title,/v0602a, NODOUB=nodoub $
;                        else cldy_wrin, params, fluxfil, OUTFIL=infil, TITLE=title, NODOUB=nodoub
                    endif
                    
                    if keyword_set(INDIR) then $
                       spawn, 'mv '+infil+' '+indir
                    
                    count=count+1
                    if count LT NPROC AND TOT_COUNT NE NTOT then sbkg = '&' else begin
                       sbkg = ''
                       count = 0
                    endelse
                    
                    ;; Run Cloudy
                    spawn, CLOUDYDIR+'/cloudy.exe < '+indir+infil+'> '+outdir+outfil+sbkg
                 endfor
              endfor
          endfor
      endfor
  endfor

end
