;+ 
; NAME:
; wfccd_combspec
;    Version 1.0
;
; PURPOSE:
;   Combines multiple exposures of the same slit
;
; CALLING SEQUENCE:
;   
;   wfccd_combspec, wfccd, mask_id, exp_id
;
; INPUTS:
;   wfstrct     - WFCCD structure
;
; RETURNS:
;
; OUTPUTS:
;   wfspec      -  WFCCD fspec structure (fits file)
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   wfccd_combspec, wfccd, mask_id, exp_id
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   23-Apr-2002 Written by JXP
;-
;------------------------------------------------------------------------------

pro wfccd_combspec_out, wffspec

  ;; Output to ASCII
  printf, 56, wffspec.field+': '+$
    strtrim(wffspec.slit_id,2)+strtrim(wffspec.obj_id,2)
  for j=0L,wffspec.nexp-1 do begin
      printf, 56, FORMAT='(10x,f7.1,1x,2f10.3,1x,a25)',$
        wffspec.texp[j], $
        wffspec.wvmnx[j,0], $
        wffspec.wvmnx[j,1], $
        wffspec.obj_fil[j]
  endfor
  return
end

;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------
;------------------------------------------------------------------------------

pro wfccd_combspec, wfccd, mask_id, exp_id, SILENT=silent, PATCH=patch, $
                    OBJ_NM=OBJ_NM, LIST=list

;
  if  N_params() LT 2 AND not keyword_set( LIST ) then begin 
    print,'Syntax - ' + $
      'wfccd_combspec, wfccd, mask_id, [exp_id], /SILENT, OBJ_NM=, /PATCH [v1.0]'
    return
  endif 
  if keyword_set( LIST ) then $
    print, 'wfccd_combspec: Combining Obj files from this list -- ', list

;  Optional Keywords

; Set exp
  if not keyword_set( LIST ) then begin
      allexp = where(wfccd.type EQ 'OBJ' AND wfccd.flg_anly NE 0 AND $
                     wfccd.mask_id EQ mask_id)
;     want to allow exp_id==0!
      if n_elements(exp_id) gt 0 then exp = allexp[exp_id] else exp=allexp
      nexp = n_elements(exp)
  endif else begin
      mask_id = 99L
      readcol, list, files, FORMAT='A'
      nexp = n_elements(files)
  endelse
      


;;;;;;;;;;
; Open Files

  if not keyword_set( SILENT ) then $
    print, 'wfccd_combspec: Loading up the files...'

; Open Obj files
  for q=0L,nexp-1 do begin
      ; Read
      if not keyword_set( LIST ) then $
        tmp = xmrdfits(wfccd[exp[q]].obj_fil, 1, STRUCTYP='specobjstrct', /silent) $
      else $ ;; LIST
        tmp = xmrdfits(files[q], 1, STRUCTYP='specobjstrct', /silent)

      ; OBJ_NM
      if keyword_set(OBJ_NM) then begin
          gdobj = x_getobjnm(tmp, obj_nm)
          ;; Checks
          if gdobj EQ -1 then begin
              print, 'wfccd_combspec: No obj ', obj_nm, ' in exp ', exp[q]
              continue
          endif
          if tmp[gdobj].flg_flux EQ 0 then begin
              print, 'wfccd_combspec: Data not good or not fluxed!'
              return
          endif 
      endif else begin
          ; Keep only good obj (require fluxed)
          gdobj = where(tmp.flg_flux NE 0 AND tmp.flg_anly NE 0, ngd)
          if ngd EQ 0 then begin
              print, 'wfccd_combspec: Data not good or not fluxed!'
              return
          endif
      endelse
      ; Add to total structure
      if q EQ 0 then wfobj = tmp[gdobj] else wfobj = [wfobj, tmp[gdobj]]
  endfor

  nspec = n_elements(wfobj)


;; ASCII OUTPUT ;;;

  if mask_id LT 10L then txtfil = 'Extract/Fspec_0'+$
    string(mask_id,format='(i1)')+'.txt' $
  else txtfil = 'Extract/Fspec_'+string(mask_id,format='(i2)')+'_comb.txt'
  close, /all
  if not keyword_set( OBJ_NM ) then openw, 56, txtfil

;;; CREATE FINAL ;;;

  if mask_id LT 10L then outfil = 'Extract/Fspec_0'+$
    string(mask_id,format='(i1)')+'.fits' $
  else outfil = 'Extract/Fspec_'+string(mask_id,format='(i2)')+'.fits'

  if not keyword_set( OBJ_NM ) then begin
      tmp = { wfccdfspecstrct }
      wffspec = replicate(tmp, nspec)
  endif else begin
      wfccd_wrfspec, wffspec, outfil, /read
  endelse

;;; FIND UNIQUE SLITS ;;;;

  if not keyword_set( OBJ_NM ) then begin
      uni_slit = uniq(wfobj.slit_id, sort(wfobj.slit_id))
      allslit = wfobj[uni_slit].slit_id
      nslit = n_elements(allslit)

  ; Loop on Slits

      cnt = 0L
      for q=0L, nslit-1 do begin
          slit = allslit[q]
          
          ;; Science first
          sciobj = where(wfobj.slit_id EQ slit AND $
                         wfobj.obj_id EQ 'a' AND wfobj.flg_anly NE 0, nsci)
          if nsci NE 0 then begin
;; MRB-- only use ones where there is actual good data to use
              isgd=bytarr(n_elements(wfobj))+1
              for isci=0L,nsci-1 do begin
                  npix=wfobj[sciobj[isci]].npix
                  gdwv= $
                    where(wfobj[sciobj[isci]].wave[0:npix-1] GT 4500. AND $
                          wfobj[sciobj[isci]].wave[0:npix-1] LT 7000. AND $
                          wfobj[sciobj[isci]].sig[0:npix-1] GT 0., ngd)
                  if(ngd eq 0) then isgd[sciobj[isci]]=0
              endfor
              sciobj = where(wfobj.slit_id EQ slit AND $
                             wfobj.obj_id EQ 'a' AND $
                             wfobj.flg_anly NE 0 AND $
                             isgd gt 0, nsci)
              if nsci eq 0 then cnt=cnt+1L
          endif

          if nsci NE 0 then begin
;; Copy
              wffspec[cnt].nexp = nsci
              for i=0L,nsci-1 do wffspec[cnt].texp[i] = wfobj[sciobj[i]].exp
              tmpstr = wffspec[cnt]
              copy_struct, wfobj[sciobj[0]], tmpstr, EXCEPT=["wave","fx","var"]
              wffspec[cnt] = tmpstr
              for i=0L,nsci-1 do begin
                  ipos = strpos(wfobj[sciobj[i]].spec2d_fil, 'ccd')
                  obj_fil = 'Extract/Obj_'+strmid(wfobj[sciobj[i]].spec2d_fil, $
                                                  ipos+3)
                  wffspec[cnt].obj_fil[i] = obj_fil
; add 2d spectra names
                  sp2d_fil = 'Final/f_ccd'+strmid(wfobj[sciobj[i]].spec2d_fil, $
                                                  ipos+3)
                  wffspec[cnt].spec2d_fil[i] = sp2d_fil
             endfor
              ;; NPIX
              npix = wffspec[cnt].npix
              ;; Coadd
              if nsci EQ 1 then begin
                  wffspec[cnt].wave[0:npix-1] = wfobj[sciobj].wave[0:npix-1]
                  wffspec[cnt].fx[0:npix-1] = wfobj[sciobj].flux[0:npix-1]
                  wffspec[cnt].var[0:npix-1] = $
                    double(wfobj[sciobj].sig[0:npix-1])^2
              endif else begin
                  var = double(wfobj[sciobj].sig[0:npix-1])^2
                  for kk=0L,nsci-1 do begin
                      a = where(wfobj[sciobj[kk]].sig[0:npix-1] LT 0., na)
                      if na GT 0 then var[a,kk] = -1.
                  endfor
                  x_combspec, wfobj[sciobj].flux[0:npix-1], var, $
                    fflux, fvar, WAVE=wfobj[sciobj[0]].wave[0:npix-1], $
                    NRMFLUX=[4500., 7000.]
                  wffspec[cnt].wave[0:npix-1] = wfobj[sciobj[0]].wave[0:npix-1]
                  wffspec[cnt].fx[0:npix-1] = temporary(fflux[0:npix-1])
                  wffspec[cnt].var[0:npix-1] = temporary(fvar[0:npix-1])
              endelse
              ;; Output
              wfccd_combspec_out, wffspec[cnt]
              cnt = cnt+1
          endif
                  
          ;; SERENDIP
          sdpobj = where(wfobj.slit_id EQ slit AND $
                         wfobj.obj_id NE 'a' AND wfobj.flg_anly NE 0, nsdp)
          if nsdp EQ 0 then continue
;; MRB-- only use ones where there is actual good data to use
          isgd=bytarr(n_elements(wfobj))+1
          for isdp=0L,nsdp-1 do begin
              npix=wfobj[sdpobj[isdp]].npix
              gdwv= $
                where(wfobj[sdpobj[isdp]].wave[0:npix-1] GT 4500. AND $
                      wfobj[sdpobj[isdp]].wave[0:npix-1] LT 7000. AND $
                      wfobj[sdpobj[isdp]].sig[0:npix-1] GT 0., ngd)
              if(ngd eq 0) then isgd[sdpobj[isdp]]=0
          endfor
          sdpobj = where(wfobj.slit_id EQ slit AND $
                         wfobj.obj_id NE 'a' AND $
                         wfobj.flg_anly NE 0 AND $
                         isgd gt 0, nsdp)
          if nsdp EQ 0 then continue

          srt = sort(wfobj[sdpobj].ycen)

          jj = 0L
          kk= 0L
          while(jj LT nsdp) do begin

              ;; Find all within 2pix
              kp = where( abs(wfobj[sdpobj].ycen-wfobj[sdpobj[srt[jj]]].ycen) $
                          LT 2., ngd)
              jj = jj+ngd
              kk = kk+1
              ;; Reset
              gd = sdpobj[kp]
              ;; Copy
              wffspec[cnt].nexp = ngd
              for i=0L,ngd-1 do wffspec[cnt].texp[i] = wfobj[gd[i]].exp
              tmpstr = wffspec[cnt]
              copy_struct, wfobj[gd[0]], tmpstr, EXCEPT=["wave","fx","var"]
              wffspec[cnt] = tmpstr
              for i=0L,ngd-1 do begin
                  ipos = strpos(wfobj[gd[i]].spec2d_fil, 'ccd')
                  obj_fil = 'Extract/Obj_'+strmid(wfobj[gd[i]].spec2d_fil, ipos+3)
                  wffspec[cnt].obj_fil[i] = obj_fil
              endfor
              ;; OBJ name
              wffspec[cnt].obj_id = x_objnumid(kk)
              ;; NPIX
              npix = wffspec[cnt].npix
              ;; Coadd
              if ngd EQ 1 then begin
                  wffspec[cnt].wave[0:npix-1] = wfobj[gd].wave[0:npix-1]
                  wffspec[cnt].fx[0:npix-1] = wfobj[gd].flux[0:npix-1]
                  gdpix = where(wfobj[gd].sig[0:npix-1] GT 0., COMPLEMENT=badpix, $
                                NCOMPLEMENT=nbad)
                  wffspec[cnt].var[gdpix] = double(wfobj[gd].sig[gdpix])^2
                  if nbad NE 0 then $
                    wffspec[cnt].var[badpix] = double(wfobj[gd].sig[badpix])
              endif else begin
                  var = double(wfobj[gd].sig[0:npix-1])^2
                  for kk=0L,ngd-1 do begin
                      a = where(wfobj[gd[kk]].sig[0:npix-1] LT 0., na)
                      if na GT 0 then var[a,kk] = -1.
                  endfor
                  x_combspec, wfobj[gd].flux[0:npix-1], var, $
                    fflux, fvar, WAVE=wfobj[gd[0]].wave[0:npix-1], $
                    NRMFLUX=[4500., 7000.]
                  wffspec[cnt].wave[0:npix-1] = wfobj[gd[0]].wave[0:npix-1]
                  wffspec[cnt].fx[0:npix-1] = temporary(fflux[0:npix-1])
                  wffspec[cnt].var[0:npix-1] = temporary(fvar[0:npix-1])
              endelse
              cnt = cnt+1
              ;; Output
              wfccd_combspec_out, wffspec[cnt]
          endwhile
      endfor
  endif else begin         ;;;;;; OBJ_NM ;;;;;;
      slen = strlen(obj_nm)
      obj = strmid(obj_nm,slen-1)
      slit = long(strmid(obj_nm,0,slen-1))
      gd = where(wfobj.slit_id EQ slit AND $
                 wfobj.obj_id EQ obj AND wfobj.flg_anly NE 0, ngd)
      if ngd EQ 0 then begin
          print, 'wfccd_combspec: No spectra found!'
          return
      endif else begin
          ;; Set cnt
          cnt = where(wffspec.slit_id EQ slit AND wffspec.obj_id EQ obj )
          ;; Copy
          wffspec[cnt].nexp = ngd
          for i=0L,ngd-1 do wffspec[cnt].texp[i] = wfobj[gd[i]].exp
          tmpstr = wffspec[cnt]
          copy_struct, wfobj[gd[0]], tmpstr, EXCEPT=["wave","fx","var"]
          wffspec[cnt] = tmpstr
          for i=0L,ngd-1 do begin
              ipos = strpos(wfobj[gd[i]].spec2d_fil, 'ccd')
              obj_fil = 'Extract/Obj_'+strmid(wfobj[gd[i]].spec2d_fil, ipos+3)
              wffspec[cnt].obj_fil[i] = obj_fil
          endfor
          ;; NPIX
          npix = wffspec[cnt].npix
          ;; Coadd
          if ngd EQ 1 then begin
              wffspec[cnt].wave[0:npix-1] = wfobj[gd].wave[0:npix-1]
              wffspec[cnt].fx[0:npix-1] = wfobj[gd].flux[0:npix-1]
              gdpix = where(wfobj[gd].sig[0:npix-1] GT 0., COMPLEMENT=badpix, $
                            NCOMPLEMENT=nbad)
              wffspec[cnt].var[gdpix] = double(wfobj[gd].sig[gdpix])^2
              if nbad NE 0 then $
                wffspec[cnt].var[badpix] = double(wfobj[gd].sig[badpix])
          endif else begin
              var = double(wfobj[gd].sig[0:npix-1])^2
              for kk=0L,ngd-1 do begin
                  a = where(wfobj[gd[kk]].sig[0:npix-1] LT 0., na)
                  if na NE 0 then var[a,kk] = -1.
              endfor
              x_combspec, wfobj[gd].flux[0:npix-1], var, $
                fflux, fvar, WAVE=wfobj[gd[0]].wave[0:npix-1], $
                NRMFLUX=[4500., 7000.]
              wffspec[cnt].wave[0:npix-1] = wfobj[gd[0]].wave[0:npix-1]
              wffspec[cnt].fx[0:npix-1] = temporary(fflux)
              wffspec[cnt].var[0:npix-1] = temporary(fvar)
          endelse

          cnt = cnt+1
          ;; Output
;          wfccd_combspec_out, wffspec[cnt]
      endelse

  endelse


;;;; OUTPUT  ;;;;



  if not keyword_set( OBJ_NM ) then wfccd_wrfspec, wffspec[0:cnt-1], outfil $
  else wfccd_wrfspec, wffspec, outfil


  close, /all

  print, 'wfccd_combspec:  All done!'


  return
end
  

