;+ 
; NAME:
; fuse_calccolm
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;   
;   lowzovi_prsdat, stucture, filename
;
; INPUTS:
;
; RETURNS:
;   structure      - IDL structure
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  LIST - File
;  ION - Input ionic column densities
;  NOELM - Supress inputting Elemental values
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calccolm, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_calccolm, strct_fil, instr_list, abs_list, NSIG=nsig

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'fuse_calccolm, strct, instr_list, abs_list (v1.0)' 
    return
  endif 

  if not keyword_set( NSIG ) then nsig = 3.

;  Read instrument file list
  readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'

  nlist = n_elements(instr_fil)
  if nlist LT 7 then stop

;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)
  msk = lonarr(n_elements(strct)) 

; Read abs files
  readcol, abs_list, abs_fil, format='a'
  nablist = n_elements(abs_fil)

  for nn=0L,nablist-1 do begin

      msk[*] = -1L
; Open abslist
      close, /all
      openr, 1, abs_fil[nn]

      ;; Redshift
      readf, 1, zabs, format='(f12.7)'
      
      ;; Set by hand
      readf, 1, nhand, format='(i3)'
      if nhand NE 0 then stop
      
      ;; Calculate
      readf, 1, ncalc, format='(i3)'
      
      wv_lin = dblarr(ncalc)
      flg_lin = lonarr(ncalc)
      dumd = 0.d
      dumi = 0L
      ;; Read list
      for i=0L,ncalc-1 do begin
          readf, 1, dumd, dumi
          wv_lin[i] = dumd
          ;; Mask
          a = where(abs(strct.zabs-zabs) LT 0.0002 AND $
                    abs(strct.wrest-wv_lin[i]) LT 0.003, na)
          if na NE 1 then begin
              print, 'na', na
              print, zabs, wv_lin[i]
              stop
          endif
          msk[a] = a[0]
          flg_lin[i] = dumi
      endfor
      
      ;; Create colm arrays
      all_N = dblarr(ncalc, 10)
      all_sN = dblarr(ncalc, 10)
      
      ;; Subset of structure
      gd = where(msk GE 0L)
      
;  Loop on instruments
      for qq=0L,nlist-1 do begin
          ;; Find lines
          lin = where(strct[gd].instr MOD 2^(qq+1) GT (2^qq-1), nlin)
          if nlin NE 0 then begin
              print, 'fuse_calcen: Reading ', instr_fil[qq]
              ;; Open data
              if qq LE 6 then $
                fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, NPIX=npix, inflg=3)$
              else begin        ; STIS
                  spos = strpos(instr_fil[qq], 'f.fits')
                  sig_fil = strmid(instr_fil[qq], 0, spos)+'e.fits'
                  fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, NPIX=npix, $
                                  fil_sig=sig_fil, inflg=0)
              endelse
              ;; Shift wave
              wave = wave + inst_dw[qq]*inst_w0[qq]/wave
              ;; Sort?
              srt = sort(wave)
              wave = wave[srt]
              fx = fx[srt]
              sig = sig[srt]
              ;; dwv
              dwv = wave - shift(wave,1)
          endif
          
          for ii=0L,nlin-1 do begin
              ;; Index
              jj = gd[lin[ii]]
              
              ;; Find pixmnx
              mn = min(abs(strct[jj].wv_lim[0] - wave), pmin)
              mx = min(abs(strct[jj].wv_lim[1] - wave), pmax)
              
              ;; Calculate N
              x_aodm, wave[pmin:pmax], fx[pmin:pmax], sig[pmin:pmax], $
                strct[jj].wrest, clm, sig_clm
              
              all_N[lin[ii],qq] = clm
              all_sN[lin[ii],qq] = sig_clm
          endfor
      endfor
      
      ;; Combine N
      for ii=0L,ncalc-1 do begin
          a = where(all_sN[ii,*] GT 0., na)
          jj = gd[ii]
          case na of
              0: stop
              1: begin          ; One value
                  weight_N = all_N[ii,a[0]]
                  weight_sN = all_sN[ii,a[0]]
              end
              else: begin       ; Weighted mean
                  twgt = total(1./all_sN[ii,a]^2)
                  weight_N = $
                    total(all_N[ii,a]/all_sN[ii,a]^2) / twgt
                  weight_sN = 1./sqrt(twgt)
                  ;; Fill up structure
                  ;; Check chi2 for divergent lines
;              chi2 = total( (strct[jj].EW[0]-strct[jj].EW[a+1])^2 / $
;                            strct[jj].sigEW[a+1]^2 ) 
;              if chisqr_pdf(chi2, na-1) GT 0.9 then $
;                    print, 'Trouble with: ', strct[jj].wrest, strct[jj].zabs
              end
          endcase  

          ;; Nsig
          if weight_N LT nsig*weight_sN AND strct[jj].flg NE 3 then begin
              weight_N = nsig*weight_sN
              ;; Adjust the flag
              if strct[jj].flg MOD 8 LT 3 then strct[jj].flg = strct[jj].flg + 4
          endif

          ;; Log
          if weight_N GT 0. then colm = alog10( weight_N ) else colm = weight_N
          lgvar = ((1.d / (alog(10.0)*weight_N))^2)*(weight_sN^2)
          sig_colm = sqrt(lgvar)
          ;; Fill up
          strct[jj].Ncolm = colm
          strct[jj].sigNcolm = sig_colm

          ;; Check flags
          if strct[jj].flg NE flg_lin[ii] AND flg_lin[ii] NE 1  AND $
            strct[jj].flg MOD 2 EQ 1 then begin
              print, 'fuse_calccolm: The flags dont match!', strct[jj].wrest
              print, flg_lin[ii], strct[jj].flg
              stop 
          endif
      endfor
      printcol, strct[gd].ion, strct[gd].Ncolm, strct[gd].sigNcolm, strct[gd].flg
  endfor

  close, /all
      
      ;; Overwrite
      if not keyword_set( OUTFIL ) and NOT keyword_set( NOSV ) then $
        mwrfits, strct, strct_fil, /create
  return
end
