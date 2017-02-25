;+ 
; NAME:
; fuse_calccolm
;  V1.1
;
; PURPOSE:
;    Given a list of DLA base files, fill up the structure ;
; CALLING SEQUENCE:
;  fuse_calccolm, strct_fil, instr_list, abs_list, NSIG=nsig
;
; INPUTS:
;
; RETURNS:
;  strct_fil  -- FUSE absorption line structure
;  instr_list -- List of instruments (expected to be the 
;                FUSE channels than STIS)
;  abs_list   -- List of files guiding the column density measurments
;                (formatting is key here)
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  NSIG=  -- Number of sigma for signficance [default: 3.]
;  ROOTDIR - pre-pend string to instr_list file names
;  SUFFIX - for STIS-like files, suffix to search and replace
;           (default: ['f.fits','e.fits'])
;;
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
;    6-Jun-2007 Sort lines in structure (to match with *.clm), KLC
;-
;------------------------------------------------------------------------------
pro fuse_calccolm, strct_fil, instr_list, abs_list, NSIG=nsig,  $
                   ROOTDIR=rootdir, SUFFIX=suffix

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'fuse_calccolm, strct, instr_list, abs_list, NSIG=  [v1.1]' 
    return
  endif 

  if not keyword_set( NSIG ) then nsig = 3.
  if not keyword_set( SUFFIX ) then suffix = ['f.fits','e.fits'] 

;  Read instrument file list
  readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'
  if keyword_set(rootdir) then instr_fil = rootdir+instr_fil

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
      if nhand NE 0 then begin
          hnd_i = lonarr(100,3)
          hnd_N = fltarr(100)
          hnd_s = fltarr(100)
          i1 = 0
          i2 = 0
          i3 = 0
          f1 = 0.
          f2 = 0.
          for jj=0L,nhand-1 do begin
              readf, 1, i1, i2, i3, f1, f2
              hnd_i[jj,0] = i1
              hnd_i[jj,1] = i2
              hnd_i[jj,2] = i3
              hnd_N[jj] = f1
              hnd_s[jj] = f2
          endfor
      endif
      
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
          a = where(abs(strct.zabs-zabs) LT 0.002 AND $
                    abs(strct.wrest-wv_lin[i]) LT 0.003, na)
;          if na NE 1 then begin
;              print, 'na', na
;              print, zabs, wv_lin[i]
;              stop
;          endif
          if na eq 0 then stop,'fuse_calccolm: no matching line',zabs,wv_lin[i]
          if na gt 1 then begin
             dum = min(strct[a].zabs-zabs,imn,/absolute)
             a = a[imn]
          endif 
          msk[a] = a[0]
          flg_lin[i] = dumi
      endfor                    ;end read lines from *.clm
      
      ;; Create colm arrays
      all_N = dblarr(ncalc, nlist)
      all_sN = dblarr(ncalc, nlist)
      
      ;; Subset of structure
      gd = where(msk GE 0L)
      srt = sort(strct[gd].wrest) ;assumes abs_lin ordered by wrest
      gd = gd[srt]              ;to match abs_lin
      
;  Loop on instruments
      for qq=0L,nlist-1 do begin
          ;; Find lines
          lin = where(strct[gd].instr MOD 2^(qq+1) GT (2^qq-1), nlin)

          ;; Test
          test = file_search(instr_fil[qq],count=ntest)
          if ntest eq 0 then continue

          ;; Grab the spectrum
          if nlin NE 0 then begin
              print, 'fuse_calcen: Reading ', instr_fil[qq]
              ;; Open data
              if qq LE 6 then $
                fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, $
                                NPIX=npix, inflg=3)$
              else begin        ; STIS
                  spos = strpos(instr_fil[qq], suffix[0])
                  sig_fil = strmid(instr_fil[qq], 0, spos)+suffix[1]
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
              
              ;; Limit to 'expected' flux range -sigma < flux <
              ;; 1+sigma, KLC
              ;; Similar to fuse_calcewn
              rng = lindgen(pmax-pmin+1) + pmin
              nfx = where(fx[rng] lt -abs(sig[rng]),nnfx)
              pfx = where(fx[rng] gt 1.+abs(sig[rng]),npfx)

              nwfx = fx[rng]
              if nnfx ne 0 then nwfx[nfx] = -2*abs(sig[rng[nfx]])
              if npfx ne 0 then nwfx[pfx] = 1+2*abs(sig[rng[pfx]])

              ;; Calculate N
;              x_aodm, wave[pmin:pmax], fx[pmin:pmax], sig[pmin:pmax],$
              x_aodm, wave[pmin:pmax], nwfx, sig[pmin:pmax], $
                strct[jj].wrest, clm, sig_clm
              
              all_N[lin[ii],qq] = clm
              all_sN[lin[ii],qq] = sig_clm
          endfor                ;end measure ncolm for lines from same instr
      endfor                    ;end instr loop
      
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
              if strct[jj].flg MOD 8 LT 3 then begin
                  strct[jj].flg = strct[jj].flg + 4
                  print,'fuse_calccolm: upper limit now for ',$
                    strct[jj].zabs,strct[jj].wrest
              endif 
          endif

          ;; Log
          if weight_N GT 0. then colm = alog10( weight_N ) else colm = weight_N
          lgvar = ((1.d / (alog(10.0)*weight_N))^2)*(weight_sN^2)
          sig_colm = sqrt(lgvar)

          ;; Fill up
          strct[jj].Ncolm = colm
          strct[jj].sigNcolm = sig_colm

          ;; By hand!
          if nhand NE 0 then begin
              getion, strct[jj].wrest, ion, Z=zval
              for vv=0L,nhand-1 do begin
                  if hnd_i[vv,1] EQ ion AND hnd_i[vv,0] EQ zval then begin
                      strct[jj].Ncolm = hnd_N[vv]
                      strct[jj].sigNcolm = hnd_s[vv]
                  endif  
              endfor
          endif

          ;; Check flags
          if strct[jj].flg NE flg_lin[ii] AND flg_lin[ii] NE 1  AND $
            strct[jj].flg MOD 2 EQ 1 then begin
              print, 'fuse_calccolm: The flags dont match!', strct[jj].wrest
              print, flg_lin[ii], strct[jj].flg
              stop 
          endif
      endfor                    ;end loop line (combining N)
      printcol, strct[gd].ion, strct[gd].Ncolm, strct[gd].sigNcolm, strct[gd].flg
  endfor                        ;end current *.clm (abs_fil[nn])

  close, /all
      
      ;; Overwrite
      if not keyword_set( OUTFIL ) and NOT keyword_set( NOSV ) then $
        mwrfits, strct, strct_fil, /create
  return
end
