;+ 
; NAME:
; fuse_calcewn
;  V1.1
;
; PURPOSE:
;    Given a FUSE structure file and a set of instrument files (each
;  of which specifies absorption lines), the code measures EW and fills
;  up the FUSE structure.  If multiple EW values are measured, then it
;  calculates the weighted mean if the two are consistent or otherwise
;  the average and sets a large uncertainty.
;
; CALLING SEQUENCE:
;   fuse_calcewn, strct_fil, instr_list
;
; INPUTS:
;  strct_fil  -- IDL FUSE file
;  instr_list -- Instrument list containing absorption line lists
;
; RETURNS:
;
; OUTPUTS:
;  strct_fil  -- The EW tags are filled up by this routine
;
; OPTIONAL KEYWORDS:
;  MODLIM - modify wavelength limits with wavelength shift
;  ROOTDIR - pre-pend string to instr_list file names
;  SUFFIX - for STIS-like files, suffix to search and replace
;           (default: ['f.fits','e.fits'])
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Sep-2003 Added metallicity sturcture
;    2-Feb-2006 added /modlim and /ewgauss options, KLC
;    3-Apr-2006 added calculation of zsig
;   24-Jan-2007 removed /ewgauss and modify EW summing when flux<0 or flux>1
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro fuse_calcewn, strct_fil, instr_list, modlim=modlim, rootdir=rootdir, $
                  suffix=suffix

; lowzovi_prsdat -- Reads in DLA data to a structure

if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
      'fuse_calcewn, files, strct, [v1.1]' 
    return
endif 

;  Read instrument file list
readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'
if keyword_set(rootdir) then instr_fil = rootdir+instr_fil
if not keyword_set(suffix) then $
   suffix = ['f.fits','e.fits'] ; instrument >= 128 

nlist = n_elements(instr_fil)
if nlist LT 7 then stop

;  Open structure
strct = xmrdfits(strct_fil, 1, /silent)
zarr = dblarr(20,2)

;  Loop on instruments

for qq=0L,nlist-1 do begin
    ;; Find lines
    lin = where(strct.instr MOD 2^(qq+1) GT (2^qq-1), nlin)
    if nlin NE 0 then begin
        print, 'fuse_calcewn: Reading ', instr_fil[qq]

        test = file_search(instr_fil[qq],count=ntest)
        if ntest eq 0 then continue ;no such file

        ;; Open data
        if qq LE 6 then $
          fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, NPIX=npix, inflg=3)$
        else begin              ; STIS
            spos = strpos(instr_fil[qq], suffix[0])
            sig_fil = strmid(instr_fil[qq], 0, spos)+suffix[1]
            fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, NPIX=npix, $
                            fil_sig=sig_fil, inflg=0)
        endelse
        ;; Shift wave
        wave = wave + inst_dw[qq]*inst_w0[qq]/wave

          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
        ;;  Modify wv_lim to account for shifts, added by KLC          
        if keyword_set(modlim) then begin
            strct[lin].wv_lim[0] = strct[lin].wv_lim[0] + $
              inst_dw[qq]*inst_w0[qq]/strct[lin].wv_lim[0]
            strct[lin].wv_lim[1] = strct[lin].wv_lim[1] + $
              inst_dw[qq]*inst_w0[qq]/strct[lin].wv_lim[1]
        endif 
          ;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;


        ;; Sort?
        srt = sort(wave)
        wave = wave[srt]
        fx = fx[srt]
        sig = sig[srt]
        ;; dwv
        dwv = wave - shift(wave,1)
        dwv[0] = dwv[1]
    endif

    for ii=0L,nlin-1 do begin
        ;; Find pixmnx
        mn = min(abs(strct[lin[ii]].wv_lim[0] - wave), pmin)
        mx = min(abs(strct[lin[ii]].wv_lim[1] - wave), pmax)

        ;; Calculate EW
        ;; Limit to 'expected' flux range -sigma < flux < 1+sigma, KLC
        rng = lindgen(pmax-pmin+1) + pmin
        nfx = where(fx[rng] lt -abs(sig[rng]),nnfx)
        pfx = where(fx[rng] gt 1.+abs(sig[rng]),npfx)

        wgt = 1-fx[rng]
        if nnfx ne 0 then wgt[nfx] = 1+2*abs(sig[rng[nfx]])
        if npfx ne 0 then wgt[pfx] = -2*abs(sig[rng[pfx]])

        ew = total( wgt*dwv[rng] )
        sigew = sqrt(total( (sig[rng] * dwv[rng])^2 ))
        
        ;; Redshift
        gd = where(fx[rng] gt 0. and fx[rng] le 1.,ngd)
        if ngd eq 0 then begin 
            ;; Flux-weighted
            wgtwv = total((wgt>0)*wave[rng])/total((wgt>0))        
            strct[lin[ii]].zabs = wgtwv / strct[lin[ii]].wrest - 1.
            strct[lin[ii]].zsig = mean(dwv[rng])/wgtwv
        endif else begin
            ;;"Optical depth-weighted", KLC
            wgt = replicate(0d,pmax-pmin+1)
            wgt[gd] = alog(1./fx[rng[gd]])
            bd=where(fx[rng] le 0,nbd)
            if nbd ne 0 then wgt[bd] = 1e6*max(wgt[gd]) ;some really large num
            wgtwv = total(wgt*wave[rng])/total(wgt)
            strct[lin[ii]].zabs = wgtwv / strct[lin[ii]].wrest - 1.
            strct[lin[ii]].zsig = mean(dwv[rng])/wgtwv
        endelse 

        ;;Flux-weighted -- use wgt from above
;        wgtwv = total((wgt>0)*wave[rng])/total((wgt>0))        
;;        wgtwv = total( (1.-(fx[pmin:pmax]<1.))*wave[pmin:pmax] ) / $
;;          total(1.-(fx[pmin:pmax]<1.))
;        strct[lin[ii]].zabs = wgtwv / strct[lin[ii]].wrest - 1.
;        strct[lin[ii]].zsig = mean(dwv[pmin:pmax])/wgtwv
;;        strct[lin[ii]].zsig = (strct[lin[ii]].zabs+1)*$
;;          sqrt(total(sig[rng]^2*wave[rng]^2)/total(wgt)^2 + $
;;               total(sig[rng])^2/total(wgt)^2)

        zarr[qq+1,*] = [strct[lin[ii]].zabs,strct[lin[ii]].zsig]
        if not finite(strct[lin[ii]].zabs) then stop


        ;; Convert to rest EW
        strct[lin[ii]].EW[qq+1] = ew / (1.+(strct[lin[ii]].zabs > 0.))
        strct[lin[ii]].sigEW[qq+1] = sigew / (1.+(strct[lin[ii]].zabs > 0.))
;          if abs(strct[lin[ii]].wrest - 1031.9261) LT 0.001 and $
;            strct[lin[ii]].zabs GT 0.361 then stop
    endfor
endfor

;; Combine EW
nlin = n_elements(strct)
for ii=0L,nlin-1 do begin
    a = where(strct[ii].sigEW[1:*] GT 0., na)
    case na of
        0:                   ; No value (shouldnt get here, I suspect)
        1: begin                ; One value
            strct[ii].EW[0] = strct[ii].EW[a+1]
            strct[ii].sigEW[0] = strct[ii].sigEW[a+1]
        end
        else: begin             ; Weighted mean
            twgt = total(1./strct[ii].sigEW[a+1]^2)
            strct[ii].EW[0] = $
              total(strct[ii].EW[a+1]/strct[ii].sigEW[a+1]^2) / twgt
            strct[ii].sigEW[0] = 1./sqrt(twgt)
            ;; Check chi2 for divergent lines
            chi2 = total( (strct[ii].EW[0]-strct[ii].EW[a+1])^2 / $
                          strct[ii].sigEW[a+1]^2 ) 
            if chisqr_pdf(chi2, na-1) GT 0.9 then $
              print, 'Trouble with: ', strct[ii].wrest, strct[ii].zabs
        end
    endcase  
endfor

;; mA
strct.EW = strct.EW * 1000
strct.sigEW = strct.sigEW * 1000
writecol, 'tmp.dat', strct.wrest, strct.zabs, strct.EW[0], strct.sigEW[0],  $
  strct.EW[4], strct.EW[6], strct.EW[7], FMT='(f7.2,1x,f9.5,1x, 5f8.2)'

;; Overwrite
mwrfits, strct, strct_fil, /create

return
end
