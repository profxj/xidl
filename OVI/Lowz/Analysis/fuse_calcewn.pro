;+ 
; NAME:
; fuse_calcewn
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
;   fuse_calcewn, struct, fil_instr
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Sep-2003 Added metallicity sturcture
;-
;------------------------------------------------------------------------------
pro fuse_calcewn, strct_fil, instr_list

; lowzovi_prsdat -- Reads in DLA data to a structure

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'fuse_calcewn, files, strct, (v1.0)' 
    return
  endif 

;  Read instrument file list
  readcol, instr_list, instr_fil, inst_dw, inst_w0, format='a,f,f'

  nlist = n_elements(instr_fil)
  if nlist LT 7 then stop

;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)

;  Loop on instruments

  for qq=0L,nlist-1 do begin
      ;; Find lines
      lin = where(strct.instr MOD 2^(qq+1) GT (2^qq-1), nlin)
      if nlin NE 0 then begin
          print, 'fuse_calcen: Reading ', instr_fil[qq]
          ;; Open data
          if qq LE 6 then $
            fx = x_readspec(instr_fil[qq], SIG=sig, wav=wave, NPIX=npix, inflg=3)$
          else begin  ; STIS
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
          ;; Find pixmnx
          mn = min(abs(strct[lin[ii]].wv_lim[0] - wave), pmin)
          mx = min(abs(strct[lin[ii]].wv_lim[1] - wave), pmax)

          ;; Calculate EW
          ew = total( (1. - fx[pmin:pmax])*dwv[pmin:pmax] )
          sigew = sqrt(total( (sig[pmin:pmax] * dwv[pmin:pmax])^2 ))

          ;; Redshift
          if inst_dw[qq] EQ 0. then begin
              wgtwv = total( (1.-(fx[pmin:pmax]<1.))*wave[pmin:pmax] ) / $
                total(1.-(fx[pmin:pmax]<1.))
              strct[lin[ii]].zabs = wgtwv / strct[lin[ii]].wrest - 1.
          endif

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
          0:  ; No value (shouldnt get here, I suspect)
          1: begin ; One value
              strct[ii].EW[0] = strct[ii].EW[a+1]
              strct[ii].sigEW[0] = strct[ii].sigEW[a+1]
          end
          else: begin ; Weighted mean
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
  if not keyword_set( OUTFIL ) then mwrfits, strct, strct_fil, /create

  return
end
