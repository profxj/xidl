;+ 
; NAME:
; fuse_calcion
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
;   11-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_calcion, strct_fil, zabs, ion_fil, NHI=nhi

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'fuse_calccolm, strct, zabs, ion_fil, NHI= (v1.0)' 
    return
  endif 

  if not keyword_set( NSIG ) then nsig = 3.

; Open ionfil for writing
  close, /all
  openw, 11, ion_fil

  if not keyword_set( NHI ) then stop $
  else printf, 11, ' 1  1  ', NHI[0], NHI[1], '1', $
    FORMAT='(a6,f7.4,1x,f7.4,2x,a1)'
  
;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)

  ;; Grab key lines
  a = where(abs(strct.zabs-zabs) LT 0.0002 AND $
           strct.flg NE 0, nsys)
      
  strct = strct[a]

  ;; Convert all N back to linear  
  b = where(strct.Ncolm LT -1., nb, COMPLEMENT=A)
  strct[a].Ncolm = 10^strct[a].Ncolm
  strct[a].sigNcolm = strct[a].sigNcolm * alog(10.0) * strct[a].Ncolm
  if nb NE 0 then $  ; Only the error was put in log
    strct[b].sigNcolm = strct[b].sigNcolm * alog(10.0) * abs(strct[b].Ncolm)

  ;; Get Z, ion values
;  nm = strtrim(strmid(strtrim(strct.ion,2),0,2),2)
  atmval = intarr(nsys)
  ionval = intarr(nsys)

  for i=0L,nsys-1 do begin
      getion, strct[i].wrest, ion, celm
      getabnd, celm, atm, abnd
      atmval[i] = atm
      ionval[i] = ion
  endfor

  ;; Create final arrays
  dum = 10*atmval + ionval
  b = dum[ uniq(dum, sort(dum)) ]
  nions = n_elements(b)

  fcolm = fltarr(nions)
  fsig  = fltarr(nions)
  fflg  = intarr(nions)

; LOOP

  cnt = 0L
  for nZ=2L,100 do begin
      for nion = 1L, 10L do begin

          ;; Lines
          glin = where(atmval EQ nZ AND ionval EQ nion $
                       AND strct.flg MOD 2 EQ 1, ngd)
          if ngd EQ 0 then continue else cnt = cnt+1

          case ngd of
              1:  begin  ; Set value
                  ;; Check n sigma
                  if strct[glin].Ncolm LT nsig*strct[glin].sigNcolm $
                    AND strct[glin].flg NE 3 then begin
                      ;; Upper limit
                      fcolm[cnt] = alog10(nsig*strct[glin].sigNcolm)
                      fsig[cnt] = strct[glin].sigNcolm / alog(10) / 10^fcolm[cnt]
                      fflg[cnt] = 5
                  endif else begin
                      fcolm[cnt] = alog10(strct[glin].Ncolm)
                      fsig[cnt] = strct[glin].sigNcolm / alog(10) $
                        / strct[glin].Ncolm
                      fflg[cnt] = strct[glin].flg
                  endelse
              end
              else: begin  ; Multiple values
                  ;; Any good values?
                  ival = where(strct[glin].flg EQ 1 AND strct[glin].Ncolm $
                               GT nsig*strct[glin].sigNcolm,  nval)
                  case nval of 
                      0: begin ;  Only limits
                          ilow = where(strct[glin].flg EQ 3, nlow)
                          if nlow NE 0 then begin  ;; Lower
                              best = max( strct[glin[ilow]].Ncolm, ibst)
                          endif else begin ;; Upper
                              ilow = where(strct[glin].flg EQ 5, nlow)
                              if nlow EQ 0 then stop
                              best = min( strct[glin[ilow]].Ncolm, ibst)
;                              best = (best > 0.) + $
;                                nsig*strct[glin[ilow[ibst]]].sigNcolm
                          endelse
                          fcolm[cnt] = alog10(best)
                          fsig[cnt] = strct[glin[ilow[ibst]]].sigNcolm / $
                            alog(10) / best
                          fflg[cnt] = strct[glin[ilow[0]]].flg
                      end
                      1:  begin ; One good value
                      ;; Check n sigma
                          if strct[glin[ival]].Ncolm LT $
                            nsig*strct[glin[ival]].sigNcolm then begin
                              fcolm[cnt] = alog10(nsig*strct[glin[ival]].sigNcolm)
                              fsig[cnt] = strct[glin[ival]].sigNcolm $
                                / alog(10) / 10^fcolm[cnt]
                              fflg[cnt] = 5
                          endif else begin
                              fcolm[cnt] = alog10(strct[glin[ival]].Ncolm)
                              fsig[cnt] = strct[glin[ival]].sigNcolm / alog(10) $
                                / strct[glin[ival]].Ncolm
                              fflg[cnt] = strct[glin[ival]].flg
                          endelse
                      end
                      else: begin ; Multiple values
                          ;; Weighted mean
                          twgt = total(1./strct[glin[ival]].sigNcolm^2)
                          weight_N = $
                            total(strct[glin[ival]].Ncolm/ $
                                  strct[glin[ival]].sigNcolm^2) / twgt
                          weight_sN = 1./sqrt(twgt)
                          fcolm[cnt] = alog10(weight_N)
                          fsig[cnt] = weight_sN / alog(10) / weight_N
                          fflg[cnt] = strct[glin[ival[0]]].flg
                      endelse
                  endcase
              end
          endcase
          printf, 11, nZ, nion, fcolm[cnt], fsig[cnt], fflg[cnt], $
            FORMAT='(i2,1x,i2,1x, f7.4,1x,f7.4,1x,i2)'
          print, nZ, nion, fcolm[cnt], fsig[cnt], fflg[cnt], $
            FORMAT='(i2,1x,i2,1x, f7.4,1x,f7.4,1x,i2)'
      endfor
  endfor

  close, /all
  

  return
end
