;+ 
; NAME:
; fuse_calcion
;  V1.1
;
; PURPOSE:
;   For FUSE observations, calculate ionic column densities for all
;   ions in the FUSE structure.  The code allows for limits and does a
;   weighted mean for multiple transitions.
;
; CALLING SEQUENCE:
;   fuse_calcion, strct_fil, zabs, ion_fil, NHI=, HAND=
;
; INPUTS:
;   strct_fil -- FITS file for the FUSE structure
;   zabs -- Absorption redshift of the system
;
; RETURNS:
;
; OUTPUTS:
;   ion_fil -- Output ion file
;
; OPTIONAL KEYWORDS:
;  NHI=  -- NHI value and error [required at present!]
;  HAND= -- Input file of ionic column densities used to override the
;           values that would otherwise be derived by this program.
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   fuse_calcion, strct_fil, zabs, 'PKS0405_z495.ion'
;
; PROCEDURES CALLED:
;  getabnd
;  getion
;
; REVISION HISTORY:
;   11-Sep-2003 Written by JXP
;-
;------------------------------------------------------------------------------
pro fuse_calcion, strct_fil, zabs, ion_fil, NHI=nhi, HAND=hand

  if (N_params() LT 3) then begin 
    print,'Syntax - ' + $
             'fuse_calcion, strct, zabs, ion_fil, NHI=, HAND= (v1.1)' 
    return
  endif 

  if not keyword_set( NSIG ) then nsig = 3.

  if keyword_set( HAND ) then begin
      print, 'fuse_calcion:  Overriding values using the file in '+$
        strtrim(hand,2)+' !!!'
      readcol, hand, h_za, h_z, h_i, h_f, h_N, h_s, FORMAT='F,I,I,I,F,F'
  endif

; Open ionfil for writing
  close, /all
  openw, 11, ion_fil

  if not keyword_set( NHI ) then stop $
  else printf, 11, ' 1  1  ', NHI[0], NHI[1], '1', $
    FORMAT='(a6,f7.4,1x,f7.4,2x,a1)'
  
;  Open structure
  strct = xmrdfits(strct_fil, 1, /silent)

  ;; Grab key lines
;  a = where(abs(strct.zabs-zabs) LT 0.0002 AND $
  a = where(abs(strct.zabs-zabs) LT 0.0005 AND $
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

          ;; By hand
          flg_h = 0
          if keyword_set( HAND ) then begin
              for kk = 0L,n_elements(h_za)-1 do begin
                  a = where(abs(h_za - zabs) LT 0.0002 AND $
                            h_Z EQ nz AND h_i EQ nion, na)
                  case na of 
                      0: 
                      1: begin
                          cnt = cnt + 1
                          fflg[cnt] = h_f[a]
                          fcolm[cnt] = h_N[a]
                          fsig[cnt] = h_s[a]
                          flg_h = 1
                      end
                      else: stop
                  endcase
                  if flg_h EQ 1 then break
              endfor
          endif
          ;; Lines
          glin = where(atmval EQ nZ AND ionval EQ nion $
                       AND strct.flg MOD 2 EQ 1, ngd)
;          if nz EQ 14 then stop

          if flg_h NE 1 then begin
          if ngd EQ 0 then continue else cnt = cnt+1
          case ngd of
              1:  begin  ; Set value
                  ;; Check n sigma
                  if (strct[glin].Ncolm LT nsig*strct[glin].sigNcolm OR  $
                  strct[glin].ew[0] LT nsig*strct[glin].sigew[0]) $
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
      endif
          printf, 11, nZ, nion, fcolm[cnt], fsig[cnt], fflg[cnt], $
            FORMAT='(i2,1x,i2,1x, f7.4,1x,f7.4,1x,i2)'
          print, nZ, nion, fcolm[cnt], fsig[cnt], fflg[cnt], $
            FORMAT='(i2,1x,i2,1x, f7.4,1x,f7.4,1x,i2)'
      endfor
  endfor

  close, /all
  

  return
end
