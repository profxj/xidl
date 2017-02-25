;+ 
; NAME:
; x_fineabnd
;   Version 1.1
;
; PURPOSE:
;    Given an abundance structure, identify the fine-structure lines,
;    and get 'final' values
;
; CALLING SEQUENCE:
;   x_fineabnd, strct
;
; INPUTS:
;   strct - Abundance structure,  i.e. {ionstruct}
;
; RETURNS:
;
; OUTPUTS:
;  Fills up the 'final' values and limits for fine structure transitions
;
; OPTIONAL KEYWORDS:
;  FINEDAT= -- File with the fine structure data 
;  /NOLOG   -- Columns in the DLA structure are not logarithmic 
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES CALLED:
;  dla_allclm
;
; REVISION HISTORY:
;   2-Dec-2005 Written by JXP
;-
;------------------------------------------------------------------------------

pro x_fineabnd, strct, FINEDAT=finedat, NOLOG=nolog

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'x_fineabnd, strct, FINEDAT= [V1.1]'
    return
  endif 

  ;; Read fine_strct data
  if not keyword_set( FINEDAT ) then $
    finedat = getenv('XIDL_DIR')+'/Spec/Lines/Lists/fine_strct.lst'
  readcol, finedat, fin_Z, fin_ion, fin_j, fin_wav, fin_E, fin_A, $
           FORMAT='I,I,F,D,F,F', /silent

  if keyword_set(NOLOG) then LOG = 0 else LOG = 1

  u_Z = fin_Z[ uniq(fin_Z,sort(fin_Z))]

  ;; Loop on the structure
  for nn=0L,n_elements(strct)-1 do begin
      ;; Loop on Z
      for qq=0L,n_elements(u_Z)-1 do begin
          ;; Loop on Ion
          gdz = where(fin_Z EQ u_Z[qq])
          u_i = fin_ion[gdz[ uniq(fin_ion[gdz],sort(fin_ion[gdz]))]]
          i = u_Z[qq]
          for kk=0L,n_elements(u_i)-1 do begin
              j = u_i[kk]
              gdi = where(fin_Z EQ u_Z[qq] AND fin_ion EQ u_i[kk], nwv)
              ;; Loop on Energy!
              u_E = fin_E[ gdi[uniq(fin_E[gdi],sort(fin_E[gdi]))]]
              u_E = u_E[sort(u_E)]
              u_E = [0., u_E]

              ;; Loop on excited states
              for mm=1L,n_elements(u_E)-1 do begin
                  gd = where(fin_Z EQ u_Z[qq] AND fin_ion EQ u_i[kk] AND $
                              fin_E EQ u_E[mm], nwv)

                  fine_ion = 20+j*10+mm
                  ionflg=0
                  fidx = [-1L]
                  for jj=0L,nwv-1 do begin
                      mt = where(abs((strct[nn].ion[i].state[j,*].lambda)[*] $
                                     - fin_wav[gd[jj]]) LT 0.01, nmt)
                      if nmt EQ 1 then fidx = [fidx,mt[0]]
                  endfor
                  nobs = n_elements(fidx) - 1
                  if nobs GE 1 then fidx = fidx[1:*] else continue

                  ;; Save
                  dims = size(strct[nn].ion[i].state[fine_ion,1:nobs].lambda,/dim)
                  if n_elements(dims) EQ 1 then dims = [1,1]
                  strct[nn].ion[i].state[fine_ion,1:nobs].lambda = $
                    reform( (strct[nn].ion[i].state[j,*].lambda)[fidx], dims)
                  strct[nn].ion[i].state[fine_ion,1:nobs].clm = $
                    reform( (strct[nn].ion[i].state[j,*].clm)[fidx] , dims)
                  strct[nn].ion[i].state[fine_ion,1:nobs].sigclm = $
                    reform( (strct[nn].ion[i].state[j,*].sigclm)[fidx] , dims)
                  strct[nn].ion[i].state[fine_ion,1:nobs].flgclm = $
                    reform( (strct[nn].ion[i].state[j,*].flgclm)[fidx] , dims)

                  ;; Normal
                  for jj=0L,nobs-1 do begin
                      k = fidx[jj]
                      if strct[nn].ion[i].state[j,k].flgclm LT 2 OR $
                        strct[nn].ion[i].state[j,k].flgclm EQ 9 OR $
                        strct[nn].ion[i].state[j,k].flgclm EQ 8 then begin
                          ionflg = 1
                          dla_allclm, strct, nn, i, fine_ion, $
                                      strct[nn].ion[i].state[j,k].clm, $
                                      strct[nn].ion[i].state[j,k].sigclm,$
                                      LOG=LOG
                          strct[nn].ion[i].state[fine_ion].flgclm  = 1
                      endif
                  endfor
                  ;; Saturated
                  if ionflg NE 1 then begin
                      for jj=0,nobs-1 do begin
                          k = fidx[jj]
                          if strct[nn].ion[i].state[j,k].flgclm EQ 3 OR $
                            strct[nn].ion[i].state[j,k].flgclm EQ 2 OR $
                            strct[nn].ion[i].state[j,k].flgclm EQ 10 OR $
                            strct[nn].ion[i].state[j,k].flgclm EQ 11 then begin
                              ionflg = 2
                              strct[nn].ion[i].state[fine_ion].flgclm  = 2
                              dla_allclm, strct, nn, i, fine_ion, $
                                          strct[nn].ion[i].state[j,k].clm, $
                                          strct[nn].ion[i].state[j,k].sigclm,$
                                          LOG=LOG
                          endif
                      endfor
                  endif
                  ;; Upperlimit
                  if ionflg EQ 0 then begin
                      for jj=0,nobs-1 do begin
                          k = fidx[jj]
                          if strct[nn].ion[i].state[j,k].flgclm EQ 4 OR $
                            strct[nn].ion[i].state[j,k].flgclm EQ 5 OR $
                            strct[nn].ion[i].state[j,k].flgclm EQ 12 OR $
                            strct[nn].ion[i].state[j,k].flgclm EQ 13 then begin
                              ionflg = 3
                              ;;  Allow for blends for upper limits
                              if strct[nn].ion[i].state[fine_ion].clm LE 0. then $
                                strct[nn].ion[i].state[fine_ion].clm = $
                                  10^strct[nn].ion[i].state[j,k].clm $
                              else $
                                strct[nn].ion[i].state[fine_ion].clm = $
                                  10^strct[nn].ion[i].state[j,k].clm < $
                                  strct[nn].ion[i].state[fine_ion].clm 
                              ;;
                              strct[nn].ion[i].state[fine_ion].flgclm  = 3
                          endif
                      endfor
                  endif
              endfor
          endfor
      endfor
  endfor
              
  return
end

