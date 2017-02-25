;+ 
; NAME:
; ism_sumvpfit
;  V1.0
;
; PURPOSE:
;    Parses a VPFIT output file and sums up individual components
;
; CALLING SEQUENCE:
;   ism_sumvpfit, vpfil, ionst
;
; INPUTS:
;  vpfil --  Output file from VPFIT
;
; RETURNS:
;
; OUTPUTS:
;  ionstr -- Ion structure.  Column densities have linear values
;
; OPTIONAL KEYWORDS:
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2007 Written by JXP 
;-
;------------------------------------------------------------------------------
pro ism_sumvpfit, vpfil, ionstr

; parse_ismlst -- Reads in DLA data to a structure

  if (N_params() LT 2) then begin 
    print,'Syntax - ' + $
             'ism_sumvpfit, vpfil, ism [v1.0]'
    return
  endif 

  ionstr = replicate({ionstruct},101)

  ;; Parse the VPFIT file
  x_parse_vpfit, vpfil, strct
  ncomp = n_elements(strct)

  union = strct[uniq(strct.ion,sort(strct.ion))].ion
  nuni = n_elements(union)

  ;; Loop
  for nn=0L,nuni-1 do begin
      ;; Get Z, ion
      getion, union[nn], ion, Z=z, /INM

      ;; Fine structure?
      ipos1 = strpos(union[nn],'*')
      ipos2 = strpos(union[nn],'*',/reverse_search)
      if ipos1 NE -1 then nstr = ipos2-ipos1+1 else nstr = 0
      if nstr NE 0 then begin
          ;; CII*
          if strmatch(union[nn],'CII*') then ion = 6 else stop
      endif

      comp = where(strmatch(strct.ion,union[nn]),ncomp)
      x_logclm, Nlin, siglin, strct[comp[0]].N, strct[comp[0]].Nsig, /reverse
      ionstr[Z].state[ion,0].clm = Nlin
      ionstr[Z].state[ion,0].sigclm = siglin

      ;; Add em up
      for jj=1,ncomp-1 do begin
          x_logclm, Ni, siglin, strct[comp[jj]].N, $
                    strct[comp[jj]].Nsig, /reverse

          ionstr[Z].state[ion,0].clm = ionstr[Z].state[ion,0].clm + Ni 
          ionstr[Z].state[ion,0].sigclm = $
            sqrt( ionstr[Z].state[ion,0].sigclm^2 + siglin^2)
      endfor
  endfor

  return
end
