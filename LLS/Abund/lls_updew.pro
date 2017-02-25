;+ 
; NAME:
; lls_updew
;  V1.1
;
; PURPOSE:
;    Given a list of DLA .dat files, calculate the EWs and write out
;    .EW files
;
; CALLING SEQUENCE:
;   lls_updew, list
;
; INPUTS:
;  list -- List of DLA .dat files.
;
; RETURNS:
;
; OUTPUTS:
;  Series of DLA files
;
; OPTIONAL KEYWORDS:
;  CONTI=  -- Estimate associated with the EW error
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   lls_updew, lls
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2006 Written by JXP 
;-
;------------------------------------------------------------------------------
pro lls_updew, list, FLG_PLT=flg_plt, CONTI=conti

; parse_llslst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'lls_updew, list, /FLG_PLT, CONTI= [v1.1]'
    return
  endif 


  close, /all
  ;; Parse
  lls_struct, lls, list
  nlls = n_elements(lls)

  ;; Loop
  for nn=0L,nlls-1 do begin
      print, 'LLS: ', lls[nn].qso, lls[nn].zabs
      ;; Abundances
      for sys=0L,lls[nn].nsys-1 do begin
          ;; Abundances
          ism_calcew, lls[nn].systems[sys].abndfil, all_cog, $
                      ISM=ism, NOCONTI=(not keyword_set(CONTI))

          if not keyword_set(all_cog) then continue
          
          ;; Output
          len = strlen(ism.tab_fil)
          EW_fil = strmid(ism.tab_fil,0,len-3)+'EW'
          close, 13
          openw, 13, EW_fil
          ncog = n_elements(all_cog)
          for ww=0,ncog-1 do begin
              for k=0,all_cog[ww].nlin-1 do begin
                  if all_cog[ww].wrest[k] LT 1. then stop
                  printf, 13,  $
                          all_cog[ww].wrest[k], $
                          all_cog[ww].ew[k], $
                      all_cog[ww].sigew[k], $
                          all_cog[ww].inst[k], $
                          format='(f10.4,1x,2f8.4,1x,i3)'
              endfor
          endfor
          close, 13
      endfor
  endfor

  ;; Write
;  lls_writestr, lls
  print, 'lls_updew: All done'

  return
end
