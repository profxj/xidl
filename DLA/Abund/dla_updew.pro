;+ 
; NAME:
; dla_updew
;  V1.1
;
; PURPOSE:
;    Given a list of DLA .dat files, calculate the EWs and write out
;    .EW files
;
; CALLING SEQUENCE:
;   dla_updew, list
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
;   dla_updew, dla
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   2006 Written by JXP 
;-
;------------------------------------------------------------------------------
pro dla_updew, list, FLG_PLT=flg_plt, CONTI=conti

; parse_dlalst -- Reads in DLA data to a structure

  if (N_params() LT 1) then begin 
    print,'Syntax - ' + $
             'dla_updew, list, /FLG_PLT, CONTI= [v1.1]'
    return
  endif 


  close, /all
  ;; Parse
  parse_dlalst, dla, list, /NOELM
  ndla = n_elements(dla)

  ;; Loop
  for nn=0L,ndla-1 do begin
      print, 'DLA: ', dla[nn].qso, dla[nn].zabs, strtrim(dla[nn].abndfil,2)
      ;; Abundances
      if strlen(strtrim(dla[nn].abndfil,2)) EQ 0 then begin
          print, 'dla_updew:  No Abnd file.  Skipping this one..'
          continue
      endif
      ism_calcew, dla[nn].abndfil, all_cog, ISM=ism, NOCONTI=(not keyword_set(CONTI))

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

  ;; Write
;  dla_writestr, dla
  print, 'dla_updew: All done'

  return
end
