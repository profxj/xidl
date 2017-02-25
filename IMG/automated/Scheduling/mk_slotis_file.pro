;+
; NAME:
;   mk_slotis_file
;
; PURPOSE:
;   Create a queue file for superLOTIS telescope
;
; CALLING SEQUENCE:
;   mk_slotis_file, 'path/to/list', 'DDMMMYYYY', NEXP=,OUTFIL=
;
; INPUTS:
;   list -- Formatted ASCII file containing target and observing info
;           in the following format 
;           Name RR:RR:RR.R -DD:DD:DD.D Filter ExposureTime
;
;   date -- String DDMMMYYYY
;
; OPTIONAL INPUTS:
;                
; OUTPUTS: 
;  Returns a structure describing the instrument configuration which
;  is used to guide the reduction steps.
;
; OPTIONAL INPUT KEYWORDS:
;   NEXP -- The keyword allows the user to specify the number of exposures
;           for each target. By default, three exposures are assumed.
;
; OUTFIL -- The keyword allows the user to specify the name of the output
;           file. If not specified the output file name is 
;           DDMMMYYY_listfilename.txt
;
; OPTIONAL OUTPUTS:
;   
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   10-Dec-2009  Written by JXP
;   05-Jun-2012  Revised by MF
;-
;------------------------------------------------------------------------------
;; 
pro mk_slotis_file, list, date, NEXP=nexp, OUTFIL=outfil

  if  N_params() LT 2  then begin 
    print,'Syntax - ' + $
             'mk_slotis_file, list, date (DDMMMYYYY), NEXP=, OUTFIL= [v1.0]'
    return
  endif 

  print, 'mk_slotis_file:  I hope you have already checked the weather at this website'
  print, 'http://...'

  ;; Initialize
  if not keyword_set(NEXP) then nexp = 3
  if not keyword_set(OUTDIR) then outdir = 'queue/'
  if not keyword_set(OUTFIL) then begin
      ipos = strpos(list, '/', /reverse_sear) 
      if ipos LT 0 then ipos = 0 else ipos = ipos + 1
      outfil = outdir+date+'_'+strmid(list,ipos)
  endif
  if not keyword_set(HA_OFF) then ha_off = 1.

  ;; Check the date format
  date = strtrim(date,2)
  if strlen(date) NE 9 then begin
      print, 'mk_slotis_file:  Date must have this format DDMMMYYYY, e.g. 02JUN2009'
      return
  endif

  ;; Input targets
  readcol, list, name, ras, decs, filter, texp, format='A,A,A,A,F'
  ntarg = n_elements(texp)
  x_radec, ras, decs, ra_targ, dec_targ

  ;; Grab twilight (actually for the previous night because of the way
                       ;; days are defined in JD)
  observatory, 'kpno', obs_str
  twilight = x_calctwilight(date, obs_str, STMID=stmid) ;; Decimal hours
      

  ;; Vette targets
  hmin = twilight[0]-HA_OFF
  if hmin LT 0. then hmin = hmin + 24
  hmax = twilight[1] + HA_OFF
  if hmax GT 24. then hmax = hmax - 24

  if hmin GT hmax then begin
      flg_t = -1
      gdt = where(NOT (ra_targ/15. LE hmin AND ra_targ/15 GE hmax) AND $
                  dec_targ GT -15., ngood)
  endif else begin
      flg_t = 1
      gdt = where(ra_targ/15. GE hmin AND ra_targ/15 LE hmax  AND $
                  dec_targ GT -15., ngood)
  endelse

  ;; Create file
  if ngood EQ 0 then begin
      print, 'mk_slotis_file:  No targets up in the sky'
      print, 'mk_slotis_file:  Try again later'
      return
  endif
  close, /all
  openw, 11, outfil

  for qq=0L,ngood-1 do begin
      idx = gdt[qq]
      ;; Time for observing
      time = ra_targ[idx]/15. - STMID ;; hrs

      if time LT 0 then time = 24+time 
      if ra_targ[idx]/15. LT twilight[0] and $
        ra_targ[idx]/15. GT twilight[1] then begin
          if time LT 12. then time = twilight[1]-STMID $
          else time = twilight[0]-STMID 
      endif
      
      ;; Conver to strings
      hr_obs = long(time)
      s_hr = strtrim(hr_obs,2)
      if hr_obs LT 10 then s_hr = '0'+s_hr
      min_obs = long( (time-long(time))*60 )
      s_min = strtrim(min_obs,2)
      if min_obs LT 10 then s_min = '0'+s_min

      ;; Start the line
      lin = s_hr+' '+s_min
      lin = lin+' '

      ;; RA
      lin = lin + strjoin( strsplit(ras[idx],':',/extrac) )
      lin = lin+' '

      ;; DEC
      lin = lin + strjoin( strsplit(decs[idx],':',/extrac) )
      lin = lin+' '

      ;; EPOCH
      lin = lin+'2000 '

      ;; Exp time
      lin = lin+strtrim(round(texp[idx]),2)
      lin = lin+' '

      ;; Nexp
      lin = lin+strtrim(nexp,2)
      lin = lin+' '

      ;; Filter
      lin = lin+filter[idx]
      lin = lin+' : '

      ;; Name
      lin = lin+name[idx]
      printf, 11, lin
  endfor
  close, /all

  return
end
