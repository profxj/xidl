;+ 
; NAME:
; sdss_fndqso   
;    Version 1.1
;
; PURPOSE:
;  Return the indices of the QAL structure corresponding to fiber,id
;
; CALLING SEQUENCE:
; indx = sdss_fndqso(qalstr, platfib, /ALL)
;
; INPUTS:
;   fibid -- 2 element vector of plate,fiberid
;
; RETURNS:
;   indices in qalstr matching indx
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /ALL -- Return all QSO's with the same ra,dec and redshift
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;  indx = sdss_fndqso(qalstr, [424,124], /all)
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   13-Oct-2004 Written by JXP
;-
;------------------------------------------------------------------------------

function sdss_fndqso, qalstr, platfib, ALL=all

  if  N_params() LT 2 then begin 
      print,'Syntax - ' + $
        'indx = sdss_fndqso, qalstr, platfib, /ALL [v1.1]'
      return, -1
  endif 
  

  ;; Parse the name
  indx = where(qalstr.plate EQ platfib[0] AND qalstr.fiberid EQ platfib[1],ni)
  if ni NE 1 then stop

  ;; ALL?
  if keyword_set( ALL ) then begin
      all=where(abs(qalstr.ra-qalstr[indx[0]].ra) le 1E-3 AND $
                abs(qalstr.dec-qalstr[indx[0]].dec) le 1E-3 AND $
                abs(qalstr.z_qso-qalstr[indx[0]].z_qso) le .1,nall) 
      indx = all
  endif

  return, indx
end
  
