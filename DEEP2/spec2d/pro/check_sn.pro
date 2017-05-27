;+
; NAME:
;   CHECK_SN
; 
; PURPOSE:
;   Estimate remaining number of frames required for a mask
;
; CATEGORY:
;   Quicklook
;
; CALLING SEQUENCE:
;   check_sn [, mask, nomsig=nomsig]
;
;
; OPTIONAL INPUTS:
;   mask - the full name of a DEIMOS mask (string), e.g. '2251.W'
;          If mask is not specified, the most recent mask analyzed by
;          quicklook will be checked.
;
; KEYWORD PARAMETERS:
;   nomsig - the desired signal-to-noise per extracted [1d]
;            pixel, derived from the find_object profile.  Default
;            0.9, which corresponds to a fairly bad mask (75% of masks
;            are better than this).  We might consider doing 1.0 (60%
;            of masks are better than that).
;
; RESTRICTIONS:
;    Must be run in current night's directory where quicklook is
;    running, so that science_qa.dat will be available.
; 
;
; MODIFICATION HISTORY:
;  2003jun25 JAN
;-


pro check_sn, mask, nomsig=nomsig

     if n_elements(nomsig) eq 0 then nomsig = 0.9

     fileexists =  file_test('science_qa.dat')
     
      if fileexists then readcol,  'science_qa.dat',  filename,  masks,  align,  pm1,  $
         alignsig,  seeing,  pm2,  seeingsigma,  s2narr,  pm3,  s2nsigma,  $
        skipline=1,  format='A,A,F,A,F,F,A,F,F,A,F',  /SILENT

      nfiles = n_elements(masks)
      if n_elements(mask) eq 0 then lastmask = masks[nfiles-1] else lastmask = mask

      whmask = where(masks eq lastmask,  nmask)
      print
      print, 'For mask: ', lastmask,  '  there are ', $
        string(nmask, form='(i3)'), '  files that have been processed'

      whbad = where(finite(s2nsigma) eq 0, badct)
      if badct gt 0 then s2nsigma[whbad] = 0.2
      s2ntot = sqrt(total(s2narr[whmask]^2))
      s2nerr =  sqrt( total(s2narr[whmask]^2*s2nsigma[whmask]^2)/total(s2narr[whmask]^2) )

      print, 'S/N so far:      ', string(s2ntot, form='(f6.3)'), ' +/- ', string(s2nerr, form='(f6.3)')

      s2nave = s2ntot/sqrt(nmask)
      aveerr = s2nerr/sqrt(nmask)

      print, 'S/N for 3 frames:', string(s2nave*sqrt(3.),  form='(f6.3)'), ' +/- ', $
        string(aveerr*sqrt(3.),  form='(f6.3)')
      print, 'Nominal: ', string(nomsig,  form='(f6.3)')


      nreq = (nomsig/s2nave)^2
      errreq = 2.*nreq*aveerr/s2nave
      nleft = nreq-nmask

      print, 'It is advisable to take ', string(nleft > 0,  form='(f6.2)'), $
        ' +/- ', string(errreq,  form='(f6.2)'), ' more frames like the average of these'
      print


return
end
