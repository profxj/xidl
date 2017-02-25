;; The MO on GTC-OSIRIS seems to be to take separate arc exposures for
;; different combinations of lamps.  To make wavelength calibrations
;; easier, this script combines the separate exposures into a single
;; image that can then be placed in the plan.par file.

Pro osiris_combine_arc, infiles, outfile

  ;; open first arc and store header info
  a0 = mrdfits(infiles[0], 0, hdr0)
  a1 = mrdfits(infiles[0], 1, hdr1)
  a2 = mrdfits(infiles[0], 2, hdr2)

  a1 = a1*sxpar(hdr1,'BSCALE')+sxpar(hdr1,'BZERO') - 1000 ;; 1000 is approx bias level
  a2 = a2*sxpar(hdr2,'BSCALE')+sxpar(hdr2,'BZERO') - 1000

  ncol = (size(a1, /dim))[0]
  nrow = (size(a1, /dim))[1]
  centercol = ncol/2

  ;; smash a few columns near the middle
  s1 = djs_avsigclip( a1[centercol-20:centercol+20,*], 1)
  s2 = djs_avsigclip( a2[centercol-20:centercol+20,*], 1)

  ;; scale so that coadd peaks around 10000 counts...
  a1 *= 10000/max(s1)/n_elements(infiles)
  a2 *= 10000/max(s2)/n_elements(infiles)

  ;; now look at other infiles
  for i=1, n_elements(infiles)-1 do begin
     b1 = mrdfits(infiles[i], 1, bhdr1)
     b2 = mrdfits(infiles[i], 2, bhdr2)
     b1 = b1*sxpar(bhdr1,'BSCALE')+sxpar(bhdr1,'BZERO')
     b2 = b2*sxpar(bhdr2,'BSCALE')+sxpar(bhdr2,'BZERO')
     s1 = djs_avsigclip( b1[centercol-20:centercol+20,*], 1)
     s2 = djs_avsigclip( b2[centercol-20:centercol+20,*], 1)
     b1 *= 10000/max(s1)/n_elements(infiles)
     b2 *= 10000/max(s2)/n_elements(infiles)
     a1 += b1
     a2 += b2
  endfor
  a1 = (fix(a1)+1000-sxpar(hdr1,'BZERO'))/sxpar(hdr1,'BSCALE')
  a2 = (fix(a2)+1000-sxpar(hdr2,'BZERO'))/sxpar(hdr2,'BSCALE')

  mwrfits, a0, outfile, hdr0, /create
  mwrfits, fix(a1), outfile, hdr1
  mwrfits, fix(a2), outfile, hdr2
end
