;+
; NAME:
;   mmt_slitmask
;
; PURPOSE:
;   Combines multiple arcs.  One generally observes the He separately
;   from the Hg exposures (the latter take longer).  The code requires
;   that you combine these exposures prior to running long_reduce.
;
; CALLING SEQUENCE:
; mmt_waveblue, hgfile, hefile, outfil
;
; INPUTS:
;  hgfile -- Filename of Raw image of the Hg arc exposure
;  hefile -- Filename of Raw image of the He arc exposure
;
; OUTPUTS:
;  outfil -- Output filename for the combined arc
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;   This will break if the binning is not 2x4
;
; PROCEDURES CALLED:
;
; REVISION HISTORY:
;   20-Apr-2005  Written by J. Hennawi Berkeley
;-
;------------------------------------------------------------------------------

pro mmt_waveblue, hgfile, hefile, outfil

  if  N_params() LT 3  then begin 
      print,'Syntax - ' + $
        'mmt_waveblue, hgfile, hefile, outfil [v1.0]' 
      return
  endif 


  long_proc, Hgfile, arc_HgCd, hdr = hdr1, gain = gain
  long_proc, Hefile, arc_HeNeAr, hdr = hdr2, gain = gain
  
  imag = arc_HgCd + arc_HeNeAr
  sz = size(arc_HgCd, /dimens)
  
  arc_fake = fltarr(2708L, sz[0])
  arc_fake[0:2687, *] = transpose(arc_HgCd + arc_HeNeAr)/double(gain)
  
  mwrfits, arc_fake, outfil, hdr1, /create
  spawn, 'gzip -f '+outfil
  long_proc, outfil, imag2, gain = gain

END
