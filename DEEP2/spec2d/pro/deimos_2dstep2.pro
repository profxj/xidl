;+
; NAME:
;    deimos_2dstep2
;
; PURPOSE:
;    combine the multi-HDU data for a given slit into one 2d and 1d spectrum
;
; CALLING SEQUENCE:
;    deimos_2dstep2, filein
; 
; INPUTS:
;    filein -- FITS file containing first stage reduced data
;
; OPTIONAL INPUTS:
;    something denoting sky-only spectra
;	
; KEYWORDS:
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   the second stage of spec2d needs to use sky-only spectrum, when
;   needed. we need a means to decide whether to use local sky, or
;   special sky.  Not certain yet how to decide which is best.
;
; REVISION HISTORY:
;   MD, DF, May02
;
;----------------------------------------------------------------------
