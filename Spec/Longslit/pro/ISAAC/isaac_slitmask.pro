;+
; NAME:
;   isaac_slitmask
;
; PURPOSE:
;   Make a slitmask for an isaac frame. 
;
; CALLING SEQUENCE:
;   isaac_slitmask, slitfile, tset_slits = tset_slits, slitmask = slitmask $
;                , outfile=outfile darkfile = darkfile
;
; INPUTS:
;   filename   - Image for finding the slits, which would typically
;                be a flat-field image, an arc image, or a sum of those
;   outfile    - Output file name with slit mask positions
;
; OPTIONAL INPUTS:
;   darkfile   - subtract dark to remove bias
;
;
; OUTPUTS:
;
; OPTIONAL OUTPUTS:
;   tset_slits - 2-element array of trace sets, where the first defines
;                the starting slit positions, and the second one defines
;                the ending slit positions
;   slitmask   - slitmask image
;   outfile    - output file    
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   ;
; REVISION HISTORY:
;   
;-
;------------------------------------------------------------------------------
PRO ISAAC_SLITMASK, filename, outfile = outfile, darkfile = darkfile $
                    , tset_slits = tset_slits $
                    , slitmask = slitmask


  HARDWIRE = 1
IF NOT KEYWORD_SET(HARDWIRE) THEN BEGIN
   long_slitmask, filename, minslit = 100, peakthresh = 0.7 $
                  , tset_slits = tset_slits1 $
                  , biasfile = darkfile
   
   tset_proto = $
      { func    :   tset_slits1[0].FUNC, $
        xmin    :   tset_slits1[0].XMIN, $
        xmax    :   tset_slits1[0].XMAX, $
        coeff   :    dblarr(3, 2), $
        dims    :   tset_slits1[0].DIMS $
      }
   
   
   tset_slits = replicate(tset_proto, 2)
;; left edge of slit0, right edge of slit1 are the same
   tset_slits[0].coeff[*, 0] = tset_slits1[0].coeff
   tset_slits[1].coeff[*, 1] = tset_slits1[1].coeff
;; right edge of slit0, left edge of slit1 are fixed
   oscan_rig1 = 510
   oscan_lef2 = 513
   tset_slits[0].coeff[0, 1] = oscan_rig1
   tset_slits[1].coeff[0, 0] = oscan_lef2
ENDIF ELSE BEGIN
   nx = 1024
   ny = 1024
   tset_slits = isaac_slitset(nx, ny)
ENDELSE
;; Construct this image to write to HDU #0 of the output file
   slitmask = long_slits2mask(tset_slits, nslit = nslit)
   ;; Write output file
   IF KEYWORD_SET(OUTFILE) THEN BEGIN
       splog, 'Writing output file'
       mwrfits, slitmask, outfile, hdr[*, 0], /create
       mwrfits, tset_slits, outfile
   ENDIF
   splog, 'Number of slits = ', nslit
;;   splog, 'Elapsed time = ', systime(1)-t0, ' sec'
   return
end


