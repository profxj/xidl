;+ 
; NAME:
; hires_fittflat   
;     Version 1.1
;
; PURPOSE:
;  To create a 2D solution which describes the order curvature.  This
;  solution is derived from the individual traces created by
;  hires_trcflat and saved within the Trace structure.  The 2D fitting
;  algorithm is a simple least-squares algorithm.  The code then
;  attempts to extrapolate the solution for orders which are partially
;  on the CCD.  The code then makes a guess for the physical order
;  number which is not particularly accurate right now.  Finally, the
;  order structure (a key input for the HIRES pipeline) is written to
;  disk.
;  ---  The code primarily calls x_fittflat for everything
;
; CALLING SEQUENCE:
;   
;  hires_fittflat, hires, setup, [chip], /DEBUG, INNY=, INNT=, LHEDGE=,
;                  /CLOBBER
;
; INPUTS:
;   hires     -  HIRES structure
;   setup    -  Setup identifier 
;   [chip] -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;            (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;  One normalized flat per slit width
;
; OPTIONAL KEYWORDS:
;   INNY   -- Number of coefficients for fitting in vertical
;             direction.  (Default: 7)
;   INNT   -- Number of coefficients for fitting in horizontal
;             direction.  (Default: 6)
;   /DEBUG -- Turn debugging on
;   LHEDG  -- Used primarily for debugging
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   hires_fittflat, hires, 1, 1
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   28-Feb-2003 Written by SB
;   18-Apr-2003 Revised by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro hires_fittflat, hires, setup, chip, DEBUG=debug, INNY=inny, INNT=innt, $
                   LHEDGE=lhedg, CLOBBER=clobber

;
  if  N_params() LT 2  then begin 
      print,'Syntax - ' + $
        'hires_fittflat, hires, setup, [chip], /CLOBBER, /DEBUG, INNY=, ' + $
        'INNT=  [v1.1]'
      return
  endif

; Optional keywords
   if not keyword_set( CHIP ) then chip = [1L,2L,3L]

;  Loop on chip
   for ii=0L,n_elements(chip)-1 do begin
       qq = chip[ii]
       ;; CHIP
       case qq of
           1: begin
               print, 'hires_fittflat: Tracing BLUE trace flat'
               if NOT keyword_set(INNY) then nycoeff = 7 else nycoeff=inny[ii]
               if NOT keyword_set(INNT) then ntcoeff = 6 else ntcoeff=innt[ii]
           end
           2: begin
               print, 'hires_fittflat: Tracing GREEN trace flat'
               if NOT keyword_set(INNY) then nycoeff = 7 else nycoeff=inny[ii]
               if NOT keyword_set(INNT) then ntcoeff = 6 else ntcoeff=innt[ii]
           end
           3: begin
               print, 'hires_fittflat: Tracing RED trace flat'
               if NOT keyword_set(INNY) then nycoeff = 7 else nycoeff=inny[ii]
               if NOT keyword_set(INNT) then ntcoeff = 6 else ntcoeff=innt[ii]
           end
           else:stop
       endcase
   
       ;; Order structure
       ordr_fil = hires_getfil('ordr_str', setup, CHIP=qq, CHKFIL=chkf, /name)
       if CHKF NE 0 AND $
         not keyword_set( CLOBBER ) then begin
           print, 'hires_trcflat: Order structure exists. Moving on!', ordr_fil
           continue
       endif

       ;; Trace
       trc_str = hires_getfil('tflat_str', setup, CHIP=qq, FIL_NM=trc_fil, $
                             CHKFIL=chkf)
       if CHKF EQ 0 then begin
           print, 'hires_trcflat: Trace doesnt exist. ' + $
             'Run hires_trcflat first!', trc_fil
           continue
       endif

       ;; QA
       qafil = hires_getfil('qa_fittflat', setup, CHIP=qq)
       idx = where(hires.setup EQ setup AND hires.chip EQ qq $
                   AND hires.flg_anly NE 0)
       sz = lonarr(2)
       sz[0] = round(4096. / hires[idx[0]].rowbin) + 10L
       sz[1] = round(2048. / hires[idx[0]].colbin) + 10L

       ;; FIT
       x_fittflat, trc_str, ordr_fil, QAFIL=qafil, SZ=sz, $
         NYCOEFF=nycoeff, NTCOEFF=ntcoeff
       
   endfor

   ;; All done
   print, 'hires_fittflat: All done'
   return

end

