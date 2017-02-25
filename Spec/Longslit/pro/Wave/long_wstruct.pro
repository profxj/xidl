;+
; NAME:
;   long_wstruct
;
; PURPOSE:
;   Create a structure which guides the data reduction process.  Parts
;   of the structure are initialized according to the instrument and
;   its configuration.
;
; CALLING SEQUENCE:
;
; INPUTS:
;   hdr  -- Image header
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;  Returns a structure describing the instrument configuration which
;  is used to guide the reduction steps.
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
;   10-Mar-2005  Written by S. Burles (MIT), David Schlegel (LBL), and
;                Joe Hennawi (UC Berkeley)
;   06-Jun-2014  Cleaned up extra '/' in calib_path and line_path (KG
;                Lee, MPIA) 
;-
;------------------------------------------------------------------------------
FUNCTION LONG_WSTRUCT, hdr, LINELIST=linelist, REID_FILE=REID_FILE $
                       , BIN_RATIO = bin_ratio

;;-----------
;;  Line list directory
calib_path = GETENV('XIDL_DIR') + '/Spec/Longslit/calib/linelists'
line_path = GETENV('XIDL_DIR') + '/Spec/Arcs/Lists'
;;----------
;; Initial guess for central wavelength and dispersion
telescope = strtrim(sxpar(hdr, 'TELESCOP'))
instrument = strtrim(sxpar(hdr, 'INSTRUME'))
detector = strtrim(sxpar(hdr, 'DETECTOR'))
binning = long(strsplit(sxpar(hdr, 'BINNING'), ',', /extract))
mask    = strtrim(sxpar(hdr, 'SLITNAME'))
mswave = double(strtrim(sxpar(hdr, 'MSWAVE')))
telid     = strtrim(sxpar(hdr, 'TELID'), 2)


wstruct = create_struct('INSTRUMENT', instrument $
                        , 'REID', 0L             $
                        , 'AUTOID', 0L           $
                        , 'CRUDEID', 0L            $
                        , 'PKWDTH', 0.0          $
                        , 'TOLER', 0.0           $
                        , 'FWEIGHT', 1L           $
                        , 'ICLSE', 4L           $
                        , 'THIN', 1L             $
                        , 'FORDR', 9L            $
                        , 'WAVE_CEN', 0.0        $
                        , 'BAND', ' '            $
                        , 'MODE', ' '            $
                        , 'NFIND', 0L            $
                        , 'DISP_GUESS', 0.0D     $
                        , 'FUNC', ' '            $
                        , 'NORD_CRUDE', 0L       $
                        , 'SIGREJ_CRUDE', 0.0    $
                        , 'PSIG_CRUDE', 0.0      $
                        , 'PSIG', fltarr(6)      $
                        , 'MXOFF', fltarr(6)     $
                        , 'SIGREJ', fltarr(6)    $
                        , 'NORD', fltarr(6)      $
                        , 'FLG_QUAL', lonarr(6)      $
                        , 'PSIG_REID', fltarr(6)     $
                        , 'MXOFF_REID', fltarr(6)    $
                        , 'SIGREJ_REID', fltarr(6)   $
                        , 'NORD_REID', lonarr(6)     $
                        , 'FLG_QUAL_REID', lonarr(6) $
                        , 'LINELIST_AP', ' '     $
                        , 'LINELIST',    ' '     $
                        , 'NPANIC', 0L           $
                        , 'NORD_PANIC', 0L       $
                        , 'DR_WAVE', [0.0, 0.0]  $
                        , 'MXSHIFT', 0.0  $
                        , 'REID_FILE', ' ' $
                        , 'BIN_RATIO',1. $  
                        , 'RADIUS', 0.0 $
                        , 'SIG_WPIX', 0.0 $
                        , 'PIX_MSK', [-1L,-1L])  ;; Implemented by KHRR to allow by-hand masking of edge pixels in the fit


wstruct.radius = 3.0
wstruct.SIG_WPIX = 3.0

if (instrument EQ 'LRISBLUE') then begin
    mask = strtrim(sxpar(hdr, 'SLITNAME'))
    grism = strtrim(sxpar(hdr, 'GRISNAME'))
    CASE grism OF
        '1200/3400': begin
            ;; parameters for arc peak finding
            wstruct.pkwdth        = 12.0/binning[1]
            wstruct.TOLER         = 4.0D/binning[1]
            wstruct.psig[0:4]     = [25.0, 18.0, 10.0, 7.0, 4.0]
            wstruct.MXOFF[0:4]    = [8.0, 8.0, 6.0, 6.0, 4.0]/binning[1]
            ;; paramters for wavelength soln fits
            wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:4]     = [4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 1L, 2L]
            wstruct.LINELIST      = line_path+'/lris_blue_1200.lst'
            wstruct.npanic        = 5L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1 ; reidentify using archived solutions
            wstruct.REID_FILE     = calib_path + '/lris_blue_1200.sav'
            wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
        END
        '600/4000': begin
            ;; parameters for arc peak finding
            wstruct.pkwdth        = 12.0D/binning[1]
            wstruct.TOLER         = 4.0D/binning[1]
            wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            ;; parameters for wavelength soln fits
            ;wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.0]
            wstruct.nord[0:4]     = [4L, 5L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_blue_600.lst'
            wstruct.npanic        = 5L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1 ; reidentify using archived solutions
            wstruct.REID_FILE     = calib_path + '/lris_blue_600.sav'
            wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
        END
        '400/3400': begin
           ;; Query on dichroic
            wstruct.pkwdth      = 12.00/binning[1]
            wstruct.TOLER       = 4.0D/binning[1]
            wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]  = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            wstruct.sigrej[0:4] = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:4]        = [4L, 5L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 5L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            case strtrim(sxpar(hdr,'DICHNAME'),2) of
               '560': begin
                  wstruct.LINELIST    = line_path+'/lris_blue_400_d560.lst'
                  wstruct.REID_FILE   = calib_path + '/lris_blue_400_d560.sav'
                  wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
               end
               '500': begin  ;; Not ideal, but better than the other
                  wstruct.LINELIST    = line_path+'/lris_blue_400_d560.lst'
                  wstruct.REID_FILE   = calib_path + '/lris_blue_400_d560.sav'
                  wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
               end
               else: begin
                  wstruct.LINELIST    = line_path+'/lris_blue_400.lst'
                  wstruct.REID_FILE   = calib_path + '/lris_blue_400.sav'
                  wstruct.BIN_RATIO   = binning[1]/1L
                  ;; Archived solutions are binned 1 in spectral
               end
            endcase
        END
        '300/5000': begin
            wstruct.pkwdth      = 12.0/binning[1] ; maxwidth of arc lines
            wstruct.TOLER       = 4.0D/binning[1]
            wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 2.5, 2.0]
            wstruct.nord[0:5]        = [4L, 4L, 5L, 5L, 6L, 6L]
            wstruct.FLG_QUAL[0:5]    = [2L, 2L, 2L, 3L, 3L, 3L]
            wstruct.LINELIST    = line_path+'/lris_blue_300.lst'
            wstruct.npanic      = 5L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + '/lris_blue_300.sav'
            wstruct.BIN_RATIO   = binning[1]/1L
        END
        ELSE: message, 'ERROR: Unknown grism'
    ENDCASE
ENDIF ELSE IF (instrument EQ 'LRIS') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'GRANAME'))
    CASE grating OF
        '150/7500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.nfind        = 50L
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_300.lst'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
               wstruct.disp_guess   = 2.45d0*binning[1]
               wstruct.dr_wave     = [0.9, 1.1]*(2.45d0)*binning[1]
               wstruct.AUTOID      = 1
            endif else begin
               ;; Best for newest detector (December 2010)
               wstruct.REID_FILE   = calib_path + '/lris_red_150_7500_bin1.sav'
               wstruct.dr_wave     = [0.9, 1.1]*(1.59d0)*binning[1]
               wstruct.AUTOID      = 0
               wstruct.REID        = 1
            endelse
        END
        '300/5000': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.nfind        = 50L
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_300.lst'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
               wstruct.disp_guess   = 2.45d0*binning[1]
               wstruct.dr_wave     = [0.9, 1.1]*(2.45d0)*binning[1]
               wstruct.AUTOID      = 1
            endif else begin
               ;; Best for newest detector (December 2010)
               wstruct.REID_FILE   = calib_path + '/lris_red_300_5000_d560_Ar.sav'
               wstruct.dr_wave     = [0.9, 1.1]*(1.59d0)*binning[1]
               wstruct.AUTOID      = 0
               wstruct.REID        = 1
            endelse
        END
        '400/8500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 2.0/binning[1]
            wstruct.psig[0:5]  = [7.0, 7.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5] = [3.0, 3.0, 3.0, 3.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 2L, 2L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_red_300.lst'
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
               wstruct.REID        = 1
               wstruct.REID_FILE   = calib_path + $
                                     '/lris_red_400_grangle_24.32_multi.sav'
;;           These are AUTOID paramaters that should be used for a
;;           a different grangle.
               wstruct.nfind        = 30L
               wstruct.disp_guess   = 1.86d0*binning[1]
               wstruct.FUNC  = 'CHEBY'
               wstruct.dr_wave       = [0.9, 1.1]*(1.86d0)*binning[1]
               wstruct.AUTOID      = 0
            endif else begin
               ;; Best for newest detector (December 2010)
               ;;wstruct.REID_FILE   = calib_path +
               ;;'/lris_400_8500.sav'
               ;; JFH 12-29-2013, new archive
               wstruct.REID        = 1
               wstruct.REID_FILE   = calib_path + '/lris_red_400_d560.sav'
               wstruct.disp_guess   = 1.19d*binning[1]
               wstruct.FUNC  = 'CHEBY'
               wstruct.dr_wave     = [0.9, 1.1]*(1.19d)*binning[1]
               wstruct.AUTOID      = 0
               wstruct.bin_ratio = binning[1] 
            endelse
        END
        '600/5000': BEGIN
           wstruct.pkwdth       = 6.0/binning[1]
           wstruct.TOLER        = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [4L, 4L, 5L, 5L, 6L, 6L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            ;; AUTOID paramaters
            wstruct.nfind        = 30L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
               ;; Currently using autoid but we should archive multislit soln's
               wstruct.disp_guess   = 1.28d0*binning[1]
               wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
               wstruct.AUTOID      = 1
               wstruct.REID        = 0
            endif else begin
               wstruct.REID_FILE   =  calib_path + '/lris_red_600_5000.sav'
               ;;calib_path + '/lris_newred_600_5000.sav'
               ;; TESTING
               wstruct.disp_guess   = 0.82d0*binning[1]
               wstruct.dr_wave     = [0.9, 1.1]*(0.82d0)*binning[1]
               wstruct.AUTOID      = 0
               wstruct.REID        = 1
               wstruct.BIN_RATIO   = float(binning[1])/2L
            endelse
            wstruct.FUNC  = 'CHEBY'
        END
        '600/7500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
            ;wstruct.sigrej[0:5]   = [4.0, 4.0, 4.0, 4.0, 4.0, 4.0]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:5]     = [3L, 4L, 4L, 5L, 6L, 6L]
            ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            case strtrim(sxpar(hdr,'DICHNAME'),2) of
               '560': begin
                  wstruct.LINELIST        = line_path+'/lris_red_600_d560.lst'
               end
               else: begin
                  wstruct.LINELIST        = line_path+'/lris_red_600.lst'
               end
            endcase
            ;; Kludge for case of red tilt observed with 7500. JFH 10.13.2016
            if float(sxpar(hdr, 'MSWAVE')) GE 8660.0 THEN  $
               wstruct.LINELIST = line_path+'/lris_red_600.lst'
            
            wstruct.npanic      = 20L
            wstruct.nord_panic  = 2L
            ;; AUTOID paramaters
            wstruct.nfind        = 30L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
               ;; Currently using autoid but we should archive multislit soln's
                wstruct.disp_guess   = 1.28d0*binning[1]
                wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
                ;wstruct.AUTOID      = 1
                ;wstruct.REID        = 0
                wstruct.AUTOID      = 0
                wstruct.REID        = 1
                wstruct.REID_FILE   = calib_path + '/lris_oldred_600_7500.sav'
                wstruct.BIN_RATIO   = binning[1]
             endif else begin
               wstruct.REID_FILE   = calib_path + '/lris_red_600_7500.sav'
               wstruct.disp_guess   = 0.82d0*binning[1]
               wstruct.dr_wave     = [0.9, 1.1]*(0.82d0)*binning[1]
               wstruct.AUTOID      = 0
               wstruct.REID        = 1
               wstruct.BIN_RATIO   = float(binning[1])/2L
               ;;arxiv soln binned spectrally
            endelse
            wstruct.FUNC  = 'CHEBY'
        END
        '600/10000': BEGIN
            wstruct.pkwdth       = 7.5/binning[1]
            wstruct.TOLER        = 2.0/binning[1]
            wstruct.FUNC  = 'CHEBY'
;            wstruct.psig[0] = 5.0
;            wstruct.MXOFF[0] = 3.0/binning[1]
;            wstruct.sigrej[0] = 3.0
;            wstruct.nord[0] = 5
;            wstruct.FLG_QUAL[0] = 2L
             wstruct.psig[0:5]     = [10.0, 10.0, 7.0, 7.0, 5.0, 5.0]
             wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 2.0, 2.0]/binning[1]
             wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
             wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
             wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
                wstruct.disp_guess   = 1.28d0*binning[1]
                wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
                wstruct.AUTOID      = 0
                wstruct.REID        = 1
                wstruct.REID_FILE   = calib_path + '/lris_red_600-10000_d680_mswave_8051.sav' ;; old chip
            endif else begin
                wstruct.disp_guess   = 0.82d0*binning[1]
                wstruct.dr_wave     = [0.9, 1.1]*(0.82d0)*binning[1]
                wstruct.AUTOID      = 0
                wstruct.REID        = 1
                wstruct.REID_FILE   = calib_path + '/lris_red_600_10000.sav' ;; New chip
                wstruct.BIN_RATIO   = binning[1]/1L

            endelse
            wstruct.nfind        = 30L
        END
        '831/8200': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [4.0, 4.0, 4.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:5]     = [3L, 3L, 3L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.npanic      = 5L ; Should be 10L
            wstruct.nord_panic  = 2L
            ;; Currently using autoid but we should archive multislit soln's
            if strmid(sxpar(hdr,'DATE'),10) LT '2013-1-1' then begin
                wstruct.REID_FILE   = calib_path + $
                  '/lris_red_831_8200.sav'
            endif else begin
                ;; new grating
                wstruct.REID_FILE   = calib_path + $
                  '/lris_red_831_8200_gold_d560.sav'
            endelse
            ;; AUTOID paramaters
            wstruct.nfind        = 30L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
                wstruct.disp_guess   = 1.28d0*binning[1]
                stop
                wstruct.dr_wave     = [0.9, 1.1]*(1.28d0)*binning[1]
                wstruct.AUTOID      = 1
                wstruct.REID        = 0
            endif else begin
                wstruct.disp_guess   = 0.58*binning[1]
                wstruct.dr_wave     = [0.9, 1.1]*(0.58d0)*binning[1]
                wstruct.AUTOID      = 0
                wstruct.REID        = 1
            endelse
            wstruct.FUNC  = 'CHEBY'
        END
        '900/5500': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.nfind        = 30L
            wstruct.disp_guess   = 0.53d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
            wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_600.lst'
            wstruct.REID        = 1
            ;wstruct.REID_FILE   = calib_path + '/lris_red_900_5500.sav'
            wstruct.REID_FILE   = calib_path + '/lris_red_900_5500_d460.sav'
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.dr_wave     = [0.9, 1.1]*(0.53d0)*binning[1]
            wstruct.AUTOID      = 0
        END
        '1200/7500': BEGIN
            wstruct.pkwdth        = 8.0/binning[1] ; width of peaks to center on
            wstruct.TOLER         = 3.0/binning[1] ; minimum centroiding diff
            wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 4.0, 3.0]
            wstruct.MXOFF[0:5]    = [4.0, 4.0, 4.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [2L, 2L, 2L, 2L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_red_1200.lst'
            wstruct.npanic      = 4L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            ;;MF Jan 2011. Update calib for redside upgrade
            if strmid(sxpar(hdr,'DATE'),10) LT '2010-12-10' then begin
                wstruct.REID_FILE   = calib_path + $
                  '/lris_red_1200-7500_d460_mswave_5137.sav'
            endif else begin
                ;;this is a very red setup
                wstruct.REID_FILE   = calib_path + $
                  '/lris_red_1200-7500_d560_mswave_9200.sav'
                wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
            endelse

            ;IF strmatch(mask, '*long*') AND mswave GT 6000.0D THEN BEGIN
;            wstruct.nfind       = 50L ; nlines for arc_pairs
;            wstruct.AUTOID      = 1L
;            wstruct.disp_guess   = 0.636d0*binning[1] ; dispersion
;            wstruct.dr_wave      = [0.95, 1.05]*(0.636d0)*binning[1]
;            wstruct.FUNC  = 'CHEBY' ; function to use for firts
            ;;ENDIF ELSE BEGIN
            ;;ENDELSE
        END
        '831/8200': BEGIN
            wstruct.pkwdth       = 6.0/binning[1]
            wstruct.TOLER        = 1.5/binning[1]
            wstruct.FUNC  = 'CHEBY'
;            wstruct.psig[0] = 5.0
;            wstruct.MXOFF[0] = 3.0/binning[1]
;            wstruct.sigrej[0] = 3.0
;            wstruct.nord[0] = 5
;            wstruct.FLG_QUAL[0] = 2L
             wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
             wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 2.0, 2.0]/binning[1]
             wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
             wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
             wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/lris_red_821.lst'
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            if strmid(sxpar(hdr,'DATE'),10) LT '2009-07-01' then begin
                wstruct.disp_guess   = 0.928d0*binning[1]
                wstruct.dr_wave     = [0.9, 1.1]*(0.928d0)*binning[1]
                wstruct.AUTOID      = 1
                wstruct.REID        = 0
            endif else begin
                wstruct.disp_guess   = 0.59d0*binning[1]
                wstruct.dr_wave     = [0.9, 1.1]*(0.59d0)*binning[1]
                wstruct.AUTOID      = 0
                wstruct.REID        = 1
            endelse
            wstruct.REID_FILE   = calib_path + '/lris_red_821_8200.sav'  ;; New chip
            wstruct.nfind        = 30L
         END
        '1200/9000': BEGIN
            ;;BPH Feb 2013. New red grating
            wstruct.pkwdth        = 8.0/binning[1] ; width of peaks to center on
            wstruct.TOLER         = 3.0/binning[1] ; minimum centroiding diff
            wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 4.0, 3.0]
            wstruct.MXOFF[0:5]    = [4.0, 4.0, 4.0, 4.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [2L, 2L, 2L, 2L, 2L, 2L]
            wstruct.LINELIST      = line_path+'/lris_red_1200.lst'
            wstruct.npanic      = 4L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1

            wstruct.REID_FILE   = calib_path + $
                  '/lris_red_1200_900_d560.sav'

          END

        ELSE: message, 'ERROR: Unknown grating'
     ENDCASE
ENDIF ELSE IF strmatch(instrument,'FORS2*') THEN BEGIN
   ind_det =  WHERE(stregex(hdr, 'HIERARCH ESO DET CHIP1 NAME', /bool))
   detector = strcompress(strmid(hdr[ind_det], 30, 14))
   grism = esopar(hdr, 'HIERARCH ESO INS GRIS1 NAME')
   binx = esopar(hdr, 'HIERARCH ESO DET WIN1 BINX')
   biny = esopar(hdr, 'HIERARCH ESO DET WIN1 BINY')
   binning = [biny, binx] ;; images are transposed
   ;stop
   CASE grism OF
      'GRIS_300V': BEGIN
          ;; Changes made here by KHRR, 6.9.2011
         wstruct.pkwdth       = 6.0/(binning[1]/2.0d)
         wstruct.TOLER        = 2.0/(binning[1]/2.0d)
         wstruct.nfind        = 40L
         wstruct.FUNC  = 'CHEBY'
         wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
         wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/(binning[1]/2.d)
         wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
         ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
         wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 5L]
         wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
         ;wstruct.LINELIST =  '/Users/joe/FORS_testbed/fors2_GRIS_300V+10_1.lst'
         wstruct.LINELIST =  line_path + '/fors_300V+10.lst'
         wstruct.npanic      = 8L
         wstruct.nord_panic  = 2L
         wstruct.disp_guess  = 1.68d*(binning[1]/2.0d)
         wstruct.dr_wave     = [0.9, 1.1]*(1.68d)*(binning[1]/2.0d)
         wstruct.AUTOID      = 0
         wstruct.REID        = 1
         wstruct.REID_FILE   = calib_path + '/fors_300V+10_mswave.sav'
      END
      'GRIS_600B': BEGIN
         wstruct.pkwdth       = 6.0/(binning[1]/2.0d)
         wstruct.TOLER        = 2.0/(binning[1]/2.0d)
         wstruct.nfind        = 40L
         wstruct.FUNC  = 'CHEBY'
         wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
         wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/(binning[1]/2.d)
         wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
         ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
         wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 5L]
         wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
         wstruct.LINELIST =  line_path + '/fors_600B+22.lst'
         wstruct.npanic      = 8L
         wstruct.nord_panic  = 2L
         wstruct.disp_guess  = 0.75d*(binning[1]/2.0d)
         wstruct.dr_wave     = [0.9, 1.1]*(0.75d)*(binning[1]/2.0d)
         wstruct.AUTOID      = 0
         wstruct.REID        = 1
         wstruct.REID_FILE   = calib_path + '/fors_600B+22_mswave.sav'
      END
      'GRIS_1200B': BEGIN        
          ;; Changes made here by KHRR, 24.3.2012
         wstruct.pkwdth       = 6.0/(binning[1]/2.0d)
         wstruct.TOLER        = 2.0/(binning[1]/2.0d)
         wstruct.nfind        = 40L
         wstruct.FUNC  = 'CHEBY'
         wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
         wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/(binning[1]/2.d)
         wstruct.sigrej[0:5]   = [3.0, 3.0, 2.5, 2.5, 2.0, 2.0]
         wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
         ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
         wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
         ;wstruct.LINELIST =  line_path + '/fors_all.lst'
         wstruct.LINELIST =  line_path + '/fors_1200B.lst'
         wstruct.npanic      = 8L
         wstruct.nord_panic  = 2L
         wstruct.disp_guess  = 1.68d*(binning[1]/2.0d)
         wstruct.dr_wave     = [0.9, 1.1]*(1.68d)*(binning[1]/2.0d)
         wstruct.PIX_MSK      = [60L,0L]
         wstruct.AUTOID      = 0
         wstruct.REID        = 1
         wstruct.REID_FILE   = calib_path + '/fors_1200B+97_mswave_2015jul08.sav'
     END
     'GRIS_600V': BEGIN
         wstruct.pkwdth       = 6.0/(binning[1]/2.0d)
         wstruct.TOLER        = 2.0/(binning[1]/2.0d)
         wstruct.nfind        = 40L
         wstruct.FUNC  = 'CHEBY'
         wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
         wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/(binning[1]/2.d)
         wstruct.sigrej[0:5]   = [3.0, 3.0, 2.5, 2.5, 2.0, 2.0]
         ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
         wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 5L]
         wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
         wstruct.LINELIST =  line_path + '/fors_600V.lst'
         wstruct.npanic      = 8L
         wstruct.nord_panic  = 2L
         wstruct.disp_guess  = 0.75d*(binning[1]/2.0d)
         wstruct.dr_wave     = [0.9, 1.1]*(0.75d)*(binning[1]/2.0d)
         wstruct.AUTOID      = 0
         wstruct.REID        = 1
         ;wstruct.REID_FILE   = calib_path + '/fors_600B+22_mswave.sav'
         wstruct.REID_FILE   = '/Users/rubin/Research/FORS2/2012nov14/fors_600V+94_mswave.sav'
      END
     'GRIS_600RI': BEGIN
         wstruct.pkwdth       = 6.0/(binning[1]/2.0d)
         wstruct.TOLER        = 2.0/(binning[1]/2.0d)
         wstruct.nfind        = 40L
         wstruct.FUNC  = 'CHEBY'
         wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
         wstruct.MXOFF[0:5]    = [6.0, 6.0, 6.0, 5.0, 3.0, 2.0]/(binning[1]/2.d)
         wstruct.sigrej[0:5]   = [3.0, 3.0, 2.5, 2.5, 2.0, 2.0]
         ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
         wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 5L]
         wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
         wstruct.LINELIST =  line_path + '/fors_all.lst'
         wstruct.npanic      = 8L
         wstruct.nord_panic  = 2L
         wstruct.disp_guess  = 0.75d*(binning[1]/2.0d)
         wstruct.dr_wave     = [0.9, 1.1]*(0.75d)*(binning[1]/2.0d)
         wstruct.PIX_MSK      = [50L,0L]
         wstruct.AUTOID      = 0
         wstruct.REID        = 1        
         wstruct.REID_FILE   = calib_path + '/fors_600RI+19_mswave.sav'
      END
      ELSE: message, 'ERROR: Unknown grism'
   ENDCASE
ENDIF ELSE IF strmatch(telescope,'LBT-SX') THEN BEGIN
    binning = [1, 1]
    ;grating = strtrim(sxpar(hdr, 'GRATNAME'),2)
    grating = strtrim(sxpar(hdr, 'GRATINFO'),2)
    grating = strmid(grating,0,7)
    mask = strtrim(sxpar(hdr, 'MASKNAME'),2)
    ;; Blue side
    CASE grating OF
       '400l/mm': BEGIN
        ;'G450L': BEGIN
           ;; parameters for arc peak finding
            wstruct.pkwdth        = 12.0D/binning[1]
            wstruct.TOLER         = 4.0D/binning[1]
            wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
            wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
            ;; parameters for wavelength soln fits
            ;wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:4]     = [4L, 5L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
            ;wstruct.LINELIST      = line_path+'/lris_blue_600.lst'
            wstruct.LINELIST      = line_path+'/mods_blue_400.lst'
            wstruct.npanic        = 5L
            wstruct.nord_panic    = 2L
            wstruct.AUTOID        = 0L
            wstruct.REID          = 1 ; reidentify using archived solutions
            ;wstruct.REID_FILE     = '/home/rubin/data/MODS/Jan12/mods_blue_400ms.sav'
            wstruct.REID_FILE     = calib_path+'/mods_blue_400ms.sav'
            wstruct.BIN_RATIO     = binning[1] ;arxiv soln unbinned spectrally
         END
       '250l/mm': BEGIN
        ;'G670L': BEGIN
           wstruct.pkwdth       = 6.0/binning[1]
           wstruct.TOLER        = 1.5/binning[1]
           wstruct.psig[0:5]     = [12.0, 9.0, 7.0, 5.0, 4.0, 4.0]
           wstruct.MXOFF[0:5]    = [5.0, 5.0, 5.0, 4.0, 3.0, 3.0]/binning[1]
           wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
            ;wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
           wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
           wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
           ;wstruct.LINELIST        = line_path+'/lris_red_600_d560.lst'
           ;wstruct.LINELIST      =  '/home/rubin/data/MODS/Jan12/mods_red_670.lst'
           wstruct.LINELIST      = line_path+'/mods_red_670.lst'
           wstruct.npanic      = 20L
           wstruct.nord_panic  = 2L
            ;; AUTOID paramaters
           wstruct.nfind        = 30L
           wstruct.disp_guess   = 0.847d0*binning[1]
           wstruct.dr_wave     = [0.9, 1.1]*(0.847d0)*binning[1]
           wstruct.AUTOID      = 0L
           wstruct.REID        = 1
           ;wstruct.REID_FILE   = calib_path + '/lris_oldred_600_7500.sav'
                                ;wstruct.REID_FILE   =
                                ;'/home/rubin/data/MODS/Jan12/mods_red_670.sav'
           if STRMATCH(mask,'*LS*') then begin
               wstruct.REID_FILE   = calib_path+'/mods_red_670.sav'
           endif else begin
               wstruct.REID_FILE   = calib_path+'/mods_red_670ms.sav'
           endelse
           wstruct.BIN_RATIO   = binning[1]
           ;;arxiv soln unbinned
           wstruct.FUNC  = 'CHEBY'
        END
        ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
 ENDIF ELSE IF strmatch(instrument, 'DEIMOS*') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'GRATENAM'), 2)
   gpos = strtrim(sxpar(hdr, 'GRATEPOS'), 2)
   wave = sxpar(hdr, 'G'+strtrim(gpos, 2)+'TLTWAV')
   CASE grating OF
      '1200G': BEGIN
         wstruct.pkwdth       = 5.0/binning[1]
         wstruct.TOLER        = 2.0/binning[1]
         wstruct.nfind        = 50L
         wstruct.FUNC  = 'CHEBY'
         wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic      = 15L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.MXSHIFT     = 50
            if wave LT 7000 then begin
                wstruct.LINELIST  = line_path+'/deimos_blue.lst'
                wstruct.REID_FILE   = calib_path + '/deimos_1200_blue.sav'
            endif else stop
        END
        '600ZD': BEGIN
            wstruct.pkwdth       = 5.0/binning[1]
            wstruct.TOLER        = 2.0/binning[1]
            wstruct.nfind        = 50L
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 2.5, 2.5]
            wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic      = 15L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            wstruct.MXSHIFT     = 150
            wstruct.LINELIST  = line_path+'/deimos_blue.lst'
            wstruct.REID_FILE   = calib_path + '/deimos_600.sav' ;; This file needs tuning
        END
        ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
endif else if (stregex(instrument,'.*kast.*',/boolean,/fold_case) eq 1) OR $
  (stregex(sxpar(hdr,'VERSION'),'kast*', /bool, /fold_case) EQ 1) then begin
    binning = [1,1]
    ;; Blue side
    if strtrim(sxpar(hdr,'SPSIDE'),2) EQ 'blue' OR $
      strmid(sxpar(hdr,'VERSION'),0,5) EQ 'kastb' then begin
        grating = strtrim(sxpar(hdr,'GRISM_N'),2)
        CASE grating OF
            '830/3460': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 7L
                wstruct.nord_panic  = 2L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_blue3.lst'
                wstruct.REID_FILE   = calib_path + '/kast_830_3460.sav'
            END
            '600/4310': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [6.0, 6.0, 6.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 4L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 2L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 200
                wstruct.LINELIST  = line_path+'/kast_blue2.lst'
                wstruct.REID_FILE   = calib_path + '/kast_600_4310.sav'
            END
            '452/3306': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 2L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                mjd_newccd = 2454534L
                dateobs = sxpar(hdr, 'DATE-OBS')
                mjd = x_setjdate(strmid(dateobs,0,10))

                wstruct.LINELIST  = line_path+'/kast_blue3.lst'
                if mjd LT 2454534L then $
                  wstruct.REID_FILE   = calib_path + '/kast_old_452_3306.sav'$
                else wstruct.REID_FILE   = calib_path + '/kast_452_3306.sav'
            END
            ELSE: message, 'ERROR: Unknown grating'
        ENDCASE
    endif else begin
        ;; Kast
        grating = strtrim(sxpar(hdr,'GRATNG_N'),2)
        CASE grating OF
            '1200/5000': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.disp_guess   = 1.17d0*binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst'
                wstruct.REID_FILE   = calib_path + '/kast_1200_5000.sav'
            END
            '1200 5000': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.disp_guess   = 1.17d0*binning[1]
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst'
                wstruct.REID_FILE   = calib_path + '/kast_1200_5000.sav'
            END
            '600/5000': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 2.35d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 4L, 4L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 100
                wstruct.LINELIST  = line_path+'/kast_red.lst'
                ;wstruct.AUTOID      = 1
                wstruct.dr_wave     = [0.9, 1.1]*(2.35d0)*binning[1]
                wstruct.REID_FILE   = calib_path + '/kast_600_5000_d46.sav'
            END
            '600/7500': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 2.32d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 4L, 4L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 100
                wstruct.LINELIST  = line_path+'/kast_red.lst'
;                wstruct.AUTOID      = 1
                wstruct.dr_wave     = [0.9, 1.1]*(2.32d0)*binning[1]
                wstruct.REID_FILE   = calib_path + '/kast_600_7500_d55.sav'
            END
            '830/8460': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 1.70d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [10.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 2L, 3L, 3L, 3L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 100
                wstruct.LINELIST  = line_path+'/kast_red.lst'
;                wstruct.AUTOID      = 1
                wstruct.dr_wave     = [0.9, 1.1]*(1.70d0)*binning[1]
                wstruct.REID_FILE   = calib_path + '/kast_830_8460_d55.sav'
            END
            '300/7500': BEGIN
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 2.32d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 50
                wstruct.LINELIST  = line_path+'/kast_red.lst'
;                wstruct.AUTOID      = 1
                wstruct.dr_wave     = [0.9, 1.1]*(2.32d0)*binning[1]
                wstruct.REID_FILE   = calib_path + '/kast_300_7500.sav'
            END
            ELSE: message, 'ERROR: Unknown grating'
        ENDCASE
    endelse
ENDIF ELSE IF strmatch(instrument,'DIS') THEN BEGIN
    binning = [1, 1]
    ;; Blue side
    grating = strtrim(sxpar(hdr, 'GRATING'), 2)
    CASE detector OF
       'blue': BEGIN
          case grating of
             'B1200': begin
                wstruct.pkwdth        = 3.0/binning[1]
                wstruct.TOLER         = 2.0/binning[1]
                wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 5.0]
                wstruct.MXOFF[0:5]    = [10.0, 8.0, 6.0, 4.0, 3.0, 3.0]/binning[1]
                wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 2.5]
                wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
                wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic        = 10L
                wstruct.nord_panic    = 2L
                wstruct.REID          = 1
                wstruct.MXSHIFT       = 200
                wstruct.LINELIST      = line_path+'/dis_henear.lst'
                wstruct.REID_FILE     = calib_path + '/dis_B1200.sav'
             END
             ELSE : BEGIN
                wstruct.pkwdth        = 3.0/binning[1]
                wstruct.TOLER         = 2.0/binning[1]
                wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 5.0]
                wstruct.MXOFF[0:5]    = [10.0, 8.0, 6.0, 4.0, 3.0, 3.0]/binning[1]
                wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 2.5]
                wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
                wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic        = 10L
                wstruct.nord_panic    = 2L
                wstruct.REID          = 1
                wstruct.MXSHIFT       = 200
                wstruct.LINELIST      = line_path+'/dis_henear.lst'
                wstruct.REID_FILE     = calib_path + '/dis_B400.sav'
             END
          endcase
       END
       'red': BEGIN
          case grating of
             'R1200': begin
                wstruct.pkwdth        = 5.0/binning[1]
                wstruct.TOLER         = 2.0/binning[1]
                wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic        = 10L
                wstruct.nord_panic    = 2L
                wstruct.REID          = 1
                wstruct.MXSHIFT       = 200
                wstruct.LINELIST  = line_path+'/dis_henear.lst'
                wstruct.REID_FILE   = calib_path + '/dis_R1200.sav'
             end
             else : begin
                wstruct.pkwdth        = 5.0/binning[1]
                wstruct.TOLER         = 2.0/binning[1]
                wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
                wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic        = 10L
                wstruct.nord_panic    = 2L
                wstruct.REID          = 1
                wstruct.MXSHIFT       = 200
                wstruct.LINELIST  = line_path+'/dis_henear.lst'
                wstruct.REID_FILE   = calib_path + '/dis_R300.sav'
             end
          ENDCASE
       end
       ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
 ENDIF ELSE IF strmatch(instrument,'*ISIS*') THEN BEGIN
    bins = sxpar(hdr, 'CCDSUM')
    binning_str = strsplit(bins, /extract)
    bin1 = long(binning_str[0])
    bin2 = long(binning_str[1])
    binning = [bin1,bin2]       
    isiswguess = float(strtrim(sxpar(hdr, 'CENWAVE')))
    grating = strtrim(sxpar(hdr, 'ISIGRAT'))
    case grating of
       'R300B': begin
          wstruct.pkwdth        = 6.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 5.0]
          wstruct.MXOFF[0:5]    = [10.0, 8.0, 6.0, 4.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 2.5]
          wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
          wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/ISIS_CuNeCuAr.lst'
          wstruct.REID_FILE     = calib_path + '/ISIS_CuNeCuAr_R300B_4200.sav'
       end
       'R600R': begin
          wstruct.pkwdth        = 6.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 5.0]
          wstruct.MXOFF[0:5]    = [10.0, 8.0, 6.0, 4.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 2.5]
          wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
          wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/ISIS_CuNeCuAr.lst'
          if(isiswguess lt 7500.) then wstruct.REID_FILE=calib_path+'/ISIS_CuNeCuAr_R600R_6400.sav' $
          else wstruct.REID_FILE=calib_path+'/ISIS_CuNeCuAr_R600R_8200.sav'
       end
       'R316R': begin
          wstruct.pkwdth        = 6.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 5.0]
          wstruct.MXOFF[0:5]    = [10.0, 8.0, 6.0, 4.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 2.5]
          wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
          wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/ISIS_CuNeCuAr.lst'
          wstruct.REID_FILE     = calib_path + '/ISIS_CuNeCuAr_R316R_6400.sav'
       end
       else: message, 'ERROR: Unknown grating'
    endcase
 ENDIF ELSE IF strmatch(telescope, 'kp4m') THEN BEGIN
    binning = [1, 1]
    grating = strcompress(sxpar(hdr, 'DISPERSE'), /rem)
    if grating eq 0 then grating = 'KPC10A'
    CASE grating OF
       'KPC10A': BEGIN
          wstruct.pkwdth        = 3.0/binning[1]
          wstruct.TOLER         = 2.0/binning[1]
          wstruct.psig[0:1]     = [10.0, 10.0]
          wstruct.MXOFF[0:1]    = [2.0, 2.0]/binning[1]
          wstruct.sigrej[0:1]   = [2.0, 2.0]
          wstruct.nord[0:1]     = [5L, 5L]
          wstruct.FLG_QUAL[0:1] = [2L, 2L]
          wstruct.FUNC  = 'CHEBY'
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/dis_henear.lst'
          wstruct.REID_FILE     = calib_path + '/kp4m_KPC10A.sav'
        END
        ELSE: message, 'ERROR: Unknown grating'
     ENDCASE
 ENDIF ELSE IF strmatch(telid, '200') THEN BEGIN
    binning = [1, 1]
    fpa = strtrim(sxpar(hdr, 'FPA'), 2)
    grating = strtrim(sxpar(hdr, 'GRATING'))
    CASE fpa OF
       'DBSP_BLUE': BEGIN
          case grating of
             '1200/5000': begin
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 2.0/binning[1]
                wstruct.nfind        = 50L
                wstruct.disp_guess   = 0.21d0*binning[1]
                wstruct.FUNC  = 'CHEBY'
                wstruct.psig[0:5]   = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
                wstruct.MXOFF[0:5] = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
                wstruct.sigrej[0:5]      = [4.0, 3.0, 3.0, 3.0, 3.0, 2.5]
                wstruct.nord[0:5]        = [3L, 3L, 4L, 4L, 4L, 4L]
                wstruct.FLG_QUAL[0:5]    = [1L, 1L, 1L, 2L, 2L, 2L]
                wstruct.npanic      = 10L
                wstruct.nord_panic  = 3L
                wstruct.REID        = 1
                wstruct.MXSHIFT     = 100
;                wstruct.AUTOID      = 1
                wstruct.LINELIST      = line_path+'/fear.lst'
                wstruct.REID_FILE     = calib_path + '/p200_DBSP_blue_fear.sav'
             end
             '600/4000': begin
                wstruct.pkwdth       = 5.0/binning[1]
                wstruct.TOLER        = 3.0/binning[1]
                wstruct.pkwdth        = 5.0/binning[1]
                wstruct.TOLER         = 3.0/binning[1]
                wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
                wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
                wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]     = [5L, 5L, 5L, 5L, 5L, 5L]
                wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
                wstruct.npanic        = 10L
                wstruct.nord_panic    = 2L
                wstruct.REID          = 1
                wstruct.MXSHIFT       = 200
                wstruct.LINELIST      = line_path+'/p200_DBSP_FeAr.lst'  
                wstruct.REID_FILE     = calib_path + '/p200_DBSP_blue_FeAr.sav'
             end
             else: begin
                wstruct.pkwdth        = 5.0/binning[1]
                wstruct.TOLER         = 3.0/binning[1]
                wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
                wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
                wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
                wstruct.nord[0:5]     = [5L, 5L, 5L, 5L, 5L, 5L]
                wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
                wstruct.npanic        = 10L
                wstruct.nord_panic    = 2L
                wstruct.REID          = 1
                wstruct.MXSHIFT       = 200
                wstruct.LINELIST      = line_path+'/p200_DBSP_blue.lst'
                wstruct.REID_FILE     = calib_path + '/p200_DBSP_blue.sav'
             end
          endcase
        END
        'DBSP_RED': BEGIN
           wstruct.pkwdth       = 5.0/binning[1]
            wstruct.TOLER        = 3.0/binning[1]
            wstruct.nfind        = 40L
            wstruct.disp_guess   = 2.12d0*binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 5.0, 5.0]
            wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.LINELIST        = line_path+'/dis_henear.lst'
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.dr_wave     = [0.9, 1.1]*(2.46d0)*binning[1]
            wstruct.AUTOID      = 1
         END
        ;;MF 2013 Added this 
        'DBSP_RED2': BEGIN
           wstruct.pkwdth       = 5.0/binning[1]
           wstruct.TOLER        = 3.0/binning[1]
           wstruct.nfind        = 40L
           wstruct.disp_guess   = 2.12d0*binning[1]
           wstruct.FUNC  = 'CHEBY'
           wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 5.0, 5.0]
           wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
           wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
           wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
           wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
           wstruct.LINELIST        = line_path+'/dis_henear.lst'
           wstruct.REID          = 1
           wstruct.REID_FILE     = calib_path + '/p200_DBSP_red_NeHeAr.sav'
           wstruct.npanic      = 10L
           wstruct.nord_panic  = 2L
           wstruct.dr_wave     = [0.9, 1.1]*(2.46d0)*binning[1]
           ;wstruct.AUTOID      = 1
        END
        ELSE: message, 'ERROR: Unknown grating'
     ENDCASE
 ENDIF ELSE IF strmatch(telescope, '*CA-3.5*') THEN BEGIN
    binning = [1, 1]
    caha_hh = caha_hdr(hdr)
    path = strcompress(sxpar(caha_hh, 'PATH'), /rem)
    grating =  (path EQ 'BLUE' ? sxpar(caha_hh, 'GRAT1_NA') : '') $
               + (path EQ 'RED'  ? sxpar(caha_hh, 'GRAT2_NA') : '')
    slitwidth = strtrim(sxpar(caha_hh, 'SLIT_WID'))
    CASE path OF
       'BLUE': BEGIN
          wstruct.pkwdth        = 5.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
          wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
          wstruct.nord[0:5]     = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/CAHA_T13_BS6800.lst'
          wstruct.REID_FILE     = calib_path + '/CAHA_T13_BS6800.sav'
          wstruct.fweight       = 0
        END
       'RED': BEGIN
          wstruct.pkwdth       = 5.0/binning[1]
          wstruct.TOLER        = 3.0/binning[1]
          wstruct.nfind        = 40L
          wstruct.disp_guess   = 2.12d0*binning[1]
          wstruct.FUNC  = 'CHEBY'
          wstruct.psig[0:5]     = [10.0, 7.0, 7.0, 5.0, 5.0, 5.0]
          wstruct.MXOFF[0:5]    = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
          wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]
          wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
          wstruct.LINELIST        = line_path+'/CAHA_T07_BS6800.lst'
          wstruct.npanic      = 10L
          wstruct.nord_panic  = 2L
          wstruct.dr_wave     = [0.9, 1.1]*(2.46d0)*binning[1]
          wstruct.AUTOID      = 1
          wstruct.fweight     = 1
       END
       ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
 ENDIF ELSE IF strmatch(telescope, '*CA-2.2*') THEN BEGIN
    binning = [1, 1]
    caha_hh = hdr
    grating = strtrim(sxpar(caha_hh,'INSGRNAM'),2)
;    slitwidth = strtrim(sxpar(caha_hh, 'SLIT_WID'))
    CASE grating OF
       'green-200': BEGIN
          wstruct.pkwdth        = 5.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
          wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
          wstruct.nord[0:5]     = [3L, 4L, 4L, 4L, 4L, 4L]
          wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 50
          wstruct.LINELIST      = line_path+'/CAHA_CAFOS.lst'
          wstruct.REID_FILE     = calib_path + '/CAHA_2.2m_CAFOS_green200.sav'
          wstruct.fweight       = 0
        END
       ELSE: message, 'ERROR: Unknown grating'
   ENDCASE
ENDIF ELSE IF strmatch(telescope, '*SOAR*') THEN BEGIN
    bins = sxpar(hdr, 'CCDSUM')
    binning_str = strsplit(bins, /extract)
    bin1 = long(binning_str[0])
    bin2 = long(binning_str[1])
    binning = [bin1,bin2]       
    grating = strtrim(sxpar(hdr, 'GRATING'))
    CASE grating OF
       'RALC_300': BEGIN
          wstruct.pkwdth        = 5.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.FUNC  = 'CHEBY'
          wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
          wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 2.5, 2.5, 2.5]
          wstruct.nord[0:5]     = [4L, 4L, 4L, 4L, 4L, 4L]
          wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = '/export/qso3/database/data/SOAR/Redux_Kate/goodman_300_hgar.lst'
          ;wstruct.LINELIST      = line_path+'/CAHA_T13_BS6800.lst'
          ;wstruct.REID_FILE     = calib_path + '/CAHA_T13_BS6800.sav'
          wstruct.REID_FILE     = '/export/qso3/database/data/SOAR/Redux_Kate/2013jan11/SOAR_Goodman_300.sav'
          wstruct.fweight       = 0
          wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally
          ;stop
        END
       
       ELSE: message, 'ERROR: Unknown grating'
    ENDCASE
 ENDIF ELSE IF strmatch(instrument, 'EFOSC') THEN BEGIN
    ind_binx = WHERE(stregex(hdr, 'HIERARCH ESO DET WIN1 BINX', /bool))
    binx = double(strmid(hdr[ind_binx], 30, 14))
    ind_biny = WHERE(stregex(hdr, 'HIERARCH ESO DET WIN1 BINY', /bool))
    biny = double(strmid(hdr[ind_binx], 30, 14))
    binning = [binx, biny]
    ind_gris = $
       WHERE(stregex(hdr, 'HIERARCH ESO INS GRIS1 NAME' $
                     , /bool))
    grism = strcompress(repstr(strmid(hdr[ind_gris], 30, 14), "'", ''), /rem)
    grism = repstr(grism, '#', '')
    CASE grism OF
       'Gr11': BEGIN
          wstruct.pkwdth        = 6.0/binning[1]
          wstruct.TOLER         = 3.0/binning[1]
          wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
          wstruct.MXOFF[0:5]    = [4.0, 4.0, 3.0, 3.0, 3.0, 3.0]/binning[1]
          wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
          wstruct.nord[0:5]     = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.FLG_QUAL[0:5] = [3L, 3L, 3L, 3L, 3L, 3L]
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'/NTT_EFOSC_G11.lst'
          wstruct.REID_FILE     = calib_path + '/NTT_EFOSC_G11_KR.sav'
          ;wstruct.REID_FILE     = '/Users/rubin/Research/CGM_DLA/EFOSC/Night1/NTT_EFOSC_G11_calib.sav'
          wstruct.fweight       = 0
        END
       ELSE: message, 'ERROR: Unknown grism'
    ENDCASE
ENDIF ELSE IF strmatch(telescope, 'mmt') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'DISPERSE'))
    detector = sxpar(hdr, 'DETECTOR')
    n1 = strtrim(sxpar(hdr, 'NAXIS1'))
    n2 = strtrim(sxpar(hdr, 'NAXIS2'))
    binning = [round(2700./n1), 1]
    ;; MMT red channel?
    IF strmatch(detector, '*ccd34*') THEN grating = 'R150'
    case grating of
        '800GPM': begin
            wstruct.pkwdth      = 7.00/binning[1]
            wstruct.TOLER       = 2.0D/binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.psig[0:4]   = [12.0, 5.0, 2.0, 2.0, 2.0]
            wstruct.mxoff[0:4]  =  [7.0, 7.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:4] = [3.0, 2.5, 2.5, 2.0, 1.75]
            wstruct.nord[0:4]        = [3L, 3L, 4L, 4L, 4L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 5L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            if strmatch(detector, '*mmtbluechan*') then begin
               ;; KHRR -- for new detector, 25 Oct 2014
               ;wstruct.AUTOID = 1
               wstruct.REID_FILE   = calib_path + '/mmt_blue_800_2015jan22.sav'
               wstruct.LINELIST    = line_path + '/mmt_blue_800_14oct25.lst'
            endif else begin
               if strmatch(sxpar(hdr,'COMPLAMP'), 'HeNeAr*') then begin
                  wstruct.REID_FILE   = calib_path + '/mmt_blue_800_nocu.sav'
                  wstruct.LINELIST    = line_path+'/henear_low.lst'
               endif else begin
                  wstruct.REID_FILE   = calib_path + '/mmt_blue_800.sav'
                  wstruct.LINELIST    = line_path+'/cuhenear.lst'
               endelse
            endelse 
        END
        '300GPM': begin
            wstruct.pkwdth      = 8.0D/binning[1]
            wstruct.TOLER       = 3.0D/binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.sigrej[0:4] = [3.0, 3.0, 2.5, 2.5, 2.0]
            wstruct.nord[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            if strmatch(detector, '*mmtbluechan*') then begin
               ;; KHRR -- for new detector, 25 April 2014
               wstruct.psig[0:4]   = [12.0, 5.0, 2.0, 2.0, 2.0]
               wstruct.mxoff[0:4]  = [7.0, 7.0, 3.0, 2.0, 2.0]/binning[1]
               wstruct.REID_FILE   = calib_path + '/mmt_blue_300_14apr24.sav'
               wstruct.LINELIST    = line_path + '/mmt_blue_300_14apr24.lst'
            endif else begin
               wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 5.0, 5.0]
               wstruct.mxoff[0:4]  = [7.0, 7.0, 3.0, 3.0, 3.0]/binning[1]
               wstruct.REID_FILE   = calib_path + '/mmt_blue_300.sav'
               wstruct.LINELIST    = line_path + '/mmt_blue_300.lst'
            endelse
         END
        '500GPM': begin
            wstruct.pkwdth      = 7.0D/binning[1]
            wstruct.TOLER       = 2.0D/binning[1]
            wstruct.FUNC  = 'CHEBY'
            wstruct.sigrej[0:4] = [3.0, 2.5, 2.5, 2.0, 1.8]
            wstruct.nord[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.REID        = 1
            ; wstruct.AUTOID      = 1
            if strmatch(detector, '*mmtbluechan*') then begin
               ;; KHRR -- for new detector, 25 April 2014
               wstruct.psig[0:4]   = [12.0, 5.0, 2.0, 2.0, 2.0]
               wstruct.mxoff[0:4]  = [7.0, 7.0, 3.0, 2.0, 2.0]/binning[1]
               wstruct.REID_FILE   = calib_path + '/mmt_blue_500.sav'
               wstruct.LINELIST    = line_path + '/mmt_blue_500.lst'
            endif 
         END
        '1200GPM': begin
           wstruct.pkwdth        = 7.0/binning[1]
           wstruct.TOLER         = 2.0D/binning[1]
           wstruct.FUNC  = 'CHEBY'
           wstruct.psig[0:4]     = [12.0, 5.0, 2.0, 2.0, 2.0]
           wstruct.MXOFF[0:4]    = [7.0, 7.0, 3.0, 2.0, 1.0]/binning[1]
            ;; paramters for wavelength soln fits
           wstruct.sigrej[0:4]   = [3.0, 3.0, 2.5, 2.5, 2.5]
           wstruct.nord[0:4]     = [4L, 4L, 4L, 4L, 4L]
           wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
           wstruct.LINELIST      = line_path+'/mmt_blue_1200.lst'
           wstruct.nord_panic    = 2L
           wstruct.AUTOID      = 0
           wstruct.REID          = 1 ; reidentify using archived solutions
           wstruct.REID_FILE     = calib_path + '/mmt_blue_1200_4575_HeArNe.sav'
        END 
;; Case for upgraded red channel gratings, using sky lines
        '1200-7700': begin
            wstruct.pkwdth      = 8.0D/binning[1]
            wstruct.TOLER       = 3.0D/binning[1]
            wstruct.FUNC  = 'LEGEND'
            wstruct.psig_reid[0:4]   = [12.0, 8.0, 7.0, 5.0, 5.0]
            wstruct.mxoff_reid[0:4]  = [7.0, 7.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej_reid[0:4] = [3.0, 3.0, 2.5, 2.5, 2.5]
            wstruct.nord[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.nord_reid[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.FLG_QUAL_REID[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 10L
            wstruct.nord_panic  = 2L
            wstruct.AUTOID      = 0
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + '/mmt_red_1200_9200_NeAr.sav'
            wstruct.LINELIST    = line_path + '//mmt_red_1200_9200_NeAr.lst'
        END
        '600-6310': begin
            wstruct.pkwdth      = 8.0D/binning[1]
            wstruct.TOLER       = 3.0D/binning[1]
            wstruct.FUNC  = 'LEGEND'
            wstruct.psig_reid[0:4]   = [12.0, 8.0, 7.0, 5.0, 5.0]
            wstruct.mxoff_reid[0:4]  = [7.0, 7.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej_reid[0:4] = [3.0, 3.0, 2.5, 2.5, 2.5]
            wstruct.nord[0:4]        = [3L, 4L, 4L, 4L, 4L]
            wstruct.nord_reid[0:4]        = [3L, 4L, 4L, 4L, 4L]
            wstruct.FLG_QUAL_REID[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 3L
            wstruct.nord_panic  = 2L
            wstruct.AUTOID      = 0
            wstruct.REID        = 1
            wstruct.REID_FILE   = calib_path + '/mmt_rcs_600_6310.sav'
            wstruct.LINELIST    = line_path + '/dis_henear.lst'
        END
;; old red channel
        'R150': BEGIN
            wstruct.pkwdth        = 7.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/dis_henear.lst'
            wstruct.REID_FILE   = calib_path + '/mmt_red_150.sav'
        END
        ELSE: STOP
    ENDCASE
ENDIF ELSE IF strmatch(instrument, '*IMACS*') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'DISPERSR'))
    binning = long(strsplit(sxpar(hdr, 'BINNING'), 'x', /extract))
    CASE grating OF
        'Gri-200-15.0': BEGIN
            wstruct.pkwdth        = 7.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 5L, 5L]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 10L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/dis_henear.lst'
            wstruct.REID_FILE   = calib_path + '/imacs_red_grism200.sav'
        END
        'Gra-300-4.3': BEGIN
            wstruct.pkwdth        = 7.0/binning[1]
            wstruct.TOLER         = 2.0/binning[1]
            wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 6.0, 4.0]
            wstruct.MXOFF[0:5]    = [12.0, 10.0, 8.0, 6.0, 6.0, 6.0]/binning[1]
            wstruct.sigrej[0:5]   = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.FLG_QUAL[0:5] = [1L, 1L, 1L, 2L, 2L, 2L]
            file = strtrim(sxpar(hdr,'FILENAME'),2)
            len = strlen(file)
            ichip = long(strmid(file,len-1))
            if ichip LE 2 or (ichip EQ 5 or ichip EQ 6) then begin
               wstruct.nord[0:5]     = [3L, 3L, 4L, 4L, 4L, 4L]  ;; Need lower order as we go bluer on the chip
               wstruct.npanic        = 9L
            endif else if ichip EQ 4 or ichip EQ 7 then begin  ;; Bluest data
               wstruct.nord[0:5]     = [1L, 1L, 1L, 1L, 1L, 1L]  ;; Need 1st order (straight line)
               wstruct.psig[0:5]     = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
               wstruct.FLG_QUAL[0:5] = [2L, 2L, 2L, 2L, 2L, 2L]
               wstruct.npanic        = 3L
            endif else begin
               wstruct.nord[0:5]     = [3L, 3L, 3L, 3L, 3L, 3L]  ;; Need lower order as we go bluer on the chip
               wstruct.psig[0:5]     = [12.0, 10.0, 7.0, 7.0, 7.0, 7.0]
               wstruct.npanic        = 6L
            endelse
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/dis_henear.lst'
            wstruct.REID_FILE   = calib_path + '/imacs_red_grating300.sav'
            wstruct.bin_ratio = binning[1] 
        END
        ELSE: STOP
    ENDCASE
ENDIF ELSE IF strmatch(instrument, 'GMOS*') THEN BEGIN
    grating = strtrim(sxpar(hdr, 'GRATING'))
    binning = long(strsplit(sxpar(hdr, 'CCDSUM'), ' ', /extrac))
    binning = shift(binning,1) ;; Swap to conform to LRIS type (i.e. transpose)
    wave = float(strcompress(sxpar(hdr, 'GRWLEN'), /rem))
    CASE strmid(grating,0,4) OF
       'B600': BEGIN
          if binning[1] eq 4 then begin
             wstruct.pkwdth        = 2*11.0/binning[1] ;; allow some fat lines
             wstruct.TOLER         = 2*3.0/binning[1]
             wstruct.psig[0:2]     = [5.0, 5.0, 4.0]
             wstruct.MXOFF[0:2]    = 2*[3.0, 3.0, 3.0]/binning[1]
             wstruct.sigrej[0:2]   = [3.0, 3.0, 3.0]
             wstruct.nord[0:2]     = [3L, 3L, 3L]
             wstruct.FLG_QUAL[0:2] = [1L, 2L, 2L] ;, 1L, 1L, 2L, 2L, 2L]
             wstruct.npanic        = 4L
             wstruct.nord_panic    = 2L
             
             ;;MF Jan 2014 added possibility of a blue tilt
             if(wave lt 490) then wstruct.REID_FILE   = calib_path + '/GMOS_CuAr_B600_480.sav' $
             else wstruct.REID_FILE   = calib_path + '/GMOS_CuAr_B600.sav'
             wstruct.BIN_RATIO     = float(binning[1])/4L ;arxiv soln binned spectrally by 4
          endif else if binning[1] eq 2 then begin
              wstruct.pkwdth        = 12.0D/binning[1]
              wstruct.TOLER         = 4.0D/binning[1]
              wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
              wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
              ;; parameters for wavelength soln fits
              wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.5]
              wstruct.nord[0:4]     = [3L, 3L, 4L, 4L, 4L]
              wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
              wstruct.npanic        = 5L
              wstruct.nord_panic    = 2L
              wstruct.REID_FILE   = calib_path + '/GMOS_CuAr_B600_bin2_495_570.sav'
              wstruct.BIN_RATIO   = float(binning[1])/2L
           endif

           wstruct.REID          = 1
           wstruct.MXSHIFT       = 200
           wstruct.LINELIST  = line_path+'/GMOS_CuAr_high.lst' ;; Might want R400
           wstruct.radius      = 2.0                           ;; special for GMOS microslits
           wstruct.fweight     = 1L
;           stop
        END
;        'R150+_G5306': BEGIN
        'R150': BEGIN
            wstruct.pkwdth        = 11.0/binning[1] ;; allow some fat lines
            wstruct.TOLER         = 3.0/binning[1]
            wstruct.psig[0:2]     = [5.0, 5.0, 4.0]
            wstruct.MXOFF[0:2]    = [3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:2]   = [3.0, 3.0, 3.0]
            wstruct.nord[0:2]     = [3L, 3L, 3L]
            wstruct.FLG_QUAL[0:2] = [1L, 2L, 2L] ;, 1L, 1L, 2L, 2L, 2L]
            wstruct.npanic        = 4L
            wstruct.nord_panic    = 2L
            wstruct.REID          = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/GMOS_CuAr_R150.lst'
            wstruct.REID_FILE   = calib_path + '/GMOS_CuAr_R150_multislit.sav'
            wstruct.radius      = 2.0 ;; special for GMOS microslits
        END
;        'R831+_G5302': BEGIN
        'R831': BEGIN
            wstruct.pkwdth      = 8.0D/binning[1]
            wstruct.TOLER       = 3.0D/binning[1]
            wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 5.0, 5.0]
            wstruct.mxoff[0:4]  = [7.0, 7.0, 3.0, 3.0, 3.0]/binning[1]
            wstruct.sigrej[0:4] = [3.0, 3.0, 2.5, 2.5, 2.5]
            wstruct.nord[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.nord_reid[0:4]        = [4L, 4L, 5L, 5L, 5L]
            wstruct.FLG_QUAL_REID[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 15L
            wstruct.nord_panic  = 2L
            wstruct.nfind        = 50L
            wstruct.dr_wave     = [0.9, 1.1]*(0.34d0)*binning[1]
            wstruct.disp_guess   = 0.34d0*binning[1]
            wstruct.REID        = 1
            wstruct.MXSHIFT       = 200
            wstruct.LINELIST  = line_path+'/GMOS_CuAr_high.lst'
            wstruct.REID_FILE   = calib_path + '/GMOS_R831.sav'
            wstruct.FUNC  = 'CHEBY'
         END
;         'R400+_G5305': BEGIN
         'R400': BEGIN
            wstruct.pkwdth        = 9.0/binning[1] ;; allow some fat lines
            wstruct.TOLER         = 3.0/binning[1]
            wstruct.psig[0:4]   = [12.0, 8.0, 7.0, 4.0, 3.0]
            wstruct.mxoff[0:4]  = ([5.0, 5.0, 3.0, 3.0, 3.0]/binning[1]) > 1.5
            wstruct.sigrej[0:4] = [3.0, 3.0, 3.0, 3.0, 3.0]
            wstruct.nord[0:4]        = [3L, 3L, 3L, 3L, 3L]
            wstruct.nord_reid[0:4]        = [3L, 3L, 3L, 3L, 3L]
            wstruct.FLG_QUAL_REID[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.FLG_QUAL[0:4]    = [1L, 1L, 1L, 2L, 2L]
            wstruct.npanic      = 15L
            wstruct.nord_panic  = 2L
            wstruct.REID          = 1
            wstruct.LINELIST  = line_path+'/GMOS_CuAr_R400.lst'
            wstruct.REID_FILE   = calib_path + '/GMOS_CuAr_R400_multislit.sav'
            wstruct.radius      = 4.0 ;; special for GMOS microslits
            wstruct.BIN_RATIO     = float(binning[1])/2L ;arxiv soln binned spectrally by 4
        END

        ELSE: STOP
    ENDCASE
; B&C spectrograph at the Bok 2.3-meter telescope; added by Moustakas
; on 2011-Jun-08
 ENDIF ELSE IF strmatch(telescope, '*bok*') THEN BEGIN
    binning = [1, 1]
    grating = strcompress(sxpar(hdr, 'DISPERSE'), /rem)
    CASE grating OF
       '400/4800': BEGIN
          wstruct.pkwdth        = 3.0/binning[1]
          wstruct.TOLER         = 2.0/binning[1]
          wstruct.psig[0:1]     = [5.0, 3.0]
          wstruct.MXOFF[0:1]    = [2.0, 2.0]/binning[1]
          wstruct.sigrej[0:1]   = [2.0, 2.0]
          wstruct.nord[0:1]     = [5L, 5L]
          wstruct.FLG_QUAL[0:1] = [2L, 2L]
          wstruct.FUNC  = 'CHEBY'
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'bok_henear.lst'
          wstruct.REID_FILE     = calib_path + 'bok_bc_400.sav'
        END
       '400/4889': BEGIN
          wstruct.pkwdth        = 3.0/binning[1]
          wstruct.TOLER         = 2.0/binning[1]
          wstruct.psig[0:1]     = [5.0, 3.0]
          wstruct.MXOFF[0:1]    = [2.0, 2.0]/binning[1]
          wstruct.sigrej[0:1]   = [2.0, 2.0]
          wstruct.nord[0:1]     = [5L, 5L]
          wstruct.FLG_QUAL[0:1] = [2L, 2L]
          wstruct.FUNC  = 'CHEBY'
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST      = line_path+'bok_henear.lst'
          wstruct.REID_FILE     = calib_path + 'bok_bc_400.sav'
        END

        ELSE: STOP
    ENDCASE
 ENDIF ELSE IF strmatch(telescope, '*DuPont*') THEN BEGIN
    binning = [1, 1]
    ;grating = strcompress(sxpar(hdr, 'DISPERSE'), /rem)
    ;CASE grating OF
    ;   '400/4800': BEGIN
          wstruct.pkwdth        = 3.0/binning[1]
          wstruct.TOLER         = 2.0/binning[1]
          wstruct.psig[0:1]     = [5.0, 3.0]
          wstruct.MXOFF[0:1]    = [2.0, 2.0]/binning[1]
          wstruct.sigrej[0:1]   = [2.0, 2.0]
          wstruct.nord[0:1]     = [5L, 5L]
          wstruct.FLG_QUAL[0:1] = [2L, 2L]
          wstruct.FUNC  = 'CHEBY'
          wstruct.npanic        = 10L
          wstruct.nord_panic    = 2L
          wstruct.REID          = 1
          wstruct.MXSHIFT       = 200
          wstruct.LINELIST  = line_path+'/dis_henear.lst'
          wstruct.REID_FILE   = calib_path + '/dupont_bc.sav'
    ;    END
    ;    ELSE: STOP
    ;ENDCASE
ENDIF ELSE IF strmatch(instrument, 'OSIRIS') then begin
   grism = strcompress(sxpar(hdr, 'GRISM'), /rem)
   binning = [2,2] ;; hardcoded right now
   case grism of
      'R1000B' : begin
         ;; parameters for arc peak finding
         wstruct.pkwdth        = 12.0D/binning[1]
         wstruct.TOLER         = 4.0D/binning[1]
         wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
         wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
         ;; parameters for wavelength soln fits
         wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.5]
         wstruct.nord[0:4]     = [4L, 5L, 5L, 6L, 6L]
         wstruct.FLG_QUAL[0:4] = [1L, 1L, 2L, 2L, 2L]
         wstruct.LINELIST      = line_path+'/GTC_osiris_R1000B.lst'
         wstruct.npanic        = 5L
         wstruct.nord_panic    = 2L
         wstruct.REID          = 1 ; reidentify using archived solutions
         wstruct.REID_FILE     = calib_path + '/GTC_osiris_R1000B.sav'
      end
      'R2000B' : begin
         ;; parameters for arc peak finding
         wstruct.pkwdth        = 12.0D/binning[1]
         wstruct.TOLER         = 4.0D/binning[1]
         wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
         wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
         ;; parameters for wavelength soln fits
         wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.5]
         wstruct.nord[0:4]     = [4L, 5L, 5L, 5L, 5L]
         wstruct.FLG_QUAL[0:4] = [1L, 1L, 2L, 2L, 2L]
         wstruct.LINELIST      = line_path+'/GTC_osiris_R2000B.lst'
         wstruct.npanic        = 5L
         wstruct.nord_panic    = 2L
         wstruct.REID          = 1 ; reidentify using archived solutions
         wstruct.REID_FILE     = calib_path + '/GTC_osiris_R2000B.sav'
      end
      'R2500R' : begin
         ;; parameters for arc peak finding
         wstruct.pkwdth        = 12.0D/binning[1]
         wstruct.TOLER         = 4.0D/binning[1]
         wstruct.psig[0:4]     = [12.0, 8.0, 7.0, 4.0, 4.0]
         wstruct.mxoff[0:4]    = [7.0, 7.0, 5.0, 5.0, 5.0]/binning[1]
         ;; parameters for wavelength soln fits
         wstruct.sigrej[0:4]   = [3.0, 3.0, 3.0, 2.5, 2.5]
         wstruct.nord[0:4]     = [4L, 5L, 5L, 5L, 5L]
         wstruct.FLG_QUAL[0:4] = [1L, 1L, 1L, 2L, 2L]
         wstruct.LINELIST      = line_path+'/GTC_osiris_R2500R.lst'
         wstruct.npanic        = 5L
         wstruct.nord_panic    = 2L
         wstruct.REID          = 1 ; reidentify using archived solutions
         wstruct.REID_FILE     = calib_path + '/GTC_osiris_R2500R.sav'
      end
      else: stop
   endcase
ENDIF ELSE IF strmatch(instrument, 'ISAAC') then begin
   mode = strcompress(esopar(hdr, 'HIERARCH ESO INS GRAT NAME'), /rem)
   order = strcompress(long(esopar(hdr, 'HIERARCH ESO INS GRAT ORDER')), /rem)
   slit_str = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
   wave_cen =  esopar(hdr, 'HIERARCH ESO INS GRAT WLEN')
   CASE slit_str OF 
      'slit_1': slit_width = 1.0d
      'slit_0.6_tilted': slit_width = 0.6d
      ELSE: message, 'Unrecognized slit'
   ENDCASE
   CASE order OF 
      '2': band = 'K'
      '3': band = 'H'
      ELSE: message, 'unrecognized order'
   ENDCASE
   plate_scale = 0.147d
   fnslit = slit_width/plate_scale
   wstruct.pkwdth = 1.3d*fnslit
   wstruct.toler = fnslit/2.0d > 2.0d
   wstruct.THIN = 1
   wstruct.FWEIGHT = 0
   wstruct.PSIG[*] = [10.0, 7.0, 5.0, 5.0, 5.0, 5.0]
   wstruct.mxoff[*] = 5.0
   wstruct.sigrej[*] = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
   wstruct.nord[*] = 4
   wstruct.linelist = line_path + $
                      '/ISAAC_modelsky_OH_linelist_' + mode + '_' +band + '.lst'
   wstruct.npanic = 3L
   wstruct.nord_panic = 2L
   wstruct.REID = 1
   wstruct.REID_FILE = calib_path + '/ISAAC_' + mode + '_' + band + '.sav'
   wstruct.radius = 5L
   wstruct.flg_qual[*] = 1L
   wstruct.wave_cen = wave_cen
   wstruct.band = band
   wstruct.mode = mode
ENDIF ELSE IF strmatch(instrument, 'Lucifer') OR strmatch(instrument, 'LUCI2') then begin
   gratname = strcompress(sxpar(hdr, 'GRATNAME'), /rem)
   slit_str   = strcompress(sxpar(hdr, 'MASKID'), /rem)
   wave_cen =  sxpar(hdr, 'GRATWLEN')
   CASE slit_str OF
      '990065': slit_width = 0.25d
      '990078': slit_width = 0.50d
      '990029': slit_width = 0.75d
      '990034': slit_width = 1.00d
      'LS0.5_300mue':  slit_width = 0.50d
      ELSE: message, 'ERROR: Unknown mask'
   ENDCASE
   plate_scale = 0.25D
   fnslit = slit_width/plate_scale
   wstruct.pkwdth = 1.3d*fnslit
   wstruct.toler = fnslit/2.0d > 2.0d
   wstruct.THIN = 1
   wstruct.FWEIGHT = 0
   wstruct.PSIG[*] = [10.0, 7.0, 5.0, 5.0, 5.0, 5.0]
   wstruct.mxoff[*] = 5.0
   wstruct.sigrej[*] = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
   wstruct.nord[*] = 4
   wstruct.linelist = line_path + $
                      '/LUCI_modelsky_OH_linelist_' + gratname + '.lst'
   wstruct.npanic = 3L
   wstruct.nord_panic = 2L
   wstruct.REID = 1
   wstruct.REID_FILE = calib_path + '/luci_' + gratname + '.sav'
   wstruct.radius = 5L
   wstruct.flg_qual[*] = 1L
   wstruct.wave_cen = wave_cen
ENDIF ELSE IF strmatch(instrument, 'SOFI') then begin
   slit_str = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 ID'), /rem)
   slit_width = double(strmid(slit_str, 10, 3))
   plate_scale = 0.288d
   fnslit = slit_width/plate_scale
   setup = strcompress(esopar(hdr, 'HIERARCH ESO INS MODE'), /rem)
   CASE setup OF
      'LONG_SLIT_RED': BEGIN
         wstruct.linelist = line_path + '/SOFI_modelsky_OH_linelist_H+K.lst' 
         wstruct.REID_FILE = calib_path + '/SOFI_H+K.sav'
         wstruct.BAND = 'H+K'
         wstruct.nord[*] = 4
      END
      'LONG_SLIT_H': BEGIN 
         wstruct.linelist = line_path + '/SOFI_modelsky_OH_linelist_H.lst'
         wstruct.REID_FILE = calib_path + '/SOFI_H.sav'
         ;;wstruct.linelist = line_path + '/oh_linelist.lst' 
         ;;wstruct.REID_FILE = calib_path + '/sofi_H.sav'
         wstruct.BAND = 'H'
         wstruct.nord[*] = [3.0, 3.0, 3.0, 4.0, 4.0, 4.0]
      END
      'LONG_SLIT_K': BEGIN 
         wstruct.linelist = line_path + '/SOFI_modelsky_OH_linelist_K.lst' 
         wstruct.REID_FILE = calib_path + '/SOFI_K.sav'
         wstruct.BAND = 'K'
         wstruct.nord[*] = 3
      END
      ELSE: message, 'unrecognized setup'
   ENDCASE
   wstruct.pkwdth = 1.3d*fnslit
   wstruct.toler = fnslit/2.0d > 2.0d
   wstruct.THIN = 1
   wstruct.FWEIGHT = 0
   wstruct.ICLSE = 2.0
   wstruct.PSIG[*] = [10.0, 7.0, 5.0, 5.0, 4.0, 4.0] ;; go down to 4-sigma 
   wstruct.mxoff[*] = 5.0
   wstruct.sigrej[*] = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
   wstruct.npanic = 3L
   wstruct.nord_panic = 2L
   wstruct.REID = 1
   wstruct.radius = 5L
   wstruct.flg_qual[*] = 1L
ENDIF ELSE IF strmatch(instrument, 'SINFONI') then begin
   optic = strcompress(esopar(hdr, 'HIERARCH ESO INS OPTI1 NAME'), /rem)
   grating = strcompress(esopar(hdr, 'HIERARCH ESO INS GRAT1 NAME'), /rem)
   slit_width = double(optic)
   plate_scale = slit_width/2.0d
   CASE optic OF
      '0.025': BEGIN
         ;; large radius for 0.025 as we need signal
         wstruct.radius = 15L
      END
      '0.10': BEGIN
         message, 'Not yet supported'
      END
      '0.25': BEGIN
         wstruct.radius = 5L
      END
      ELSE: message, 'ERROR: Unknown setup'
   ENDCASE
   fnslit = slit_width/plate_scale
   wstruct.pkwdth = 1.3d*fnslit
   wstruct.toler = fnslit/2.0d > 2.0d
   wstruct.THIN = 1
   wstruct.FWEIGHT = 0
   wstruct.PSIG[*] = [10.0, 7.0, 5.0, 5.0, 5.0, 5.0]
   wstruct.mxoff[*] = 5.0
   wstruct.sigrej[*] = [3.0, 3.0, 3.0, 3.0, 3.0, 3.0]
   wstruct.nord[*] = 4
   wstruct.REID_FILE = calib_path + $
                       '/SINFONI_' + optic + '_' + grating +'.sav'
   wstruct.linelist = line_path + $
                      '/SINFONI_modelsky_OH_linelist_'  + grating + '.lst'
   wstruct.npanic = 3L
   wstruct.nord_panic = 2L
   wstruct.REID = 1
   wstruct.flg_qual[*] = 1L
   ;;wstruct.wave_cen = wave_cen
endif else BEGIN
   message, 'ERROR: Unknown instrument!'
ENDELSE

if keyword_set(LINELIST) then begin
   wstruct.LINELIST = line_path+LINELIST
   ;; Kludge for no Ne lamps
   ipos = strpos(LINELIST,'noNe')
   if ipos GT 0 then wstruct.npanic = 5L
endif

if keyword_set(REID_FILE) then wstruct.REID_FILE = calib_path+REID_FILE
if keyword_set(BIN_RATIO) then wstruct.BIN_RATIO = BIN_RATIO

splog, 'Wavelength solution using linelist:', wstruct.LINELIST, ' and archive file:' $
       , wstruct.REID_FILE
RETURN, wstruct
END

