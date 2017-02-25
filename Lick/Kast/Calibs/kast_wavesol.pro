;+ 
; NAME:
; kast_wavesol   
;   Version 1.1
;
; PURPOSE:
;  Create a wavelength solution for a specific object using 
;     an archived solution.  The wavelength array in the object
;     structure and writes that structure out.
;
; CALLING SEQUENCE:
; kast_wavesol, kast, setup, side, obj_id, [exp], /STD,
;  /SILENT, /AINTER, /AUTO, CALIB=, CALIBFIL=, NFIN= ,ARC_SAT = 
;
; INPUTS:
;   kast  --  Kast IDL structure
;  setup  --  Setup value
;   side  --  Specific camera [blue (1) vs. red (2)]
; obj_id  --  Object value
;  [exp]  --  Exposure indices
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;  /AUTO    -- Automatically fit the arc
;  NFIN=    -- Requires nfin lines to start the fit [default: 5]
;  /SILENT
;  /AINTER  -- Peak up on the arc cross-correlation interactively
;  ARC_SAT= -- max counts for arc lines for determining which to use
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   01-Mar-2003 Written by JXP
;   12-Nov-2005 Added defaults for Blue G2 (600/4310) and Red 600/7500
;               and ARC_SAT keyword, KLC
;------------------------------------------------------------------------------
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro kast_wavesol, kast, setup, side, obj_id, exp, STD=std, $
                  SILENT=silent, AINTER=ainter, AUTO=auto, FCALIB=calib, $
                  CALIBFIL=calibfil, NFIN=nfin,ARC_SAT=arc_sat, ARCSP=arcsp, $
                  TEMPLT=templt, LINLIST=linlist

;
if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
      'kast_wavesol, kast, setup, side, obj_id, [exp], /STD, /SILENT, /AINTER'
    print, '   /AUTO, FCALIB=, CALIBFIL=, NFIN=,  ARCSP= [v1.1]'
    return
endif 

;  Optional Keywords
if not keyword_set( NFIN ) then nfin = 5L
if not keyword_set( SETUP ) then setup = 0
if not keyword_set( SIDE ) then side = 1
if not keyword_set( OBJ_ID ) then obj_id = 0
if not keyword_set(ARC_SAT) then arc_sat = 20000.

; Find objects
if not keyword_set( STD ) then begin
    indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                 kast.side EQ side AND kast.setup EQ setup AND $
                 kast.obj_id EQ obj_id AND kast.type EQ 'OBJ', nindx)
    if nindx EQ 0 then begin
        print, 'kast_wavesol: No images to find obj for!', obj_id
        return
    endif
endif else begin                ; STANDARD STAR
    indx = where(kast.flg_anly NE 0 AND kast.mode EQ 1 AND $
                 kast.side EQ side AND kast.setup EQ setup AND $
                 kast.type EQ 'STD', nindx)
    if  N_params() LT 4  OR n_elements(OBJ_ID) NE 1 then begin 
        print,'Syntax - ' + $
          'kast_wavesol, kast, setup, side, EXP, /STD   [v1.0]'
        return
    endif 
    indx = indx[obj_id]
    nindx = 1L
endelse

;  Exposures
if not keyword_set(exp) then exp = lindgen(nindx)

if not keyword_set(templt) and not keyword_set(linlist) then begin
    if side EQ 1 then begin
        case strtrim(kast[indx[0]].grising,2) of 
            '452/3306': begin
                templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g1.idl'
                linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue3.lst'
            end
            '600/4310': begin
                templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g2.idl'
                linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue2.lst'
            end
            '830/3460': begin
                templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_g3.idl'
                linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_blue2.lst'
            end 
            else: stop
        endcase
    endif else begin
        linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/kast_red.lst'
        case strtrim(kast[indx[0]].grising,2) of 
            '300/7500': $
              templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_3007500.idl'
            '600/7500': $
              templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_6007500.idl'
            '1200/5000': $
              templt = getenv('XIDL_DIR')+'/Lick/Kast/Calibs/kastfit_12005000.idl'
            else: stop
        endcase
    endelse
endif

;  Loop
for q=0L,n_elements(exp)-1 do begin

    ;; Open objfil
    objfil = kast[indx[exp[q]]].obj_fil 
    if x_chkfil(objfil+'*') EQ 0 then begin
        print, 'kast_extract: No Obj file ', objfil
        continue
    endif
    objstr = xmrdfits(objfil, 1, STRUCTYP='specobjstrct', /silent)

    ;; Nobj
    nobj = n_elements(objstr)
    if keyword_set( STD) then nobj = 1

    ;; Read IMG+VAR
    imgfil = 'Final/f_'+kast[indx[exp[q]]].img_root
    if x_chkfil(imgfil+'*') EQ 0 then begin
        print, 'kast_extract: No Final file ', imgfil, '  Continuing...'
        continue
    endif
    img = xmrdfits(imgfil, /silent)
    var = xmrdfits(imgfil, 1, /silent)
    sz = size(img, /dimensions)

    ;; Read arc
    arc_fil = kast[indx[exp[q]]].arc_fil
    if x_chkfil(arc_fil+'*') EQ 0 then begin
        print, 'kast_extract: No Arc file ', arc_fil, '  Continuing...'
        continue
    endif
    arc = x_readimg(arc_fil, /fscale)

    for kk=0L,nobj-1 do begin
        ;; Extract
        arcsp = x_extract(arc, [objstr[kk].ycen-7.,objstr[kk].ycen+7.], $
                          objstr[kk].trace[0:sz[0]-1], CAPER=objstr[kk].xcen)

        ;; CALIBRATE
        if keyword_set( AUTO ) then begin
            ;; Grab template
            restore, templt
            ;; Cross-correlate
            step = lindgen(400L) - 200L
            corr = c_correlate((aspec<ARC_SAT), (arcsp<ARC_SAT), step, /double)
            mx = max(corr, imx)
            imx = step[imx]
            print, 'kast_wavesol: Offseting ', strtrim(imx,2), ' pixels'
            ;; Line list
            x_arclist, linlist, Lines
            ;; ID
            x_templarc, arcsp, lines, calib, SHFT=imx
            ;; Peak up
                                ;stop
            x_arcpeakup, arcsp, lines, FFIT=ffit, WV=wave, INTER=ainter, $
              NFIN=nfin
        endif else begin
            x_identify, arcsp, calib, LINELIST=linlist
            wave = x_calcfit(n_elements(arcsp), FITSTR=calib)
            if keyword_set( CALIBFIL ) then begin
                aspec = arcsp
                save, aspec, calib, filename=calibfil
            endif
        endelse

        ;; VACUUM
        if not keyword_set( NOVAC ) then airtovac, wave
        ;; Save
        objstr[kk].wave[0:sz[0]-1] = float(wave)
    endfor
    ;; Write
    print, 'kast_wavesol: Updating ', objfil
    mwrfits, objstr, objfil, /create
    spawn, 'gzip -f '+objfil

endfor


print, 'kast_extract: All Done!'
return
end
