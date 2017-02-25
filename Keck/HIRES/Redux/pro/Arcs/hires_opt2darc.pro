;+ 
; NAME:
; hires_opt2darc
;     Version 1.0
;
; PURPOSE:
;   Iteratively optimise the 2D wavelength solution. Can be run
;   multiple times to improve the solution. Can be used to identify
;   lines in sporadic orders where hires_fitarc fails. It performs 
;   the following steps (or more accurately, x_opt2darc does):
;   
;    1.  Identifies potential arc lines per order, assigns
;        wavelengths based on the current 2D wavelength solution
;    2.  Matches lines in LINLIST with potential arc lines per
;        order according to a deta v < MXVELOFF criterion
;    3.  Performs initial sigma rejection using 
;        djs_iterstat where sigrej = SIGPRE
;    4.  Optionally, performs iterative sigma rejection with SIGREJ
;        as the key parameter.
;    5.  Optionally, iteratively crops lines from the fit with a 
;        delta v > VELCUT criterion and refits every iteration.
;        Essentially, allows you to define a desired maximum delta v for 
;        the arc lines about your fit at the expense of fitting less
;        lines.
;    6.  Outputs QA files (overwrites those output by hires_fit2darc),
;        including additional diagnostics.
;   
;   Steps 3, 4, and 5 above are technically optional but run by default.
;   If results aren't up to scratch, try re-running this routine. If 
;   results are still suboptimal experiment with the optional parameters.
;   It may be worthwhile to backup your 2D arc solutions since this
;   routine naturally overwrites them.
;
; CALLING SEQUENCE:
;   
;  hires_opt2darc, hires, setup, chip
;
; INPUTS:
;   hires   -  HIRES structure
;   setup   -  Setup ID
;   chip    -  Blue (1), Green (2), Red (3), or multiple [1,2L]
;              (Default: [1,2,3L])
;
; RETURNS:
;
; OUTPUTS:
;   Overwrites the 2D arc solution 
;   (e.g. Arcs/Fits/Arc_B0192_fit2D.fits)
;
; OPTIONAL KEYWORDS:
;   OBJ_ID    - Object ID
;   NOCOEFF   - Number of coefficients to use in the x-direction
;               (default varies with chip, see source code)
;   NYCOEFF   - Number of coefficients to use in the y-direction
;               (default varies with chip, see source code)
;   PKWDTH    - Peak width parameter passed to x_fndpeaks in step 1
;   PKSIG     - Peak detection sigma threshold parameter passed to
;               x_fndpeaks in step 1
;   MXVELOFF  - Maximum velocity offset parameter in step 2 [m/s]
;               (default is 500.0)
;   SIGPRE    - Sigma for initial line rejection in step 3
;               (default is 2.5)
;   SIGREJ    - Sigma for iterative line rejection in step 4
;               (default is 2.0)
;   VELCUT    - Iterative velocity cut in step 5
;               (default varies with chip, see source code)
;   LINLIST   - Arc line list
;               (default is $XIDL_DIR/Spec/Arcs/Lists/hires_thar.lst)
;   DRY       - Dry run, hires_opt2darc won't write anything - useful
;               for experimenting with parameters
;   
; Optional OUTPUTS:
;    None
; 
; COMMENTS:
;  The optional parameter defaults will suffice in most cases. Note
;  they were arrived at through experimentation on a variety RED
;  cross disperser XDA and ECA 2x1 binned ThAr exposures. 
;
; EXAMPLES:
;   hires_opt2darc, hires, setup, LINLIST='data/hiredux/thar_wav_opt.dat'
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_opt2darc
;
; REVISION HISTORY:
;   23-Nov-2011 Written by Adrian L. Malec
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro hires_opt2darc, hires, setup, chip, OBJ_ID=obj_id, $
		nycoeff=nycoeff, nocoeff=nocoeff, $
		PKWDTH=pkwdth, PKSIG=pksig, $
		MXVELOFF=mxveloff, SIGPRE=sigpre, SIGREJ=sigrej, VELCUT=velcut, $
		LINLIST=linlist, DRY=dry, _EXTRA=extra
	
	; defaults
	if not keyword_set(chip) then chip = [1L,2L,3L] ; OPTk
	if not keyword_set(linlist) then $
		linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' 
	print, 'hires_opt2darc: Using line list -- ', linlist
	if not keyword_set(mxveloff) then mxveloff = 500.0
	if not keyword_set(sigpre) then sigpre = 2.5
	if not keyword_set(sigrej) then sigrej = 2.0
	
	for i=0, n_elements(chip)-1 do begin
		ii = chip[i]
		
		; chip-dependant defaults
		if not keyword_set(nocoeff) then begin
			if ii eq 1L then nocoeff = 7
			if ii eq 2L then nocoeff = 6
			if ii eq 3L then nocoeff = 4
		endif
		
		if not keyword_set(nycoeff) then begin
			if ii eq 1L then nycoeff = 6
			if ii eq 2L then nycoeff = 7
			if ii eq 3L then nycoeff = 5
		endif
		
		if not keyword_set(velcut) then begin
			if ii eq 1L then velcut = 60.0
			if ii eq 2L then velcut = 80.0
			if ii eq 3L then velcut = 100.0
		endif
		
		if keyword_set(OBJ_ID) then begin
			obj = where((hires.type EQ 'OBJ' OR hires.type EQ 'STD') $
				AND hires.flg_anly NE 0 $
				AND hires.setup EQ setup and $
				hires.obj_id EQ obj_id $
				AND hires.chip EQ ii, narc)
			dumfil = hires_getarcfil(hires, obj, raw_fil=arcfil)
			idx = obj
		endif else begin
			arcs = where(hires.type EQ 'ARC' AND hires.flg_anly NE 0 $
				AND hires.setup EQ setup $
				AND hires.chip EQ ii, narc)
			arcfil= strtrim(hires[arcs].rootpth,2) + strtrim(hires[arcs].img_root,2)
			idx = arcs
		endelse
		
		narc = n_elements(arcfil)
		
		for qq=0L, narc-1 do begin
			
			; maybe this will work.. maybe not...
			pos = strpos(arcfil[qq], '.fits')
			frame = long(strmid(arcfil[qq],pos-4,4))
			
			; hires[idx[qq]].setup , chip, frame? ??
			arc_fil = hires_getfil('arc_fil', FRAME=frame, CHIP=ii, /name, CHKFIL=chkfil)
			arc_2dfit = hires_getfil('arc_2Dfit', setup, CHIP=ii, CHKFIL=chkf, FRAME=frame)
			ordr_str = hires_getfil('ordr_str', setup, CHIP=ii, fil_nm=ordr_fil)
			qatxt = hires_getfil('qa_arc2dtxt', setup, CHIP=ii, /name, CHKFIL=chkf, FRAME=frame)
			qawv = hires_getfil('qa_arc2dwv', setup, CHIP=ii, /name, CHKFIL=chkf, FRAME=frame)
			qafil = hires_getfil('qa_arc2dfit', setup, CHIP=ii, /name, CHKFIL=chkf, FRAME=frame)
			
			if not keyword_set( DRY ) then begin
				x_opt2darc, arc_2dfit, arc_fil, ordr_str, nycoeff, nocoeff, $
					MXVELOFF=mxveloff, PKWDTH=pkwdth, PKSIG=pksig, $
					SIGPRE=sigpre, SIGREJ=sigrej, VELCUT=velcut, LINLIST=linlist, $
					OUT_STR=out_str, QATXT=qatxt, QAWV=qawv, QAFIL=qafil, _EXTRA=extra
			endif else begin
				x_opt2darc, arc_2dfit, arc_fil, ordr_str, nycoeff, nocoeff, $
					MXVELOFF=mxveloff, PKWDTH=pkwdth, PKSIG=pksig, $
					SIGPRE=sigpre, SIGREJ=sigrej, VELCUT=velcut, LINLIST=linlist, $
					OUT_STR=out_str, _EXTRA=extra ; no QA
			endelse
			
			if not keyword_set( DRY ) then begin ; a dry run, don't overwrite 2D fit file
				outfil = hires_getfil('arc_2Dfit', setup, CHIP=ii, /name, CHKFIL=chkf, FRAME=frame)
				print, 'hires_opt2darc: Writing ' + outfil
				mwrfits, out_str, outfil, /create
			endif
			
		endfor
		
	endfor
	
	print, 'hires_opt2darc: All done!!!'
	
end

