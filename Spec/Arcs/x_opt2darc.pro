;+ 
; NAME:
; x_opt2darc
;     Version 1.0
;
; PURPOSE:
;   Does the work for hires_opt2darc. See the hires_opt2darc 
;   documentation for details of how this routine works.
;
; CALLING SEQUENCE:
;
; INPUTS:
;
; RETURNS:
;
; OUTPUTS:
;
; OPTIONAL KEYWORDS:
;   
; Optional OUTPUTS:
; 
; COMMENTS:
;
; EXAMPLES:
;
; PROCEDURES/FUNCTIONS CALLED:
;  x_arclist
;  x_fndpeaks
;  djs_iterstat
;  djs_reject
;
; REVISION HISTORY:
;   23-Nov-2011 Written by Adrian L. Malec
;-
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_fit2dlines, fitln, nycoeff, nocoeff, MSK=msk, res=res, nrmt=nrmt, nrmp=nrmp

	npix = n_elements(fitln)
	t = fitln.order
	all_pix = fitln.pix
	all_wv = fitln.wave * t ; note x t (order number)
	
	; msk is a line list mask, fit all by default. Note: result still outputs fit for all lines
	if not keyword_set( MSK ) then msk = bytarr(npix) + 1B
	
	nrmp = dblarr(2)
	mnx = min(all_pix, MAX=mxx)
	nrmp[0] = 0.5 * (mnx + mxx)
	nrmp[1] = mxx - mnx
	pix_nrm = 2. * (all_pix - nrmp[0])/nrmp[1]
	
	nrmt = dblarr(2)
	mnx = min(t, MAX=mxx)
	nrmt[0] = 0.5 * (mnx + mxx)
	nrmt[1] = mxx - mnx
	t_nrm = 2. * (t - nrmt[0])/nrmt[1]
	
	invvar = replicate(1., npix) ; ^-1 variances set to 1
	
	work2d = dblarr(npix, nycoeff*nocoeff)
	worky = flegendre(pix_nrm[*], nycoeff)
	workt = flegendre(t_nrm[*], nocoeff)
	
	for i=0, nocoeff-1 do begin
		for j=0, nycoeff-1 do begin
			work2d[*, j*nocoeff+i] = worky[*, j] * workt[*, i]
		endfor
	endfor
	
	work2di = transpose(work2d[where(msk), *] * (invvar[where(msk)] # replicate(1,nocoeff*nycoeff)))
	alpha = work2di # work2d[where(msk), *]
	beta = work2di # all_wv[where(msk)]
	choldc, alpha, p, /double
	res = cholsol(alpha, p, beta, /double)
	
	wv_mod = dblarr(npix)
	wv_mod[*] = (work2d # res) / t
	
	fitln.wave_fit = wv_mod

end

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

pro x_opt2darc, arc_2dfit, arcfil, ordr_str, nycoeff, nocoeff, $
	PKSIG=pksig, PKWDTH=pkwdth, SATUR=satur, $
	MXVELOFF=mxveloff, SIGPRE=sigpre, SIGREJ=sigrej, VELCUT=velcut, $
	LINLIST=linlist, QAFIL=qafil, QATXT=qatxt, QAWV=qawv, OUT_STR=out_str
	
	c_ms = 299792458.0
	
	no = arc_2dfit.no
	ny = arc_2dfit.ny
	nrm = arc_2dfit.nrm
	nrmt = arc_2dfit.nrmt
	res = arc_2dfit.res
	
	fitln_struct = {wave: 0d, pix: 0.0, wave_fit: 0d, name: 'XX', order: 0}
	
	nordr = n_elements(ordr_str)
	gd_ordr = where(ordr_str.flg_anly NE (-1), ngd_ordr)
	head = xheadfits(arcfil, /silent)
	sz = lonarr(2)
	sz[0] = sxpar(head,'NAXIS1')
	sz[1] = sxpar(head,'NAXIS2')
	
	;; Open master line list and sort
	if not keyword_set( LINLIST ) then linlist = getenv('XIDL_DIR')+'/Spec/Arcs/Lists/hires_thar.lst' ; OPTk
	x_arclist, linlist, lines
	print, 'x_opt2darc: Using line list ', linlist, ' ['+string(n_elements(lines),format='(I5)')+' lines]'
	srt = sort(lines.wave)
	lines = lines[srt]
	
	;; Grab Arc IMG
	if x_chkfil(arcfil+'*',/silent) EQ 0 then begin
		print, 'x_opt2darc: Arc ', arcfil, ' does not exist!'
		return
	endif
	print, 'x_opt2darc: Reading arc: ', arcfil
	; Extract the arc (pseudo-boxcar)
	arc_img = xmrdfits(arcfil, 0, head, /silent)
	sz_arc = size(arc_img, /dimensions) ; should be the same as sz
	aspec = fltarr(sz_arc[1], 50)
	aspec[*,0:nordr-1] =  x_extractarc(arc_img, ordr_str)
	
	bound_nrm = 2. * ([0, round(sz[1])-1.0] - nrm[0])/nrm[1]
    bound_worky = flegendre(bound_nrm[*], ny)
	
	print, 'x_opt2darc: Identifying lines...'
	; Identify lines per order
	fitln_1 = 0
	for j=0L, ngd_ordr-1 do begin
		i = gd_ordr[j]
		ord = ordr_str[i].order
		
		; Find order bounds
		bound_t_nrm = 2. * (replicate(float(ord), 2) - nrmt[0])/nrmt[1]
		bound_work2d = dblarr(2, ny*no)
		bound_workt = flegendre(bound_t_nrm[*], no)
		for ii=0, no-1 do begin
			for jj=0, ny-1 do begin
				bound_work2d[*, jj*no+ii] = bound_worky[*, jj] * bound_workt[*, ii]
			endfor
		endfor
		bound_wv = dblarr(2)
		bound_wv[*] = bound_work2d # res / replicate(float(ord), 2)
		
		; Min and max wavelengths for order
		wvmn = min(bound_wv, max=wvmx)
		
		msk = bytarr(sz_arc[1])
		msk[where(aspec[*,i] GT 0)] = 1B

		;; SAT
		if not keyword_set( SATUR ) then satur = 45000; OPTk
		satpix = (aspec[*,i] GT SATUR)
		sat_region = 25*sz_arc[1]/2048L
		sat = where(smooth(1.0*satpix,sat_region) GT 0, nsat) 
		if nsat NE 0 then msk[sat] = 0B
		
		if not keyword_set( PKWDTH ) then pkwdth = (6L / round(3990./sz[1])) > 2L
		if not keyword_set( PKSIG ) then pksig = 7. ; OPTk
		
		x_fndpeaks, aspec[*,i], peak, NSIG = pksig, /silent, PKWDTH = pkwdth, $
			THIN = thin, MSK = msk, NORDB = 9, $
			TOLER = 2.0D, /FGAUSS
		if peak[0] EQ -1 then begin
			print, '   No peaks found for order '+string(ord,format='(I3)')+'! Skipping..' ; specify which order and continue looping
			continue
		endif
		
		; Find candidate lines for order in line list
		gdlin = where(lines.wave GT wvmn AND $
			lines.wave LT wvmx AND lines.flg_qual GT 0 AND $
			lines.flg_qual LE 99999L, ngdlin)
		
		if (ngdlin eq 0L) then begin
			print, '   No lines found for order '+string(ord,format='(I3)')+' in line list! Skipping..' ; specify which order and continue looping
			continue
		endif
		
		ordlines = lines[gdlin]
		
		; Make wave (wv) array for peaks
		npk = n_elements(peak)
		pix_nrm = 2. * (peak - nrm[0])/nrm[1]
		worky = flegendre(pix_nrm[*], ny)		

		tsub = replicate(float(ord), npk)
		t_nrm = 2. * (tsub - nrmt[0])/nrmt[1]
		
		work2d = dblarr(npk, ny*no)
		workt = flegendre(t_nrm[*], no)
		
		for ii=0, no-1 do begin
			for jj=0, ny-1 do begin
				work2d[*, jj*no+ii] = worky[*, jj] * workt[*, ii]
			endfor
		endfor
		
		wv = dblarr(npk)
		wv[*] = work2d # res / tsub
		
		if not keyword_set( MXVELOFF ) then mxveloff = 600. ; OPTk (m/s)
		
		; Identify lines - loop on line list subset for order and try to match with peaks
		; in arc spectrum according to detav < MXVELOFF criterion
		for ii=0L, ngdlin-1 do begin
			lnwv = replicate(ordlines[ii].wave, npk)
			deltav = abs(lnwv - wv)/lnwv*c_ms
			mindv = min(deltav, imin)
			if mindv le mxveloff then begin
				ordlines[ii].flg_plt = 1
				ordlines[ii].pix = peak[imin]
			endif
		endfor
		
		id = where(ordlines.flg_plt eq 1)
		if id[0] NE -1 then begin
			print, '   '+string(n_elements(id),format='(I3)')+' lines identified for order '+string(ord,format='(I3)')
			; append to final fit list
			idln = ordlines[id]
			fitln_i = replicate(fitln_struct, n_elements(id))
			fitln_i.wave = idln.wave
			fitln_i.pix = idln.pix
			fitln_i.name = idln.name
			fitln_i.order = ord
			if not fitln_1 then begin
				fitln = [fitln_i]
				fitln_1 = 1
			endif else fitln = [fitln, fitln_i]
		endif else begin
			print, '   No lines identified for order '+string(ord,format='(I3)')+'!'
		endelse
	
	endfor ; end order loop
	
	nln = n_elements(fitln)
	print, 'x_opt2darc: '+string(nln,format='(I4)')+' lines identified in total.'
	x_fit2dlines, fitln, nycoeff, nocoeff
	vrms = sqrt(total(((fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms)^2)/nln)
	print, '   Delta velocity rms (m/s) of '+ string(vrms,format='(F10.1)')
	
	if keyword_set( SIGPRE ) then begin
		print, 'x_opt2darc: Preliminary sigma rejection'
		djs_iterstat, (fitln.wave_fit - fitln.wave), sigrej=sigpre, mask=sigmsk
		gdsig = where(sigmsk EQ 1B)
		fitln = fitln[gdsig]
		nln = n_elements(fitln)
		vrms = sqrt(total(((fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms)^2)/nln)
		print, 'x_opt2darc: '+ string(nln,format='(I4)') + ' lines with delta velocity rms (m/s) of '+string(vrms,format='(F10.1)')
	endif
	
	if keyword_set( SIGREJ ) then begin
		print, 'x_opt2darc: Iterative sigma rejection'
		iiter = 0
		maxiter = 50
		prev_tnln = 0
		while (NOT keyword_set(qdone) AND iiter LE maxiter) do begin
			qdone = djs_reject(fitln.wave, fitln.wave_fit, outmask=outmask, upper=sigrej, lower=sigrej)
			x_fit2dlines, fitln, nycoeff, nocoeff, msk=outmask
			iiter = iiter + 1
			
			tln = fitln[where(outmask)]
			tnln = n_elements(tln)
			vrms = sqrt(total(((tln.wave_fit - tln.wave)/tln.wave_fit*c_ms)^2)/tnln)
			print, '   ' + string(iiter,format='(I3)') + ': ' + string(tnln,format='(I4)') + ' lines with delta velocity rms (m/s) of '+ string(vrms,format='(F10.1)')
			
			if tnln eq prev_tnln then begin
				print, 'x_opt2darc: Breaking... no lines rejected in previous iteration.'
				break
			endif
			if iiter eq maxiter then begin
				print, 'x_opt2darc: Breaking... maximum iterations reached.'
				break
			endif
			if tnln lt nycoeff*nocoeff*2 then begin
				print, 'x_opt2darc: Breaking... number of lines < 2 x fitted parameters.'
				break
			endif
			if tnln lt nycoeff*ngd_ordr*2 then begin
				print, 'x_opt2darc: Breaking... number of lines < 2 x N_ord x nycoeff.'
				break
			endif
			
			prev_tnln = tnln
		endwhile
	
		fitln = fitln[where(outmask)]
		nln = n_elements(fitln)
		
		vrms = sqrt(total(((fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms)^2)/nln)
		print, 'x_opt2darc: '+ string(nln,format='(I4)') + ' lines with delta velocity rms (m/s) of '+string(vrms,format='(F10.1)')
	endif
	
	if keyword_set( VELCUT ) then begin
		; iteratively crop points with delta v > velcut and refit, cross your fingers and hope the thing doesn't run away
		print, 'x_opt2darc: Iterative velocity cut'
		iiter = 0
		prev_nln = 0
		while (1) do begin
			cutmask = abs(fitln.wave_fit-fitln.wave)/fitln.wave_fit*c_ms lt velcut
			x_fit2dlines, fitln, nycoeff, nocoeff, msk=cutmask
			iiter = iiter + 1
			
			fitln = fitln[where(cutmask)]
			nln = n_elements(fitln)
			vrms = sqrt(total(((fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms)^2)/nln)
			print, '   ' + string(iiter,format='(I3)') + ': ' + string(nln,format='(I4)') + ' lines with delta velocity rms (m/s) of '+ string(vrms,format='(F10.1)')
			
			if nln eq prev_nln then begin
				print, 'x_opt2darc: Breaking... no lines rejected in previous iteration.'
				break
			endif
			if nln lt nycoeff*ngd_ordr then begin
				print, 'x_opt2darc: Stopping... number of lines < N_ord x nycoeff.'
				stop
			endif
			
			prev_nln = nln
		endwhile
	
		; fitln = fitln[where(cutmask)]
		; nln = n_elements(fitln)
		
		vrms = sqrt(total(((fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms)^2)/nln)
		print, 'x_opt2darc: '+ string(nln,format='(I4)') + ' lines with delta velocity rms (m/s) of ' + string(vrms,format='(F10.1)')

	endif
	
	; for k=0L, nln-1 do begin
	; 	print, string(fitln[k].wave_fit,format='(F12.6)'), string(fitln[k].wave,format='(F12.6)'), string(fitln[k].order,format='(I4)'), string(fitln[k].pix,format='(F12.3)'), string(fitln[k].wave_fit-fitln[k].wave,format='(F12.5)'), string(fitln[k].pix,format='(F12.3)'), string((fitln[k].wave_fit-fitln[k].wave)/fitln[k].wave_fit*c_ms,format='(F12.5)')
	; endfor
	
	; repeat final fit, retrieve fit parameters and coefficients
	x_fit2dlines, fitln, nycoeff, nocoeff, res=res, nrmt=nrmt, nrmp=nrmp
	
	; Arc 2D Fit OUTPUT Structure
	out_str = { $
		nrm: nrmp, $
		nrmt: nrmt, $
		ny: nycoeff, $
		no: nocoeff, $
		res: res }
	
	; Convenience structure for qa files 
	qa_struct = {wave: dblarr(round(sz[1])), pix: dblarr(round(sz[1])), order: 0, vrms: 0}
	qastats = replicate(qa_struct, ngd_ordr)
	; Print summary and QA files
	if keyword_set( QATXT ) then begin
		print, "x_opt2darc: Writing ", qatxt
	    openw, 1, qatxt
		printf, 1, '#ID  ORD V_RMS  NL    WV_MIN   WV_MAX'
	endif
	if keyword_set( QAWV ) then begin
		print, "x_opt2darc: Writing ", qawv
	    openw, 2, qawv
	endif
	print, 'x_opt2darc: SUMMARY'
	print, 'ID   ORD V_RMS  NL    WV_MIN   WV_MAX'
	for j=0L, ngd_ordr-1 do begin
		i = gd_ordr[j]
		ord = ordr_str[i].order
		ord_l = where(fitln.order eq ord, olcount)
		
		nonzi = where(not aspec[*,i] eq 0.0)
		nzmin = min(nonzi, max = nzmax)
		
		; Find real order bounds
		bound_pix = dindgen(nzmax-nzmin+1.0)+nzmin
		bound_nrm = 2. * (bound_pix - nrmp[0])/nrmp[1]
		bound_npix = n_elements(bound_nrm)
		bound_worky = flegendre(bound_nrm[*], nycoeff)
		
		bound_t_nrm = 2. * (replicate(float(ord), bound_npix) - nrmt[0])/nrmt[1]
		bound_work2d = dblarr(bound_npix, nycoeff*nocoeff)
		bound_workt = flegendre(bound_t_nrm[*], nocoeff)
		for ii=0, nocoeff-1 do begin
			for jj=0, nycoeff-1 do begin
				bound_work2d[*, jj*nocoeff+ii] = bound_worky[*, jj] * bound_workt[*, ii]
			endfor
		endfor
		bound_wv = dblarr(bound_npix)
		bound_wv[*] = bound_work2d # res / ord
		
		wvmn = min(bound_wv, max=wvmx)
		
		qastats[j].order = ord
		qastats[j].wave = bound_wv
		qastats[j].pix = bound_pix
		
		if olcount ne 0 then begin
			ord_vrms = sqrt(total(((fitln[ord_l].wave_fit - fitln[ord_l].wave)/fitln.wave_fit[ord_l]*c_ms)^2)/olcount)
			qastats[j].vrms = ord_vrms
			qtline = string([j, ord, ord_vrms, olcount, wvmn, wvmx], FORMAT='(I2,2X,I3,2X,F5.1,2X,I3,2X,F7.1,2X,F7.1)')
			if keyword_set( QATXT ) then printf, 1, qtline   
			print, qtline
		endif else begin
			qtline = string([j, ord, 0, wvmn, wvmx], FORMAT='(I2,2X,I3,2X,"-----",2X,I3,2X,F7.1,2X,F7.1)')
			if keyword_set( QATXT ) then printf, 1, qtline
			print, qtline
		endelse
		
		if keyword_set( QAWV ) then begin
			printf, 1, ii, bound_pix, bound_wv, format = '((i3)' + string(n_elements(bound_pix)) + '(1x, f8.1)' + string(n_elements(bound_wv)) + '(1x, f11.5))'
		endif
		
	endfor
	
	vrms = sqrt(total(((fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms)^2)/nln)
	qtline = string([vrms, nln], FORMAT='(9X,F5.1,2X,I3,18X)')
	if keyword_set( QATXT ) then begin
		printf, 1, '-------------------------------------'
		printf, 1, qtline
		close, 1
	endif
	print, '-------------------------------------'
	print, qtline
	
	if keyword_set( QAWV ) then close, 2
	
	; QA PS
	if keyword_set( QAFIL ) then begin
		print, "x_opt2darc: Writing ", qafil
		x_psopen, qafil, /maxs
		clr = getcolor(/load)
		
		; Main plot
		!p.multi=[0,1,1]
		qs = qastats.wave ; this is here, like this, because IDL is poorly designed
		qs = qs[where(qs gt 0.0)] ; IGNORE 0 as min!!
		ptmn = min([qs, fitln.wave], max=ptmx)
		; qspix = qastats.pix
		; pxmx = max(qspix)
		; xlbl = pxmx*0.05
		plot, [0.], [0.], color=clr.black, $
			background=clr.white, charsize=1.5, xrange=[ptmn, ptmx], $
			yrange=[0.,sz[1]], xstyle=1, ystyle=1, ytitle='Row', $
			xtitle='Wavelength', xmargin=[11,2], ymargin=[5,2], /nodata
		
		resfct = 4000.0 ; was 500
		for j=0L, ngd_ordr-1 do begin
			ord = qastats[j].order
			qavalid = where(qastats[j].wave gt 0.0)
			qa_wv = qastats[j].wave[qavalid]
			qa_px = qastats[j].pix[qavalid]
			oplot, qa_wv, qa_px, color=clr.blue
			pts = where(fitln.order eq ord, npts)
			if npts ne 0 then begin
				sres = fitln[pts].wave_fit - fitln[pts].wave
				oplot, [fitln[pts].wave_fit - sres*resfct], [fitln[pts].pix], psym=1, color=clr.black
			endif
			; Label
			xyouts, 0.5, 0.96, 'Arc 2D FIT (Res x'+string(resfct, FORMAT='(I4)')+') nx='+strtrim(nocoeff,2)+ $
				' ny='+strtrim(nycoeff,2)+' ' + 'RMS='+string(vrms,format='(f7.2)')+' m/s', $
				color=clr.black, charsize=1.5, /normal, alignment=0.5
		endfor

		; Individual fit plots
		!p.multi = [0, 3, 2]
		for j=0L, ngd_ordr-1 do begin
			ord = qastats[j].order
			qavalid = where(qastats[j].wave gt 0.0, qav_cnt)
			qa_wv = qastats[j].wave[qavalid]
			qa_px = qastats[j].pix[qavalid]
			ptmn = min(qa_wv, max=ptmx)
			plot, qa_px, qa_wv, color=clr.black, $
				background=clr.white, charsize=1.5, yrange=[ptmn, ptmx], $
				xrange=[0.,sz[1]], xstyle=1, ystyle=1, xtitle='Row', $
				ytitle='Wavelength (res x 500)', xmargin=[11,2], ymargin=[5,1], /nodata
			
			; Fit
			oplot, qa_px, qa_wv, color=clr.blue
			
			; Data
			pts = where(fitln.order eq ord, npts)
			pixrms = 9.99
			if npts NE 0 then begin
				sres = fitln[pts].wave_fit - fitln[pts].wave
				oplot, [fitln[pts].pix],  [fitln[pts].wave_fit - sres*500.0], $ ; was *100
					psym=1, color=clr.black
				pixrms = sqrt(total(sres^2)/float(npts))
			endif
			
			; Label
			ylbl = ptmx - (ptmx-ptmn)*0.08*(findgen(5)+1)
			xyouts, sz[1]*0.05, ylbl[0], 'Order = '+strtrim(ord,2), color=clr.black, charsize=1.5
			dwv = abs((qa_wv[0]-qa_wv[qav_cnt-1])/float(qa_px[0]-qa_px[qav_cnt-1]))
			xyouts, sz[1]*0.05, ylbl[2], '!9Dl!X = '+string(dwv,format='(f6.4)'), color=clr.black, charsize=1.5
			xyouts, sz[1]*0.05, ylbl[1], 'RMS(pix) = '+string(pixrms/dwv,format='(f4.2)'), color=clr.black, charsize=1.5
		endfor

		; Residual plots
		!p.multi = [0, 3, 2]
		for j=0L, ngd_ordr-1 do begin
			ord = qastats[j].order
			pts = where(fitln.order eq ord, npts)
			all_res = (fitln.wave_fit - fitln.wave)/fitln.wave_fit*c_ms
			max_res = max(abs(all_res)) ; used for residual plot max and min limits
			qavalid = where(qastats[j].wave gt 0.0)
			qa_wv = qastats[j].wave[qavalid]
			ptmn = min(qa_wv, max=ptmx)
			xlbl = ptmn+5.0
			if npts NE 0 then begin
				plot, [fitln[pts].wave_fit], [(fitln[pts].wave_fit - fitln[pts].wave)/fitln[pts].wave_fit*c_ms], $
					psym=1, charsize=1.8, background=clr.white, color=clr.black, xtitle='Wavelength (A)', ytitle='Residual (m/s)', $
					xmargin=[11,1], ymargin=[5,1], yrange=[-max_res*1.05, max_res*1.05], $
					xrange=[ptmn,ptmx], xstyle=1, ystyle=1
				xyouts, xlbl, max_res*0.9, 'Order = '+string(ord, FORMAT='(i3)')+$
				'   Nlin = '+string(npts, FORMAT='(i3)'), charsize=1.3
				xyouts, xlbl, max_res*0.8, 'vRMS = '+string(qastats[j].vrms, FORMAT='(f7.1)')+' m/s', charsize=1.3
			endif else begin
				plot, [0], [0], psym=1, charsize=1.8, background=clr.white, color=clr.black, $
					xtitle='Wavelength (A)', ytitle='Residual (m/s)', $
					xmargin=[11,1], ymargin=[5,1], yrange=[-max_res*1.05, max_res*1.05], $
					xrange=[ptmn,ptmx], xstyle=1, ystyle=1, /NODATA
				xyouts, xlbl, max_res*0.9, 'Order = '+string(ord, FORMAT='(i3)')+'   Nlin = 0', charsize=1.3
			endelse

			oplot, [-9e9, 9e9], [0., 0.], linestyle=1, color=clr.blue
		endfor
		
		x_psclose
		!p.multi=[0,1,1]
		
		replace_title = '"' + '%%Title: '+qafil + ' ' +systime() + '"'
		ps_replacetitle, replace_title, qafil
		spawn, 'gzip -f '+qafil
	endif
	
	return

end
