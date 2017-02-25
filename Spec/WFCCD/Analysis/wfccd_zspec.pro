;--------------------------------------------------------------------
;+
; NAME:
;   wfccd_zspec
; 
; PURPOSE:
;    neat-o interface for WFCCD quality control
;
; CALLING SEQUENCE:
;    wfccd_zspec,'Extract/Fspec_10.fits.gz', SERENDIP=serendip
;
; INPUTS:
;    final output from wccd_zfind
;    serendip = plot serendips in each slit along with objects
;
; 
; OUTPUTS:
;        plots 1D spectrum + model
;        postage stamp of 2D spectrum +trace       
;
;
; COMMENTS:
;       to do-  scale greyscale images better
;               input just mask number?
;               trace for 2nd image frame
;               qa control (add bitmask flg)
;
; MODIFICATION HISTORY:
;   mg 4/04
;- 
;--------------------------------------------------------------------

PRO wfccd_zspec_event, ev
	COMMON slacker_specs_common, slacer_state,wffspec
 
	widget_control, ev.id, get_uvalue=value
	if (n_elements(value) eq 0) then value = ''
	name = strmid(tag_names(ev, /structure_name), 7, 4)
	
	case (name) of
		'BUTT': handle_button, ev
		'TEXT': handle_text, ev
    		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
  	ENDCASE
END




PRO handle_button, ev
	COMMON slacker_specs_common, slacker_state,wffspec

	WIDGET_CONTROL, ev.ID, GET_UVALUE=uval

	case (uval) of
		'Next': show_object, ev, slacker_state.ix+1
		'Back': show_object, ev, slacker_state.ix-1
		'sky': toggle_twod, ev
		'smooth1d': toggle_smoothing, ev
		'nexp': toggle_nexp, ev
		'DONE': WIDGET_CONTROL, ev.TOP, /DESTROY
	ENDCASE
    END

	
PRO handle_text, ev

	widget_control, ev.id, get_uvalue=uval

	case (uval) of
		'smooth1d': oned_smooth, ev
		'jump_to_slit': jump_to_slit, ev
	end

END

PRO toggle_twod, ev
	COMMON slacker_specs_common, slacker_state,wffspec

	slacker_state.show_twod = not slacker_state.show_twod
	show_twod

    END

PRO toggle_nexp, ev
	COMMON slacker_specs_common, slacker_state,wffspec

	slacker_state.show_nexp = not slacker_state.show_nexp
        read_file2d
	show_twod

END


PRO toggle_oned, ev
	COMMON slacker_specs_common, slacker_state,wffspec

	slacker_state.show_oned = not slacker_state.show_oned
	slacker_state.show_ktrs_img = 0

	show_object, ev, slacker_state.ix
END



PRO show_oned
	COMMON slacker_specs_common, slacker_state,wffspec
	COMMON npk_specs_common, npk_state
	COMMON splot_state,state,graphkeys

	tellam = [6860,7600]

  	 slit = slacker_state.ix
         spec = wffspec[slit].fx

         zans=wffspec[slit].zans
         synth = synthspec(zans, loglam=alog10(wffspec[slit].wave))

      ; set limits not so affected by bad points
        q=where(wffspec[slit].fx ge 0 and wffspec[slit].wave ge 3900)

        medspec = median(spec[q])
        stdspec = stddev(spec[q])
        maxs = medspec + 3.*stdspec 
        mins = medspec - 2.*stdspec
         
      ; smooth synthetic spectrum
	if slacker_state.smooth1d[0] eq 1 then begin
		spec = poly_smooth(spec, slacker_state.smooth1d[1])
		synth= poly_smooth(synth, slacker_state.smooth1d[1])
                str = string('smoothing = ',slacker_state.smooth1d[1]) 
                print,str
        endif

	if (size(state))[0] eq 0 then begin
		tely = 1.5
		splot, wffspec[slit].wave,spec,yrange=[mins,maxs],$
                              xrange=[3500,9000] 
	endif else begin
		if slacker_state.persistent eq 0 then $
		     splot, wffspec[slit].wave,spec,yrange=[mins,maxs], $
                              xrange=[3500,9000] $
		else $
		     splot, wffspec[slit].wave,spec, $
			xrange=[3500,9000], yrange=state.yrange
	endelse
        
      
	if slacker_state.smooth1d[0] eq 1 then begin
	        sxyouts, 8000, maxs - (maxs-mins)*0.1, str,charthick=2
        endif

        
        soplot,wffspec[slit].wave,synth,color=4
        soplot,wffspec[slit].wave,sqrt(wffspec[slit].var),color='red'


	tely=state.yrange[1]
	sxyouts, tellam[0], tely, 'B', color=1,charthick=3
	sxyouts, tellam[1], tely, 'A', color=1,charthick=3
	overplot_lines

END


PRO overplot_lines
	COMMON slacker_specs_common, slacker_state,wffspec

        slit = slacker_state.ix
        n = n_elements(slacker_state.linenames)
        zlines = wffspec[slit].zans.z * slacker_state.linewvl  + $
                   slacker_state.linewvl
        for i = 0, n-1 do begin
                wl = zlines[i]
                wn = slacker_state.linenames[i]

                soplot, [wl, wl], [-500, 5000], color=2,linestyle=1
        endfor

END



PRO show_twod
	COMMON slacker_specs_common, slacker_state,wffspec

      ; SET 2D SPECTRUM Sky/NoSKY, NEXP = 1/2
        if (slacker_state.show_nexp and slacker_state.show_twod) then begin
          twod=slacker_state.spec2d_n1
          wave2d = slacker_state.wave2d_n1
        endif
        if (not slacker_state.show_nexp and slacker_state.show_twod) then begin
          twod=slacker_state.spec2d_n2
          wave2d = slacker_state.wave2d_n2
        endif
        if (slacker_state.show_nexp and not slacker_state.show_twod) then begin
          twod=slacker_state.subsky2d_n1
          wave2d = slacker_state.wave2d_n1
        endif
       if (not slacker_state.show_nexp and not slacker_state.show_twod) then begin
          twod=slacker_state.subsky2d_n2
          wave2d = slacker_state.wave2d_n2
        endif


      ; Allow toggle between raw and sky subtracted frame
        slit = slacker_state.ix
        spec = twod[*,wffspec[slit].ycen-24:wffspec[slit].ycen+24]
        
        xatv,spec,/align

      ; Slight hack to cut off trace at reasonable location  
        yt = ceil(wffspec[slit].trace)

        q=where((yt - wffspec[slit].ycen + 25) gt 0 and $
                  yt - wffspec[slit].ycen + 25 le 50,sz)

      ; plot trace
        x=findgen(sz)
        xatvplot,[x],[yt[q] - wffspec[slit].ycen + 25],thick=0.5,color='green'
 
      ; wavelength array along trace
        wave=fltarr(sz)
        for j=0,sz-1 do begin
            wave[j] = wave2d[j,yt(q(j))]
        endfor

      ; plot location of important lines
        zlines = wffspec[slit].zans.z * slacker_state.linewvl  + $
                                        slacker_state.linewvl
        tmp =findgen(10)
        xline= fltarr(10)
        yline1 = tmp
        yline2 = tmp+40
        for j=0,n_elements(slacker_state.linenames)-1 do begin
            diff = min(abs(wave - zlines[j]),q)
            if (q ne 0) then xatvxyouts,q,50,slacker_state.linenames[j],charsize=2
            xline[0:*] = q
            xatvplot,[xline],[yline1],color='red'
            xatvplot,[xline],[yline2],color='red'
        endfor

        
      ; Add wavelength soln to 2D plot
        wtmp = wave2d[*,wffspec[slit].ycen]
        wlines= (findgen(11)/2. + 3.5)*1000
        wstr = string(ceil(wlines))
        wstr = strcompress(wstr,/remove_all)
        for j=0,10 do begin
            diff = min(abs(wtmp - wlines[j]),q)
            xatvxyouts,q,-15,wstr[j],charsize=2,color='blue',alignment=0.5
            xline[0:*] = q
            xatvplot,[xline],[yline1],color='blue'
        endfor


END



PRO show_object, ev, ix
	COMMON slacker_specs_common, slacker_state, wffspec
	common splot_state, state, graphkeys

	if ix lt 0 or ix gt n_elements(wffspec.slit_id)-1 then begin
		null = dialog_message(['There are not more objects!'])
		return
	endif


	slacker_state.ix = ix
	show_twod
	show_oned
	describe_spectrum, ev


END



PRO show_message, ev, msg
	WIDGET_CONTROL, ev.TOP, GET_UVALUE=textwid
	WIDGET_CONTROL, textwid, SET_VALUE=msg
END

PRO describe_spectrum, ev
	COMMON slacker_specs_common, slacker_state, wffspec
        get_info,msg
	show_message, ev, msg
    END



PRO get_info,msg
	COMMON slacker_specs_common, slacker_state, wffspec
        
        c=2.9979e5

        slit = slacker_state.ix
        id = wffspec[slit].slit_id

        idstr = 'Slit ID ='   + string(wffspec[slit].slit_id,format='(i7)')

      ; NOTE IF OBJECTS IS A SERENDIP
        if (wffspec[slit].obj_id eq 'b') then begin
               idstr = 'Slit ID ='  + string(wffspec[slit].slit_id,'srndp')
               idstr = strcompress(idstr)
        endif

      ; Slit number and frame info
        snum  = strcompress('Slit # = '+string(slit,format='(i5)')+' of '+$
                           string(slacker_state.nslit,format='(i5)'))
        frame = 'Frame = '+slacker_state.frame

      ; Mask info
        mask  = 'Mask ID =  ' + wffspec[slit].field
        pstr  = 'wfccd_pri = '+ string(wffspec[slit].wfccd_pri,format='(i4)')
        magstr= 'rmag    = '  + string(wffspec[slit].rmag,format='(f7.2)')
        zstr  = 'zfind z = '  + string(wffspec[slit].zans.z*c,format='(f10.3)')
        vgroup= 'group z = '  + string(wffspec[slit].vgroup,format='(f9.2)')
        chi2  = 'chi2    = '  + string(wffspec[slit].zans.rchi2,format='(f9.2)')
       
        subc  = wffspec[slit].zans.class +'  '+wffspec[slit].zans.subclass


        x='             '
        line='--------------------------------------------------------------'

        msg = string(format=$
           '(%"%19s%s%s\n%s\n%19s%s%s\n%19s%s%s\n%19s%s%s\n%19s%s%s\n")',$
                     snum,x,frame,$
                     line,$
                     idstr,x,mask,$
                     magstr,x,pstr,$
                     chi2,x,subc,$
                     zstr,x,vgroup)


      ; Add extra line if sdss redshift exist
        if (wffspec[slit].zsdss gt 0) then begin
            
            zsdss = 'SDSS redshift = ' + string(wffspec[slit].zsdss,$
                                                format='(f9.2)')
            b = ''
            if (wffspec[slit].berlind eq 1) then berlind = 'berlind galaxy'
            msg = string(format=$
           '(%"%19s%s%s\n%s\n%19s%s%s\n%19s%s%s\n%19s%s%s\n%19s%s%s\n\n%s%s%s")',$
                     snum,x,frame,$
                     line,$
                     idstr,x,mask,$
                     magstr,x,pstr,$
                     chi2,x,subc,$
                     zstr,x,vgroup,$
                     zsdss,x,berlind)

         endif
                
END

PRO toggle_smoothing, ev
	COMMON slacker_specs_common, slacker_state, wffspec
	
	slacker_state.smooth1d[0] = not slacker_state.smooth1d[0]
	show_object, ev, slacker_state.ix
END


PRO oned_smooth, ev
	COMMON slacker_specs_common, slacker_state

	widget_control, ev.id, get_value=val
	slacker_state.smooth1d[1] = fix(val)
	widget_control, ev.id, set_value = $
		strcompress(string(slacker_state.smooth1d[1]),/remove_all)
        show_oned

    END


PRO jump_to_slit, ev
	COMMON slacker_specs_common, slacker_state

	widget_control, ev.id, get_value=val
	slacker_state.ix = fix(val)
	widget_control, ev.id, set_value = $
		strcompress(string(slacker_state.ix),/remove_all)
        show_object,ev,slacker_state.ix

END

PRO read_file2d
	COMMON slacker_specs_common, slacker_state,wffspec

        if (slacker_state.show_nexp) then begin
            nexp=1
            if not (slacker_state.read_n1) then begin
                file2d = wffspec[0].spec2d_fil[nexp]        
                slacker_state.spec2d_n1=xmrdfits(file2d,1)
                slacker_state.wave2d_n1=xmrdfits(file2d,2)
                slacker_state.subsky2d_n1=xmrdfits(file2d,3)
                slacker_state.read_n1 = 1
            endif
         endif

        if not (slacker_state.show_nexp) then begin
            nexp=2
            if not (slacker_state.read_n2) then begin
                file2d = wffspec[0].spec2d_fil[nexp]        
                slacker_state.spec2d_n2=xmrdfits(file2d,1)
                slacker_state.wave2d_n2=xmrdfits(file2d,2)
                slacker_state.subsky2d_n2=xmrdfits(file2d,3)
                slacker_state.read_n2 = 1
            endif
         endif

END

Pro wfccd_zspec,fspec_fil,SERENDIP=serendip

	COMMON slacker_specs_common, slacker_state,wffspec

        !p.font =6
        c=2.9979e5


      ; Important lines, convert these to vacuum wvls
        lines = [3726.16,3728.9,$   
                 3933.7,3968.5,$
                 4861.32,4958.9, 5006.84,$
                 6548.1,6562.80, 6583.6]
        airtovac,lines


        nlines= ['[OII]',' ',$
                 ' ','CaH/K',$
                 'H!7b!6',' ',' ',$
                 ' ','H!7a!6',' ']

	slacker_state = {$
		ix: 0, $		; the object # being examined
                nslit:0.0,$             ; number of good slits
                frame: ' ',$            ; frame          
                smooth1d: [1,5], $	; [smooth y/n, smoothing in pixels]
		persistent: 0,$		; Are splot options persistent
		show_twod: 1, $		; Show the twod ATV plot
		show_nexp: 0, $		; Show the nexp 
		read_n1: 0, $		; has nexp1 been read
		read_n2: 0, $		; has nexp2 been read
		linewvl: lines,$	;  list of names to display lines for
		linenames: nlines,$	;  list of names to display lines for
                spec2d_n1:fltarr(2048,2047),$
                wave2d_n1:fltarr(2048,2047),$
                subsky2d_n1:fltarr(2048,2047),$
                spec2d_n2:fltarr(2048,2047),$
                wave2d_n2:fltarr(2048,2047),$
                subsky2d_n2:fltarr(2048,2047)$
	}


      ; READ FINAL FSPEC, EXAMINE ONLY PRIMARY OBJECTS
        wfccd_wrfspec, wffspec, fspec_fil, /read
        gslit = where(wffspec.flg_anly NE 0 and wffspec.obj_id eq 'a',ngd) 

      ; ADD SERENDIPS
        if (keyword_set(serendip)) then gslit = where(wffspec.flg_anly NE 0,ngd) 

        wffspec=wffspec[gslit]
        slacker_state.nslit = ngd
        slacker_state.frame = fspec_fil
        

      ; PRINT INITIAL MESSAGE ON MASK
        q = where(abs(c*wffspec.zans.z - wffspec.vgroup) le 1500,ngrp)
        print,'Number of slits within 1500km/s of group = ',string(ngrp)
        q = where(wffspec.berlind eq 1,nbld)
        print,'Number of berlind galaxies',nbld
        q = where(wffspec.obj_id eq 'b',nser)
        print,'Number of serendip objects',nser 


      ; READ 2D SPECTRUM FRAME
        read_file2d

      ; get info on first slit
        get_info,msg
        

	; HIGHEST LEVEL
  	base = WIDGET_BASE(/column)
	WIDGET_CONTROL, /MANAGED, base

	; 	Level 1

	b1 = WIDGET_BASE(base, /frame, /row)
  	button = WIDGET_BUTTON(b1, VALUE='Back', UVALUE='Back')
  	button = WIDGET_BUTTON(b1, VALUE='Next', UVALUE='Next')
  	text = WIDGET_LABEL(b1, value='jump to slit =')
  	text = WIDGET_TEXT(b1, XSIZE=5,value=string(slacker_state.ix),$
                           /EDITABLE,uvalue='jump_to_slit')

	b2 = WIDGET_BASE(base, /frame, /row)


	;		Level 2 Middle
	l2_right = WIDGET_BASE(b2, /column, /frame)
	l2_right_top = WIDGET_BASE(l2_right,/row)
	button = WIDGET_BUTTON(l2_right_top, VALUE='Smooth1d', uvalue='smooth1d')
	text2 = WIDGET_TEXT(l2_right_top, XSIZE=5, value='5', $
		/EDITABLE, uvalue='smooth1d')

	l2_right_mid = WIDGET_BASE(l2_right, /row)
	button = WIDGET_BUTTON(l2_right_mid, value='Nexp', uvalue='nexp')


	l2_right_bot = WIDGET_BASE(l2_right, /row)
	button = WIDGET_BUTTON(l2_right_bot, value='Sky/noSky', uvalue='sky')



	;	Level 3
  	text = WIDGET_TEXT(b2, XSIZE=80,value=msg)
  	button4 = WIDGET_BUTTON(base, VALUE='Done', UVALUE='DONE')

        show_twod
        show_oned
        
  	WIDGET_CONTROL, base, SET_UVALUE=text
  	WIDGET_CONTROL, base, /REALIZE
  	XMANAGER, 'wfccd_zspec', base



    end






