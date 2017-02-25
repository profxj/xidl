;+
; NAME:
;    long_match_objs
;
; PURPOSE:
;   Creates a postscript file that shows which objects,
;   in a set of reduced LRIS longslit exposures, correspond to which 
;   other objects.
;
; CALLING SEQUENCE:
;    long_match_objs, infile, OUTFILE=OUTFILE, WAVE=WAVE
;
; INPUTS:
;     infile  - String array specifying the file names of the reduced exposures.
;               These files should have been created by the long_reduce code.
;
; OPTIONAL INPUTS:
;      WAVE   - Wavelength at which the trace position should be considered.
;               If not specified, default is to take the mean x-position along
;               the entire trace.  If the exposures include red and blue data,
;               WAVE should be an overlapping wavelength.
;
;      OUTFILE - Name of the outfile to be written.  If not specified, the postscript
;                will be named match_objs.ps.
;
;      CLOSECRIT - Minimum separation (in arcseconds) two objects  must have to be
;                  considered different objects.  Defaults to 1 arcsec.
;
; COMMENTS:
;     This code allows for the possibility of different binnings on red and blue sides.
;     This code does NOT take into account the chip gap on the blue side!
;
; EXAMPLES:
;     long_match_objs, ['sci-lblue0050.fits.gz', '../../LongslitRed_YesFlex/Science/sci-lred0042.fits.gz'], wave=7000.
;
;
; BUGS:
;     This code assumes that the brightest object on the blue side is the same as the
;     brightest object on the red side.  If that is not the case, the code will produce
;     a plot that looks and is wrong.
;
;
; PROCEDURES CALLED:
;
;
; REVISION HISTORY:
;   1-October-2007   Written by L. K. Pollack
;-

PRO long_match_objs, infile, OUTFILE=OUTFILE, WAVE=WAVE, CLOSECRIT=CLOSECRIT

if not(keyword_set(OUTFILE)) then OUTFILE='match_objs.ps'

;find the total number of objects extracted in all exposures.
num_objs=0
for i=0,n_elements(infile)-1 do begin
   str=mrdfits(infile[i],5)
   num_objs1=n_elements(str)
   num_objs=num_objs+num_objs1
endfor

anon = {xmean:0.0, xshift:0.0, x_arcsec:0.0, peakflux:0.0, file:'', filenum:0, platescale:0.0, objnum:0, mean_obj_pos:0.0}
anon = replicate(anon, num_objs)
filenum = findgen(n_elements(infile))

;read the data into a new structure.
counter_obj=0
for i=0, n_elements(infile)-1 do begin
   str=mrdfits(infile[i],5)
   for j=0,n_elements(str)-1 do begin
      anon[counter_obj].file = infile[i]
      anon[counter_obj].filenum = filenum[i]
      anon[counter_obj].xmean = mean(str[j].xpos)

      ;Use xpos near a specific wavelength.
      if keyword_set(WAVE) then begin
         wave_here=where(str[j].wave_opt ge WAVE - 10. and str[j].wave_opt le WAVE + 10., count_wave)
         if count_wave gt 0 then anon[counter_obj].xmean = mean(str[j].xpos[wave_here]) else $
            print, 'No data was found in specified wavelength region.  Resorting to mean xpos.'
      endif

      anon[counter_obj].peakflux = str[j].peakflux
      if strpos(infile[i], 'red') ne -1 then anon[counter_obj].platescale = 0.21*str[j].binning[0] $
         else anon[counter_obj].platescale = 0.135*str[j].binning[0]

      ;BUG: Currently doesn't account for blue chip gap.
  
      counter_obj=counter_obj+1
      
   endfor
endfor


;find the total number of red exposures, and trim down filenames for plotting.
;num_redfile=0
file_plotname=strarr(n_elements(infile))
for i=0,n_elements(infile)-1 do begin
   file_plotname_helper2=strpos(infile[i],'-', /reverse_search)
   file_plotname[i]=strmid(infile[i], file_plotname_helper2+1)
   file_plotname_helper=strpos(file_plotname[i],'.')
   file_plotname[i]=strmid(file_plotname[i], 0, file_plotname_helper)
endfor


max_peakflux_pos=fltarr(n_elements(infile))
;Find the position corresponding to maximum peakflux for each exposure,
;and shift to reference position.
for i=0, n_elements(infile)-1 do begin
   thisfile=where(anon.file eq infile[i])
   where_maxpeak=where(anon[thisfile].peakflux eq max(anon[thisfile].peakflux))
   where_maxpeak=thisfile[where_maxpeak]
   max_peakflux_pos[i]=where_maxpeak

   ;Let the first file define the reference position.
   if i eq 0 then begin
      xmean_ref = anon[where_maxpeak].xmean
      ;Don't shift reference file
      anon[thisfile].xshift = anon[thisfile].xmean
      anon[thisfile].x_arcsec = (anon[thisfile].xshift - xmean_ref) * anon[thisfile].platescale
   endif else begin
      ;Shift all x-positions (and arcsec) to match reference file.
      shifter = anon[where_maxpeak].xmean - xmean_ref
      anon[thisfile].xshift = anon[thisfile].xmean - shifter
      anon[thisfile].x_arcsec = (anon[thisfile].xshift - xmean_ref) * anon[thisfile].platescale
   endelse
endfor

;Decide how many objects are identified using closeness criterion.
;Default: Objects within 1.0 arcsec of each other are defined to be the same.
if not(keyword_set(CLOSECRIT)) then closecrit=1.0
total_objs=1
last_objnum=1
sorted=sort(anon.x_arcsec)
for i=0, n_elements(anon)-2 do begin
   pos1=anon[sorted[i]].x_arcsec
   pos2=anon[sorted[i+1]].x_arcsec
   if abs(pos2 - pos1) gt closecrit then begin
      anon[sorted[i]].objnum = last_objnum
      anon[sorted[i+1]].objnum = last_objnum + 1
      total_objs=total_objs+1
      last_objnum=last_objnum+1
   endif else begin
      anon[sorted[i]].objnum = last_objnum
      anon[sorted[i+1]].objnum = last_objnum
   endelse
endfor

;Calculate the mean object position for each object.
for i=1,total_objs do begin
   here=where(anon.objnum eq i, count_here)
   mean_obj_pos = mean(anon[here].x_arcsec)
   anon[here].mean_obj_pos = replicate(mean_obj_pos, count_here)
endfor


;plot the data to the screen and to a postscript file.
!p.multi=[0,1,3]
!y.minor=-1

plot, anon.xmean, anon.filenum, psym=4, charsize=1.5, $
      yrange=[-1,n_elements(infile)], xrange=[min(anon.xmean)-50., max(anon.xmean)+50.], $
      ytickname=[' ', file_plotname, ' '], ytickinterval=1, xtitle='Pixel Positions of Traces'

plot, anon.x_arcsec, anon.filenum, psym=4, charsize=1.5, $
      yrange=[-1,n_elements(infile)], xrange=[min(anon.x_arcsec)-5., max(anon.x_arcsec)+5.], $
      ytickname=[' ', file_plotname, ' '], ytickinterval=1, xtitle='Relative Arcsec'

plot, anon.x_arcsec, anon.filenum, psym=4, charsize=1.5, $
      yrange=[-1,n_elements(infile)], xrange=[min(anon.x_arcsec)-5., max(anon.x_arcsec)+5.], $
      ytickname=[' ', file_plotname, ' '], ytickinterval=1, xtitle='Relative Arcsec', /nodata
for i=1, total_objs do begin
   here=where(anon.objnum eq i)
   xyouts, anon[here[0]].mean_obj_pos, [0], strcompress(string(anon[here[0]].mean_obj_pos), /rem), orientation=90
endfor


set_plot,'ps'
device, filename=OUTFILE,landscape=0, xsize = 7.5, $
        ysize = 6.5, /inches, xoffset = 0.5, yoffset = 0.5, /color
plot, anon.xmean, anon.filenum, psym=4, charsize=1.5, $
      yrange=[-1,n_elements(infile)], xrange=[min(anon.xmean)-50., max(anon.xmean)+50.], $
      ytickname=[' ', file_plotname, ' '], ytickinterval=1, xtitle='Pixel Positions of Traces'
plot, anon.x_arcsec, anon.filenum, psym=4, charsize=1.5, $
      yrange=[-1,n_elements(infile)], xrange=[min(anon.x_arcsec)-5., max(anon.x_arcsec)+5.], $
      ytickname=[' ', file_plotname, ' '], ytickinterval=1, xtitle='Relative Arcsec'
plot, anon.x_arcsec, anon.filenum, psym=4, charsize=1.5, $
      yrange=[-1,n_elements(infile)], xrange=[min(anon.x_arcsec)-5., max(anon.x_arcsec)+5.], $
      ytickname=[' ', file_plotname, ' '], ytickinterval=1, xtitle='Relative Arcsec', /nodata
for i=1, total_objs do begin
   here=where(anon.objnum eq i)
   xyouts, anon[here[0]].mean_obj_pos, [0], strcompress(string(anon[here[0]].mean_obj_pos), /rem), orientation=90
endfor
!p.multi=0
device,/close
set_plot,'x'


END
