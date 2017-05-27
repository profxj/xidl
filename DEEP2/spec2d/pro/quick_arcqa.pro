
;+
; NAME:
;  quick_arcqa
;
;
; PURPOSE:
;  perform QA tests on calibration files
;
;
; CATEGORY:
;
;   quicklook
;
; CALLING SEQUENCE:
;
; quick_arcqa
;
; RESTRICTIONS:
;
;  must be called in directory containing arcqa files
;
; MODIFICATION HISTORY:
;   2003may jan
;-



pro quick_arcqa

buckledmask = findfile('buckle*.txt',  count=count)

if count gt 0 then begin
      spawn, 'play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/STTNG-redalert.au'
      spawn, 'play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/STTNG-redalert.au'
     cd,current=cwd
     stringsep = strsplit(strcompress(cwd, /REMOVE), '/', /extract)
     mask = stringsep[n_elements(stringsep)-1]
;            string =  +', mask: '+mask+'- Retake arcs!!!!!'
            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2, 'Warning! Slit tracing had problems - '
            printf, 2, 'mask was likely buckled! '
            printf,2
            printf, 2, 'Mask: '+mask
            printf, 2
            printf, 2, 'Inspect the flats;'
            printf,2,'consider retaking them and rechecking for buckling'
            printf,2,'or else do not observe this mask!'
            printf,2
            printf,2,'The typical signature of buckling will be out-of-focus'
            printf,2,'slits on one end of the mask, in-focus on the other.'
            printf,2
            printf,2,'Quicklook tests will not be performed for this mask.'
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning'
      spawn, spawnstring
      return
      endif

qafiles = findfile('arc*qa*.fits',  count=count)

lamps = mrdfits(qafiles[0], 2)

for i=0, count-1 do begin   
   arcf=mrdfits(qafiles[i],0) 
   if i eq 0 then arc = arcf else arc = [[arc], [arcf]] 
endfor

lamps_on = (lamps.element)[uniq(lamps.element, sort(lamps.element))]

wh=where(total(arc,2) NE 0,ngood)
arc=arc[wh,*]

lamps=lamps[wh,*]

; test code
;lamps_on = [lamps_on, 'Zn']

nlamps = n_elements(lamps_on)

for i=0, nlamps-1 do begin
   nlines = long(total(lamps.element eq lamps_on[i]))
   print, lamps_on[i], nlines
   if nlines eq 0 then begin
      spawn, 'play_sound -v 99 -host hamoa -account deepteam /home/deepteam/sounds/doh2.au'
      spawn, 'play_sound -v 99 -host pohue -account deepteam /home/deepteam/sounds/doh2.au'
     cd,current=cwd
     stringsep = strsplit(strcompress(cwd, /REMOVE), '/', /extract)
     mask = stringsep[n_elements(stringsep)-1]
;            string =  +', mask: '+mask+'- Retake arcs!!!!!'
            openw, 2, '/home/deepteam/temp/error.txt'
            printf, 2, 'Warning! No arclines for species '+lamps_on[i]
            printf, 2, 'Mask: '+mask
            printf, 2
            printf, 2, 'If this is either Kr, Ar, Ne, or Xe,'
            printf, 2, 'you should retake this arc!'
            close, 2
      spawnstring = 'cat /home/deepteam/temp/error.txt | tkmessage -type warning'
      spawn, spawnstring
      endif

endfor

return
end
