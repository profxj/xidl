;+
; NAME:
;   deimos_bias_problem
;
; PURPOSE:
;   determine level of DEIMOS "bias problem"
;
; CALLING SEQUENCE:
;   deimos_bias_problem, flist
;
; INPUTS:
;
; OPTIONAL INPUTS:
;   flist   - list of files to examine (default is current dir)
;	
; KEYWORDS:
;   chip    - chip number (default 2)
;   all     - do all 8 chips
;   verbose - print output for every file
;
; OUTPUTS:
;   text to screen
;
; RESTRICTIONS:
;
; EXAMPLES:
;
; COMMENTS:
;   Chip 2 appears to have an intermittent fluctuating bias problem.  This
;    code examines the overscan region looking for the symptoms. 
;   For dual amp reads, we only read amp B(?)  The signal appears to
;    always be symmetric in amps A and B for dual mode. 
;
; REVISION HISTORY:
;   2002-Oct-03  Written - DPF
;
;----------------------------------------------------------------------
pro deimos_bias_problem, flist, chip=chip, all=all, verbose=verbose, $
          bias=bias, last=last

  delvarx, bias
; -------- File list
  if NOT keyword_set(flist) then begin 
     flist = findfile('d*.fits', count=ct)
	if ct lt 5 then flist = findfile('d*.fits.gz', count=ct)
        if ct lt 5 then flist = findfile('d*.fits', count=ct)
        print, ct , ' Files found'
  endif else begin
     ct = n_elements(flist) 
  endelse

  if keyword_set(last) then begin 
     flist = flist[ct-last:ct-1]
     ct = n_elements(flist) 
     print, ' Checking only last', last, ' files'
  endif 
  
; -------- Call self if keyword_set(all)
  if keyword_set(all) then begin 
     for chip=1, 8 do deimos_bias_problem, flist, chip=chip, verbose=verbose
     return
  endif 

; -------- Chip number
  if NOT keyword_set(chip) then chip = 2 ; chip 2 is usually the bad one
  print, 'Chip: ', chip

; -------- Begin loop over files
  nbad = 0
  for i=0, ct-1 do begin 
     phdr = headfits(flist[i])
     ampmode = sxpar(phdr, 'AMPMODE')
     single = strmid(ampmode, 0, 6) EQ 'SINGLE' ; is it single-amp mode?

     badboy = single ? chip : chip*2-1   ; use HDU 3 for dual amp mode
     if single OR ((single EQ 0) AND chip LE 4) then begin 
        im = mrdfits(flist[i], badboy, header, /silent)


; -------- Extract datasec from header
        datasec = sxpar(header, 'DATASEC')
        synopsis = sxpar(phdr, 'SYNOPSIS')
        
; -------- Parse datasec using a (clumsy) regexp
        x = (stregex(datasec, '\[([0-9]+):([0-9]+),([0-9]+):([0-9]+)\]', /subexp, /extract))[1:4] - 1
        
; -------- Pull out science image and overscan
        over  = im[x[1]+1:*, *]+32768
        ncol = (size(over))[1]
        
        diff = total(over-shift(over,0,1),1)/ncol
        tover = total(over, 1)/ncol
        diff = diff[20:*]    ; trim first 20 rows
        ntoss = 1               ; toss ntoss highest values
        sind = sort(abs(diff))
        ind = sind[0:n_elements(sind)-ntoss*2-1] ; times 2 because it is a shift difference
        bad = max(abs(diff[ind])) GT 3
        badstr = bad ? '  BAD!' : '      '
        mode = single ? '  SING' : '  DUAL'
        if bad OR keyword_set(verbose) then begin 
           print, flist[i], mode+'   Max:', max(abs(diff[ind])), $
             '    Sig:',stdev(diff[ind]), sxpar(phdr, 'ROTATVAL'), badstr, $
             synopsis, format='(A,A,F7.2,A,F7.3,F7.1,A,"   ",A)'
        endif 
        if bad then nbad = nbad+1
        
     bias = keyword_set(bias) ? [bias, tover] : tover
     endif 
  endfor 
     
  if nbad eq 0 then print, 'All files look fine'
end




