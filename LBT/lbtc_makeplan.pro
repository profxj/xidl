
;procedure that prepare a plan file for lris_ccdproc

;logfile  as it comes out of the telescope 
;PATH--> path to a folder where images are (Defaut ./Raw/)
;LOG--> name of the log


PRO lbtc_makeplan,  PATH=path, LOG=logfile

  
  if ~keyword_set(PATH) then path='./Raw/'
  
  
;get the list of images (exclude focus)
  cd, path, current=up_dir
  splog, 'Entering ', path

  spawn, "ls lbcb.*fits", blu_list
  spawn, "ls lbcr.*fits", red_list


  num_red=N_ELEMENTS(red_list)
  num_blu=N_ELEMENTS(blu_list)

  splog, 'Found ', num_red, ' red images.'
  splog, 'Found ', num_blu, ' blue images.'
  splog, 'Getting header information'


;get the date
  date=strtrim(STRMID(blu_list[0],7,6),2)
  if ~keyword_set(LOG) then logfile='obs_'+date+'.log'
  
  
  openw, lun, logfile, /get_lun


  for i=0, num_blu-1 do begin
     header=headfits(blu_list[i],exten=0) 
     filt=sxpar(header,"FILTER")
     time=sxpar(header,"EXPTIME")
     obj=sxpar(header,"OBJNAME",/silent)
     type=sxpar(header,"IMAGETYP")
         
     printf, lun, blu_list[i], filt, time, obj, type, 'B', format='(A30,"   ",A20,"  ",F6.1,"   ",A20,"   ",A3,"   ", "   ",A1)'     
  endfor


  for i=0, num_red-1 do begin
     header=headfits(red_list[i],exten=0) 
     filt=sxpar(header,"FILTER")
     time=sxpar(header,"EXPTIME")
     obj=sxpar(header,"OBJNAME",/silent)
     type=sxpar(header,"IMAGETYP")
     
     printf, lun, red_list[i], filt, time, obj, type, 'R', format='(A30,"   ",A20,"  ",F5.1,"   ",A20,"   ",A3,"   ","   ", A1)'

  endfor

  free_lun, lun
  spawn, 'mv '+logfile+' '+up_dir
  
  splog, 'Entering ', up_dir
  cd, up_dir
  splog, 'Created log ', logfile


  

END
