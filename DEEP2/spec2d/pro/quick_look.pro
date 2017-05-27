; script to call deimos_quicklook


; DEIMOS ID test 2001-Aug-22
; DEIMOS readiness review plots 2001-Sep-18
; D. Finkbeiner


; plate scales
; 1- 8.56  5- 8.56
; 2- 8.47  6- 8.47
; 3- 8.47  7- 8.47
; 4- 8.56  8- 8.56


  device, pseudo=8

  chiplist = 1 ;could be a list of chips to examine
  verbset, 4
  ybin = 8
  flat = 1
  anamorph = 1.6 ; anamorphic factor, which is a function of grating and tilt
  spline =  1 ; do spline fit
  maskno = 3110

; -------- Get environment variables, set paths
  deimos_data = getenv('DEIMOS_DATA')+'/'
  if deimos_data eq '/' then message, 'You need to set $DEIMOS_DATA!'

  deep_dir = getenv('DEEP_DIR')+'/'
  if deep_dir eq '/' then message, 'You need to set $DEEP_DIR!'

  maskstr = string(maskno, format='(I4.4)')
  maskdir = deimos_data+maskstr+'/'


  flatname = maskdir+'davis.1000.fits' ; flat
  arcname = maskdir+'davis.1001.fits'

; Read lamp file into a structure
;  lampfilename = '/deep0/marc/deep/spectral/lamplist.dat'
  lampfilename = deep_dir+'spec2d/etc/lamplist.dat'
  lamps = read_lampfile(lampfilename)


; perhaps things like this should be versioned???
  badmaskname = '/deep0/marc/deep/spectral/deimos_badmask.fits.Z'

; mask definition file
  mdfile = '/deep0/marc/deep/targetselection/masks/archive/' + $
    'testmask.3110.fits'
  maskdef = mrdfits(mdfile, 2, /silent)

; objectcat stuff output by Sybase?
  maskfitsfile = maskdir+'masl3110.fits'
  objectcat = deimos_tables(maskfitsfile, bluslits=slitcoords) 
; check for bogus slits
  good = (slitcoords.slitlength ge 0.1) and (slitcoords.slitwidth ge 0.1)
  w = where(good, ngood)
  if n_elements(good) NE ngood then begin 
     slitcoords = slitcoords[w]
     print, 'Removing bad entry in slitcoords!'
     print, where(good eq 0)
  endif 
  
  for i=0, n_elements(chiplist) do begin 
     chipno = chiplist[i]
     print
     print, 'Chipno: ', chipno
     

; -------- read subimage
     flatimage = deimos_read(flatname, chipno) -1600 ;remove overscan, roughly
     deimos_badchip = mrdfits(badmaskname, chipno, /silent) 
;get mask file of bad regions
     
     ivar = (flatimage LT 60000)*(deimos_badchip eq 0)
     
; read in mask file-- hardwired for the moment.  Note that the FITS
;                     table information will in future be included in
;                     data file. This routine does not yet handle
;                     multiple objects/slitlet.

; slit tops:
; note that deimos_slit_id should actually refer to objectcat and
; bluslits, not the testmask.xxxx.fits file used to construct the
; mask! TBD !!
     
     arcimage = deimos_read(arcname, chipno, header=arc_header) -1000
     arcivar = float(arcimage*0.+1./(abs(arcimage)+1))*(deimos_badchip eq 0)
;give saturated points nonzero wt.

; what is pixmask?
     
     specimage = arcimage
     specivar = arcivar

     slitlet_list = 101
     
     deimos_quicklook, chipno, maskno, arc_header, $
       flatimage, flativar, arcimage, arcivar $
       slitcoords, lamps, slitlet_list, outdir=maskdir, $
       ybin=ybin, anamorph=anamorph, /flat, spline=spline, /plot, $
       pspath= './'
  endfor

; verbosity levels

; 0 nothing
; 1  per mask
; 2  per chip
; 3  per slit
; 4  everything!
end
