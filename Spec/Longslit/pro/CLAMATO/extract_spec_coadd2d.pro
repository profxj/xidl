;+
;
; Extract 1D spectra from 2D co-added image, and flux them (using a
; sensitivity file from the local directory, which should be one level
; below the Science/ sub-directory
;
; This ingests a list of object names created by GEN_COADD2D_INPUT,
; which had been manually edited to select the desired traces.
;
; INPUTS:
;   maskname   - String with name of slitmask. This is used to search
;                for the input file (by default) and as file suffix
;
; KEYWORD INPUTS:
;   listfil    - Filename for input source list, in columns of slit
;                number, obj name, obj mag, & obj id (starting from 1)
;   scifil     - Explicitly define the Science file to open in the
;                'Science/' directory. Otherwise, will look for
;                'sci-coadd2d_' file
;   sensfunc   - Flux sensitivity file
;   sci_dir    - Science directory, where the LONG_REDUCE sci*.fits.gz
;                file will be looked for. Defaults to Science/
;   final_dir  - Output directory, dafaults to Final2d/
;   plotspec   - If set, will open up each spectrum for inspection
;
; EXAMPLE:
;   extract_spec_coadd2d, maskname, listfil=listfil, scifil=scifil, $
;                         sci_dir=sci_dir, final_dir=final_dir, $
;                         plotspec=plotspec, sensfunc=sensfuncfil
;-

pro extract_spec_coadd2d, maskname, listfil=listfil, $
                          sci_dir=sci_dir, final_dir=final_dir, $
                          plotspec=plotspec, scifil=scifil, $
                          sensfunc=sensfuncfile

if not keyword_set(listfil) then $
     listfil=findfile('coadd2d_input*'+maskname+'*txt')
print, 'Reading list of spectra from ', listfil
  
if not keyword_set(sensfuncfile) then  $
   sensfuncfile=findfile('*sens.fits')

print, 'Fluxing with ', sensfuncfile

if strlen(sensfuncfile) EQ 0 then begin
   print, 'No fluxing sensitivity file found'
   stop
endif

if not keyword_set(sci_dir) then sci_dir ='./Science/'
if not keyword_set(final_dir) then final_dir = './Final2d/'

if not file_test(final_dir, /dir) then file_mkdir, final_dir

if not keyword_set(scifil) then $
   scifil = file_search(sci_dir, 'sci-coadd_2d*.fits*') else begin
   ;; If science file is explicitly specified, check whether it
   ;; includes explicit directory path, otherwise prepend ./Science/
   if not file_test(scifil) then begin
      if file_test(sci_dir+scifil) then scifil=sci_dir+scifil else begin
         print, 'Error: Science file not found'
         stop
      endelse
   endif
endelse

if n_elements(scifil) GT 1 then begin
   scifil = scifil[sort(strlen(scifil))]
   ;; Assume shortest file name is the latest one
   print, 'More than one science file. Using '+scifil[0]
   scifil = scifil[0]
endif
infil = scifil

;; Read input list. We still want this even for 2D coadd in order to get
;; the object names
readcol, listfil, slitnum, name, mag, id1, $
         f='l,a,f,l', stringskip='#'

nspec = n_elements(name)

for ii=0, nspec-1 do begin
   print, 'Co-adding object '+name[ii]
   
   outfil_N=final_dir+maskname+'_'+name[ii]+'_N.fits'
   outfil_F=final_dir+maskname+'_'+name[ii]+'_F.fits'

   objid_tmp = [id1[ii]]
   infil_tmp = infil

   neg_id= where(objid_tmp LT 0, n_neg)
   if n_neg GT 0 then remove, neg_id, objid_tmp, infil_tmp

   long_coadd, infil_tmp, objid_tmp, outfil=outfil_N
   long_fluxcal,outfil_N,sensfuncfile=sensfuncfile,outfil=outfil_F 
   if keyword_set(plotspec) then $
      x_specplot, outfil_F, inflg=2, xsize=1000, ysize=700, /block

endfor

;stop
end
