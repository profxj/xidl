pro mike_mk_allfig, ALL=all, ARCS=arcs


  ;; Detector
  rab_suboscan, /debug

  ;; Flats
  flatfig
  fig_pca
  fig_fittrcarc
  fig_trcflat

  ;; Slit
  ;slitprofile
  order_width

  ;; Arcs
  fig_1darcfit
  fig_2darcfit
  fig_fittrcarc
  fig_wavedisp
  fig_wavedisp, /red

  multisky
  
  ;; Extraction
  fwhm_and_atmos
  fig_trcobj
  fig_profile

  ;; Fluxing
  fig_eg21
  fig_eg21, /reduce
  fig_eg21, /feige
  fig_eg21, /feige, /reduce
  ;figflux_1
  stdfig_1_b
  stdfig_2_b

  ;; Coadd
  fig_coadd
  

  return
end

  
