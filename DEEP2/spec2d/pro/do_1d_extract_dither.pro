pro do_1d_extract_dither, planfile

  if n_params() lt 3 then begin
    planfile = findfile('*.plan')
    planfile = planfile[0]
  endif

  if NOT keyword_set(planfile) then message, 'You must specify a plan file!'

  ;read the plan file to get information on the number of exposures
  read_planfile, planfile, maskname, rawdatadir, outdatadir, flatnames, arcnames, sciencenames

  objposfile= maskname + '.object_positions.txt'

  get_object_positions,maskname,objposfile
  do_extract_dither,objposfile

end