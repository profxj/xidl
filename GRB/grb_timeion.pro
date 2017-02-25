;+ 
; NAME:
; grb_timeion
;   Version 1.0
;
; PURPOSE:
;    Calculates the excitation rate as a function of Luminosity
;      and distance and abundance
;
; CALLING SEQUENCE:
;   
; 
;
; INPUTS:
;   grb  -- Structure for GRB afterglow
;
;
; RETURNS:
;   
; OUTPUTS:
;   Creates a Plot
;
; OPTIONAL OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;   
;
;
; PROCEDURES/FUNCTIONS CALLED:
;
; REVISION HISTORY:
;   May-2007 Written by JXP
;-
;------------------------------------------------------------------------------


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
pro grb_timeion, grb, LSTRT=lstrt, OUTFIL=outfil, ALSIHHEN=alsihhen, $
                 ALL=all, nH=nH, NTIME=ntime

  if  N_params() LT 1  then begin 
    print,'Syntax - ' + $
             'nphot = grb_calcphot( dE, dt, grb, Rphot=, /SILENT) [v1.0]'
    return
  endif 

  c = x_constants()

  if not keyword_set(OUTFIL) then outfil = 'timeion.idl'
  if not keyword_set(NTIME) then ntime = 1000L  ;; Time steps
  if not keyword_set(TSTEP) then tstep = 1.  ;; Seconds
;  if not keyword_set(nHe) then nHe = 0.1 
  if not keyword_set(FTAU) then ftau = 1e-4
  if not keyword_set(SIG) then sig = 1d-18
  if not keyword_set(NH) then nh = 1                      ;; H density (cm-3)
  if not keyword_set(t0) then t0 = 10.               ;; starting time (seconds)
;  if not keyword_set(NTSTP) then ntstp = 1e3  ; Number of time steps
  if not keyword_set(LSTRT) then lstrt = 1 * c.pc  ; Starting radius of gas
  if not keyword_set(LEND) then lend = 100 * c.pc  ; Ending radius

  if keyword_set(HHEN) then begin
      ions = [ [1,1], $
               [2,1], $
               [2,2], $
               [7,1], $
               [7,2], $
               [7,3], $
               [7,4], $
               [7,5]  $
             ]
      abnd = [1., 1/12., 1/12., replicate(1e-4,5)]
      neut = [0, 1, 3]
  endif
  if keyword_set(ALSIHHEN) then begin
      ions = [ [1,1], $
               [2,1], $
               [2,2], $
               [7,1], $
               [7,2], $
               [7,3], $
               [7,4], $
               [7,5], $
               [13,2], $
               [13,3], $
               [14,2], $
               [14,3], $
               [14,4]  $
             ]
      abnd = [1., 1/12., 1/12., replicate(1e-4,5), replicate(10.^(-4.5),2), $
             replicate(10.^(-3.5),3)]
      neut = [0, 1, 3,8,10]
  endif
  if keyword_set(ALL) then begin
      ions = [ [1,1], $
               [2,1], $
               [2,2], $
               [7,1], $
               [7,2], $
               [7,3], $
               [7,4], $
               [7,5], $
               [13,2], $
               [13,3], $
               [14,2], $
               [14,3], $
               [14,4], $
               [6,2], $
               [6,3], $
               [6,4], $
               [12,1], $
               [12,2]  $
             ]
      abnd = [1., 1/12., 1/12., replicate(1e-4,5), replicate(10.^(-5.5),2), $
              replicate(10.^(-4.5),3), replicate(10.^(-3.5),3), $
              replicate(10.^(-4.5),2)]
      neut = [0, 1, 3,8,10,13,16]
  endif
  sz_i = size(ions, /dimen)

  ;; Energy
  nengy = 250L  ;; Number of energy steps (linear)
  e0 = 5.
  estep  = (130. - e0) / nengy
  print, 'grb_timeion: Using energy intervals of ', estep, 'eV'
  engy = e0 + findgen(nengy)*estep ;; eV   5 eV to 130 eV
  ion_str = { $
            name: '', $
            Z: 0L, $
            ion: 0L, $
            Eth: 0., $
            sigma: dblarr(nengy), $  ;; cm^-2
            abnd: 0., $
            ftot: fltarr(1000000L) $
            }
  all_ion = replicate(ion_str, sz_i[1])
  all_ion.abnd = abnd
  

  ;; Fill up ion struct
  for qq=0L,sz_i[1]-1 do begin
      ;; Ion
      all_ion[qq].Z = ions[0,qq]
      all_ion[qq].ion = ions[1,qq]
      if ions[1,qq] EQ min(ions[1,where(ions[0,*] EQ ions[0,qq])]) then $
        all_ion[qq].ftot = 1.

      ;; Cross-section
      all_ion[qq].sigma = x_photocross(all_ion[qq].Z, $
                                       all_ion[qq].ion, $
                                       engy, Eth=Eth)
      all_ion[qq].Eth = Eth

      ;; Threshold correction
      b = where(engy LT Eth,nb)
      if nb EQ 0 then stop
      frac = (engy[b[nb-1]] + estep - Eth) / estep
      all_ion[qq].sigma[b[nb-1]] = x_photocross(all_ion[qq].Z, $
                                       all_ion[qq].ion, Eth) * frac
  endfor


  ;; Fill up photons
  tot_phot= dblarr(nengy, ntime)
  dt = fltarr(2,ntime)
  dt[0,*] = t0 + findgen(ntime)*tstep
  dt[1,*] = t0 + (findgen(ntime)+1)*tstep
  for ii=0L,nengy-1 do begin
      dE = [engy[ii],engy[ii]+estep]
      tot_phot[ii,*] = grb_calcphot(dE, dt, grb, /silen)
  endfor
  stop

  ;; Arrays
  dl = (1e16 / nH) < 1e16       ;; cm -- Optically thin demanded
  ncell = round((lend-lstrt) / dl)
  print, 'Ncell = ', ncell
  tau = fltarr(ncell)
  phi_i = fltarr(ncell)
  phi_f = fltarr(ncell)

  ;; Linear radial scheme
  rad = lstrt + dl*dindgen(ncell)

  Ni = dl * nH ; Column of H (cm^-2)

;  ;; Restricting to H for now
;  all_ion = all_ion[0]

  uniZ = all_ion[uniq(all_ion.Z, sort(all_ion.Z))].Z
  nelm = n_elements(uniZ)
  maxi = intarr(nelm)

  for ss=0L,nelm-1 do begin
      maxi[ss] = max(all_ion[where(all_ion.Z EQ uniZ[ss])].ion)
  endfor

  for qq=0L,ncell-1 do begin
      if qq mod 1000 EQ 0 then print, 'qq = ', qq

      ;; Radial drop in Photons
      if qq EQ 0 then photon = tot_phot / (4*!pi*lstrt^2)  $
      else photon = photon * (rad[qq-1]/rad[qq])^2

      ;; Loop on time
      for tt=0L,ntime-1 do begin
          ;; Photoionized?
          if total(all_ion.ftot[qq]) LT 1e-5 then break

          ;; Loop on elements
          for kk=0L,nelm-1 do begin
              ion = where(all_ion.Z EQ uniZ[kk], nion)
              ;; Loop on ions
              for ii=0L,nion-1 do begin
                  idx = ion[ii]
                  if all_ion[idx].ftot[qq] LT 1e-5 then continue
                  tau_phot = all_ion[idx].sigma * photon[*,tt]
                  tot_tau = total(tau_phot)
                  newf = exp(-1.*tot_tau)
                  
                  ;; Attenuate
                  deltaN = all_ion[idx].ftot[qq]*(1 - newf) * Ni * all_ion[idx].abnd
                  photon[*,tt] = (photon[*,tt] - deltaN*tau_phot/tot_tau) > 1.

                  ;; Propogate to next ion
                  if all_ion[idx].ion LT maxi[kk] then $
                    all_ion[ion[ii+1]].ftot[qq] = all_ion[ion[ii+1]].ftot[qq] + $
                      all_ion[idx].ftot[qq]*(1 - newf)
                  ;; Save change
;                  if qq EQ 3000 AND ii EQ 1 then stop
                  all_ion[idx].ftot[qq] = newf*all_ion[idx].ftot[qq]
;                  if kk EQ 1 then stop
              endfor
          endfor
          ;; Neutral?
;          if min(all_ion.ftot[qq]) GT 0.99 then break
      endfor
      if min(all_ion[neut].ftot[qq]) GT 0.90 then break
  endfor
;  x_splot, rad/c.pc, all_ion[0].ftot[0:ncell-1], /bloc
;  x_splot, rad/c.pc, all_ion[1].ftot[0:ncell-1], /bloc
;  x_splot, rad/c.pc, all_ion[2].ftot[0:ncell-1], /bloc
;  stop
  save, rad, all_ion, nH, filename=outfil

  print, 'All done!'
return
end

