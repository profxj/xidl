pro parse_DLA, strctnm, fil
  strctnm = {sDLA, dla_qso: '',$
				qso_ra: '',$
				qso_dec: '',$
				qso_zem: 0.0,$
				flg_QSOmag: 0,$
				qso_mag: 0.0,$
				dla_zabs: 0.0,$
				dla_NHI: 0.0,$
				dla_sigNHI: 0.0,$
				dla_abndfil: '',$
				dla_flgFe: 0,$
				dla_FeH: 0.0,$
				dla_sigFeH: 0.0,$
				dla_flgZn: 0,$
				dla_ZnH: 0.0,$
				dla_sigZnH: 0.0 }
  dumc = ''
  dumr = 0.0
  dumi = 0
  openr, 1, fil
	readf, 1,  format='(a15)', dumc
	strctnm.dla_qso = dumc
	readf, 1,  format='(a15)', dumc
	strctnm.qso_ra = dumc
	readf, 1,  format='(a15)', dumc
	strctnm.qso_dec = dumc
	readf, 1,  format='(f9.6)', dumr
	strctnm.qso_zem = dumr
	readf, 1,  format='(i2)', dumi
	strctnm.flg_QSOmag = dumi
	readf, 1,  format='(f9.6)', dumr
	strctnm.qso_mag = dumr
	readf, 1,  format='(f9.6)', dumr
	strctnm.dla_zabs = dumr
	readf, 1,  format='(f6.3)', dumr
	strctnm.dla_NHI = dumr
	readf, 1,  format='(f6.3)', dumr
	strctnm.dla_sigNHI = dumr
	readf, 1,  format='(a60)', dumc
	strctnm.dla_abndfil = dumc
	readf, 1,  format='(i2)', dumi
	strctnm.dla_flgFe = dumi
	readf, 1,  format='(f7.3)', dumr
	strctnm.dla_FeH = dumr
	readf, 1,  format='(f7.3)', dumr
	strctnm.dla_sigFeH = dumr
	readf, 1,  format='(i2)', dumi
	strctnm.dla_flgZn = dumi
	readf, 1,  format='(f7.3)', dumr
	strctnm.dla_ZnH = dumr
	readf, 1,  format='(f7.3)', dumr
	strctnm.dla_sigZnH = dumr
  close, 1
return
end
