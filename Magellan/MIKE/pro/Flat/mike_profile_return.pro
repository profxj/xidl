
function mike_profile_return, slit_frac, y, ordr_p
 
     nrow = n_elements(ordr_p.lhedg)
     ynorm      = (2.0* y- nrow)/nrow

     p0 = ordr_p.profile0
     p1 = ordr_p.profile1

     if total(finite(p1) EQ 0) GT 0 then p1[*] = 0. 

     npoints = n_elements(p0)
     xp      = (npoints-1)/2 + 100.0 * slit_frac
     inter0 = interpolate(p0, xp)
     inter1 = interpolate(p1, xp)

     profile = (inter0 + inter1 * ynorm) * (xp GE 0 AND xp LE npoints-1)
     return, profile
end

