function x_constants,CGS=cgs,KSG=ksg
c= {constants          $
    ,c:2.99792458D10   $    ; cm/s
    ,G:6.6726D-8       $    ; cm^3 g^-1 s^-2
    ,h:6.62618D-27     $    ; erg s
    ,k:1.38066D-16     $    ; erg K^-1
    ,mp:1.672623D-24   $    ; g
    ,me:9.10953D-28    $    ; g
    ,eV:1.602D-12      $    ; ergs
    ,sigma:5.6703D-5   $    ; erg s^-1 cm^-2 K^-4
    ,sigmat:6.6525D-25 $    ; cm^2
    ,Msun:1.989D33     $    ; g
    ,Lsun:3.90D33      $    ; erg s^-1
    ,pc:3.08567802D18  $    ; cm
    ,kpc:3.08567802D21 $    ; cm
    ,Mpc:3.08567802D24 $    ; cm
    ,yr:3.155815D7     $    ; s
    ,Gyr:3.155815D16   $    ; s
    ,mu:0.62           $    ; mean moleculare wieght of astrophysical gas
    ,kms2kpcGyr:0.D    $
    ,msunkpc2cm3:0.D   $
    ,erg2oort:0.D      $
    ,Gmix:0.D          $
}
c.Gmix=c.G*(c.Msun/c.kpc)*1.E-5^2  ; cm to km
c.kms2kpcGyr=(1.E-5*c.kpc)/c.Gyr
c.msunkpc2cm3=(c.kpc/c.Msun)*(c.kpc*c.mp)*c.kpc
c.erg2oort=c.Gyr/c.Msun/c.kpc*(c.Gyr/c.kpc)
if KEYWORD_SET(KSG) then begin 
    c.c=c.c*c.kpc/c.Gyr
    c.G=c.G*(c.Msun/c.kpc)*(c.Gyr/c.kpc)^2
    c.pc=1.E-3
    c.Msun=1.
    c.Lsun=c.Lsun*c.erg2oort/c.Gyr
    c.kpc=1.
    c.Mpc=1.E3
    c.yr=1.E-9
    c.Gyr=1.
endif
return,c
end
