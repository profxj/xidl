PRO GNIRS_TELLURIC_STD, type, loglam = loglam, flux = flux, V = V

IF NOT KEYWORD_SET(V) THEN V = 6.0D
gnirs_telluric_params, type, logR = logR, T = T, logg = logg, M_V = M_V

parsec = 3.086d18
R_sol = 6.96d10
; distance modulus
logd = 0.2*(V-M_V) + 1.0D
D = parsec*10.d^logd
R = R_sol*10.d^logR
; factor converts the kurucz surface flux densities to flux observed on Earth
flux_factor = (R/D)^2

T1 = 3000  + lindgen(28)*250
T2 = 10000 + lindgen(6)*500
T3 = 13000 + lindgen(22)*1000
T4 = 35000 + lindgen(7)*2500

Tk = [T1, T2, T3, T4]
loggk = dindgen(11)*0.5

minT = min(abs(Tk-T[0]), indT)
ming = min(abs(loggk-logg[0]), indg)

longslit_dir = getenv('XIDL_DIR') + '/Spec/Longslit/'
std_file = longslit_dir+ '/calib/standards/kurucz93/kp00/kp00_' + $
  strcompress(string(Tk[indT]), /rem) + '.fits.gz'
; read in standar star spectrum
std = mrdfits(std_file, 1)

loglam = alog10(std.wavelength)

CASE indg OF 
    0:  flux = std.G00
    1:  flux = std.G05
    2:  flux = std.G10
    3:  flux = std.G15
    4:  flux = std.G20
    5:  flux = std.G25
    6:  flux = std.G30
    7:  flux = std.G35
    8:  flux = std.G40
    9:  flux = std.G45
    10: flux = std.G50
    ELSE: message, 'Cant find surfact gravity value in table'
ENDCASE

flux = flux_factor[0]*flux

RETURN
END

