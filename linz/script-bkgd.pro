;
; "Background" estimation
;

;
; look at spectra on and off the disk for the second target.
;

f1map094=map_rflag_cube(dmap094,r=1000,/out)

@foxsi-setup-script

get_target_data, 2, d0,d1,d2,d3,d4,d5,d6, /good

t0 = '2-Nov-2012 17:55:00.000'

i0_disk = where( sqrt(d0.hit_xy_solar[0]^2+d0.hit_xy_solar[1]^2) lt 917  )
i0_off  = where( sqrt(d0.hit_xy_solar[0]^2+d0.hit_xy_solar[1]^2) gt 1017 )
i1_disk = where( sqrt(d1.hit_xy_solar[0]^2+d1.hit_xy_solar[1]^2) lt 917  )
i1_off  = where( sqrt(d1.hit_xy_solar[0]^2+d1.hit_xy_solar[1]^2) gt 1017 )
i2_disk = where( sqrt(d2.hit_xy_solar[0]^2+d2.hit_xy_solar[1]^2) lt 917  )
i2_off  = where( sqrt(d2.hit_xy_solar[0]^2+d2.hit_xy_solar[1]^2) gt 1017 )
i3_disk = where( sqrt(d3.hit_xy_solar[0]^2+d3.hit_xy_solar[1]^2) lt 917  )
i3_off  = where( sqrt(d3.hit_xy_solar[0]^2+d3.hit_xy_solar[1]^2) gt 1017 )
i4_disk = where( sqrt(d4.hit_xy_solar[0]^2+d4.hit_xy_solar[1]^2) lt 917  )
i4_off  = where( sqrt(d4.hit_xy_solar[0]^2+d4.hit_xy_solar[1]^2) gt 1017 )
i5_disk = where( sqrt(d5.hit_xy_solar[0]^2+d5.hit_xy_solar[1]^2) lt 917  )
i5_off  = where( sqrt(d5.hit_xy_solar[0]^2+d5.hit_xy_solar[1]^2) gt 1017 )
i6_disk = where( sqrt(d6.hit_xy_solar[0]^2+d6.hit_xy_solar[1]^2) lt 917  )
i6_off  = where( sqrt(d6.hit_xy_solar[0]^2+d6.hit_xy_solar[1]^2) gt 1017 )

bin=1.0
spex0_disk = make_spectrum( d0[i0_disk], bin=bin)
spex0_off  = make_spectrum( d0[i0_off ], bin=bin)
spex1_disk = make_spectrum( d1[i1_disk], bin=bin)
spex1_off  = make_spectrum( d1[i1_off ], bin=bin)
spex2_disk = make_spectrum( d2[i2_disk], bin=bin)
spex2_off  = make_spectrum( d2[i2_off ], bin=bin)
spex3_disk = make_spectrum( d3[i3_disk], bin=bin)
spex3_off  = make_spectrum( d3[i3_off ], bin=bin)
spex4_disk = make_spectrum( d4[i4_disk], bin=bin)
spex4_off  = make_spectrum( d4[i4_off ], bin=bin)
spex5_disk = make_spectrum( d5[i5_disk], bin=bin)
spex5_off  = make_spectrum( d5[i5_off ], bin=bin)
spex6_disk = make_spectrum( d6[i6_disk], bin=bin)
spex6_off  = make_spectrum( d6[i6_off ], bin=bin)

spex_disk = spex6_disk
spex_off  = spex6_off 
spex_disk.spec_n = spex0_disk.spec_n + spex2_disk.spec_n + spex4_disk.spec_n + $
	spex5_disk.spec_n + spex6_disk.spec_n
spex_disk.spec_p = spex0_disk.spec_p + spex2_disk.spec_p + spex4_disk.spec_p + $
	spex5_disk.spec_p + spex6_disk.spec_p
spex_off.spec_n = spex0_off.spec_n + spex2_off.spec_n + spex4_off.spec_n + $
	spex5_off.spec_n + spex6_off.spec_n
spex_off.spec_p = spex0_off.spec_p + spex2_off.spec_p + spex4_off.spec_p + $
	spex5_off.spec_p + spex6_off.spec_p


plot, spex_disk.energy_kev, spex_disk.spec_p, /xlo, /ylo, psym=10, xr=[3.,13.], /xsty, $
	yr=[1.e-1,200.], /ysty, thick=4, charsize=1.2, $
	xtit='Energy [keV]', ytit='Counts/keV', $
	title='On- and off-disk count spectra for second target, D0,2,4,5,6'
oplot, spex_disk.energy_kev, spex_disk.spec_p, psym=10, color=7, thick=4
oplot_err, spex_disk.energy_kev, spex_disk.spec_p, yerr=sqrt(spex_disk.spec_p*bin)/bin, $
	color=7, thick=4
oplot, spex_off.energy_kev, spex_off.spec_p, psym=10, color=6, thick=4
oplot_err, spex_off.energy_kev, spex_off.spec_p, yerr=sqrt(spex_off.spec_p*bin)/bin, $
	color=6, thick=4
legend,['on-disk','off-disk'],color=[7,6],line=0, /right, box=0
legend, ['Area within 50 arcsec of limb is excluded'], box=0, /left

err = spex_disk.spec_p/spex_off.spec_p * sqrt( 1./spex_disk.spec_p + 1./spex_off.spec_p )

plot_err, spex_disk.energy_kev, spex_disk.spec_p/spex_off.spec_p, yerr = err, $
	/xlo, psym=10, yr=[0.,15.], xr=[3.,13.], /xsty, thick=4, charsize=1.2, $
	xtit='Energy [keV]', ytit='Counts/keV', $
	title='Ratio of on- to off-disk count spectra for second target, D0,2,4,5,6'
legend, ['Area within 50 arcsec of limb is excluded'], box=0, /left


;
; Calculate how much area is on/off disk in 2nd target.
; As usual, exclude region 50" above/below the limb.
;

; set up array of pixels (every detector pixel)
x = reform( rebin( findgen(128), 128, 128 ), 128^2 )
y = reform( transpose( rebin( findgen(128), 128, 128 ) ), 128^2 )

xerr = [  55.4700,  81.490,   96.360,  87.8900,  48.2700,   49.550,   63.450]
yerr = [-135.977, -131.124, -130.241, -92.7310, -95.3080, -120.276, -106.360]

target_center = [-850,150]
coords = get_payload_coords( transpose([[x],[y]]), 6 )
coords[0,*] += target_center[0]
coords[1,*] += target_center[1]
coords[0,*] -= xerr[1]
coords[1,*] -= yerr[1]

print, n_elements( where( sqrt(coords[0,*]^2 + coords[1,*]^2) lt 917 ) ) / 128.^2
print, n_elements( where( sqrt(coords[0,*]^2 + coords[1,*]^2) gt 1017) ) / 128.^2

;
; Results:
; (see lab book for detail)
;
; Units are counts / sec / arcmin^2
; Energy bins are 3-13 (10 bins total)
;

bkgd_disk = [ 0.000106, 0.00201, 0.00529, 0.00624, 0.00254, 0.000952, 0.000423, 0.000423, 0.000106, 0.0]
bkgd_off  = [ 0.0, 0.000329, 0.000494, 0.000576, 0.000247, 0.000165, 0.0000823, 0.0000823, 0.000165, 0.0]

;
; Background estimates using first target.
;

spec = get_target_spectra( 1, bin=1., /good )

;plot, spec[0].energy_kev, spec[0].spec_p, psym=10, /xlo, /ylog, yr=[1.e-2,1.e1], /ysty, $
;	xr=[3.,13.], /xsty
	
i=where(spec[0].energy_kev ge 4. and spec[0].energy_kev le 15.)
print, total( spec[0].spec_p[i] ) / (15.-4.)
print, total( spec[1].spec_p[i] ) / (15.-4.)
print, total( spec[2].spec_p[i] ) / (15.-4.)
print, total( spec[3].spec_p[i] ) / (15.-4.)
print, total( spec[4].spec_p[i] ) / (15.-4.)
print, total( spec[5].spec_p[i] ) / (15.-4.)
print, total( spec[6].spec_p[i] ) / (15.-4.)

; Results:
;
; D0: 0.0475
; D1: 0.0610
; D2: 0.0413
; D3: 0.0370
; D4: 0.150
; D5: 0.0196
; D6: 0.0270
;
; Average = 0.0548 cts / sec / keV
;
;
; Scale to HPD area of ~30"
;
; Detector active area is (128x7.75")^2 = (992")^2 = 991697 sq arcsec
; HPD area is !pi*(30")^2 = 2827.43
; Ratio is 2827.43 / 984064 = 0.002851
;
; Average cts/sec/keV scaled = 0.000156
;
; Other notes:  HPD of 27" corresponds to 261.8 um on the detector, or 3.49 strips.
; Half of this (i.e. the radius) is 130.9 um.
; encircled area is 5.38e-4 cm2.  x2 for both sides --> 1.08e-3 cm2.
; Vol. is 2.69e-5 cm3.

; RHESSI bkgd counts

restore,'rhessi_spectral_fit_foxsi_flare_September2014.sav',/verbose
en1 = get_edges(ebins, /edges_1)
rhessi = rhessi_eff_area(en1, 0.25, 0) / 9.
hsi_bkgd = average(bkg_all, 2)
hsi_bkgd_cts = hsi_bkgd*rhessi
; by now, units are probably cts/kev/sec.  Just average over the energy range.
;
; Result: average is 0.64 cts / sec / kev.
;
