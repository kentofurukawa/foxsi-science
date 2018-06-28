PRO flare1_ospex, filename, det=det, factor=factor, bin=bin, erange=erange, $
		trange=trange, left=left,right=right,top=top,bottom=bottom, $
		center=center, radius=radius

COMMON FOXSI_PARAM

; default values
default, det, 6
default, factor, 1.0
default, bin, 1.0
default, erange, [5.,12.]
default, trange, [t1_pos2_start,t1_pos2_end]
default, center, flare1
default, radius, 150.

; select detector for analysis (Si only)
case det of
	0: data = data_lvl2_d0
	1: data = data_lvl2_d1
	4: data = data_lvl2_d4
	5: data = data_lvl2_d5
	6: data = data_lvl2_d6
endcase

; select portion of FOV to consider
if keyword_set(left) then dat = area_cut(data, center, radius, /xy,/left) $
	else dat = area_cut(data, center, radius, /xy)
if keyword_set(right) then dat = area_cut(data, center, radius, /xy) $
        else dat = dat
if keyword_set(top) then dat = area_cut(data, center, radius, /xy) $
        else dat = dat
if keyword_set(bottom) then dat = area_cut(data, center, radius, /xy) $
        else dat = dat

; time cuts
t1 = trange[0] + t_launch
t2 = trange[1] + t_launch
time_int = t2 - t1
ltime = factor*time_int		;livetime 
i2 = where(dat.wsmr_time gt t1 and dat.wsmr_time lt t2)
i3 = where(dat.wsmr_time gt t1 and dat.wsmr_time lt t2 and dat.error_flag eq 0)

;recenter map
map = foxsi_image_map(dat[i3],center)
flare_cen=map_centroid(map,threshold=0.1*max(map.data))

loadct,7
popen,filename+'_map', xsi=5, ysi=5
plot_map,map,center=flare_cen,fov=3.
pclose

; make energy spectrum
spec_all = make_spectrum( dat[i2], bin=bin )
spec = make_spectrum( dat[i3], bin=bin )

; select energy range
i = where( spec.energy_kev gt 3. and spec.energy_kev lt 15. )
enm = spec.energy_kev[i]
en_lo = enm - deriv(enm)*0.5
en_hi = enm + deriv(enm)*0.5
en2 = transpose( [[en_lo],[en_hi]] )

eff = total(spec.spec_p[i])/total(spec_all.spec_p[i])
;eff = spec.spec_p[i]/spec_all.spec_p[i]
;eff[ where( spec.spec_p[i] eq 0.)] = 0.
;eff[ where( spec.spec_p[i] eq 0.)] = min( eff[ where( spec.spec_p[i] ne 0.)])

; create response matrix
resp = get_foxsi_effarea( energy_arr=enm, module=det )
max_EA = max(resp.eff_area_cm2[where(resp.eff_area_cm2 ge 0.)])
diag = resp.eff_area_cm2 * eff / max_EA
ndiag = n_elements( diag )
nondiag = fltarr(ndiag,ndiag)
for j=0, ndiag-1 do nondiag[j,j]=diag[j]

fwhm  = 0.5			; energy resolution
sigma = fwhm / 2.355
toty = total(gaussian(findgen(ndiag),[0.3989*bin/sigma,round((12./bin)-1.),sigma/bin] ))
for j=0, ndiag-1 do begin
	y = diag[j]*gaussian( findgen(ndiag), [0.3989*bin/sigma,j,sigma/bin] )/toty
	nondiag[*,j] = y
endfor

livetime_array = fltarr(ndiag)
livetime_array[*] = ltime  

o = ospex(/no_gui)
o -> set, spex_data_source = 'SPEX_USER_DATA'
o -> set, spectrum = spec.spec_p[i]*bin, spex_ct_edges = en2, errors = spec.spec_p_err[i]*bin,livetime=livetime_array
o -> set, spex_respinfo = nondiag
o -> set, spex_area = max_EA

o-> set, fit_function= 'vth'
o-> set, fit_comp_params= [0.003, 1.0, 1.00000]                                
o-> set, fit_comp_free_mask= [1B, 1B, 0B]
o-> set, fit_comp_spectrum= ['full', '']
o-> set, fit_comp_model= ['chianti', '']
o-> set, spex_erange = erange
o-> set, spex_fit_time_interval = [0.,3600.]

stop

o-> set, spex_fit_manual=0
o-> set, spex_autoplot_enable=1

o-> dofit
o-> savefit, outfile = filename+'.fits'

p0 = spex_read_fit_results(filename+'.fits')
temp = string(p0.spex_summ_params[1]*11.6,format='(f20.1)')
temp = strtrim(temp,2)
em_exp = floor(alog10(p0.spex_summ_params[0]))
em_base = string(p0.spex_summ_params[0]/(10.^em_exp),'(f20.2)')
em_base = strtrim(em_base,2)
em_exp_string = strtrim(49+em_exp,2) 
temp_legend = 'T = '+temp+' MK'
em_legend = 'EM = '+em_base+'x10!U'+em_exp_string+'!N cm!U-3!N'

popen,filename, xsi=7, ysi=4
ch = 1.2
en = get_edges( p0.spex_summ_energy, /mean )
conv = p0.spex_summ_conv*p0.spex_summ_area
counts = p0.spex_summ_ct_rate/bin 
counterr = sqrt(counts/(ltime*bin))
big = where(counterr gt counts)
counterr[big] = counts[big]*0.9
fit = p0.spex_summ_ph_model*conv
linecolors
plot, en, counts, psym=10, /xlo, /ylo, xr=[4.,10.], yra=[1.e-1,1.e1], thick=8, /xsty, xth=2, yth=2, $
	xtit='Energy [keV]', ytit='Counts s!U-1!N keV!U-1!N', tit='FOXSI-2 First Microflare', charsi=ch, $
	charth=2   
oplot_err, en, counts, yerr=counterr, psym=10, thick=8
oplot, en, fit, psym=10, clip=[erange[0],1.e-1,erange[1],1.e1], col=3, thick=8
oplot, [erange[0],erange[0],erange[1],erange[1]], [1.e-1,1.e1,1.e1,1.e-1], line=1, thick=8, col=3
al_legend, [em_legend, temp_legend], /right, box=0, charsi=1.3,background_color='white'
pclose

;stop

;popen,filename+'_eff',xsi=7,ysi=4
;plot, en, eff, xtit='Energy (keV)', ytit='Efficiency',thick=8
;pclose

print, spec.spec_p[i]
print, total(spec.spec_p[i])
print, spec_all.spec_p[i]
print, eff

end
