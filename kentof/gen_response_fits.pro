PRO gen_response_fits, outfilename=outfilename, energy_bin_arr=energy_bin_arr, bin_width=bin_width


COMMON FOXSI_PARAM

; default values
default, outfilename, 'response.fits'
default, energy_bin_arr, [4.,5.,6.,7.,8.,9.,10.,11.,12.,13.,14.,15.,16.,17.,18.,19.,20.]
default, det, 6
default, bin_width, 1.0


; select detector for analysis (Si only)
case det of
	0: data = data_lvl2_d0
	1: data = data_lvl2_d1
	4: data = data_lvl2_d4
	5: data = data_lvl2_d5
	6: data = data_lvl2_d6
endcase


trange = [t1_pos2_start,t1_pos2_end]
factor = 1.

; time cuts
t1 = trange[0] + t_launch
t2 = trange[1] + t_launch
time_int = t2 - t1
ltime = factor*time_int		;livetime 
i2 = where(data.wsmr_time gt t1 and data.wsmr_time lt t2)
i3 = where(data.wsmr_time gt t1 and data.wsmr_time lt t2 and data.error_flag eq 0)

; make energy spectrum
spec_all = make_spectrum( data[i2], bin=bin_width )
spec = make_spectrum( data[i3], bin=bin_width )

; select energy range
i = where( spec.energy_kev gt 3. and spec.energy_kev lt 15. )

print, spec.spec_p[i]
print, spec_all.spec_p[i]

eff = total(spec.spec_p[i])/total(spec_all.spec_p[i])
;eff = spec.spec_p[i]/spec_all.spec_p[i]
;eff[ where( spec.spec_p[i] eq 0.)] = 0.
;eff[ where( spec.spec_p[i] eq 0.)] = min( eff[ where( spec.spec_p[i] ne 0.)])



bin_edges_arr = get_edges(energy_bin_arr)
print, energy_bin_arr
print, bin_edges_arr.mean

; create response matrix
resp = get_foxsi_effarea( energy_arr=bin_edges_arr.mean, module=det )
print, resp.eff_area_cm2
max_EA = max(resp.eff_area_cm2[where(resp.eff_area_cm2 ge 0.)])
print, max_EA
diag = resp.eff_area_cm2 * eff / max_EA
print, diag
ndiag = n_elements( diag )
nondiag = fltarr(ndiag,ndiag)
for j=0, ndiag-1 do nondiag[j,j]=diag[j]


fwhm  = 0.5; energy resolution
sigma = fwhm / 2.355
toty = total(gaussian(findgen(ndiag),[0.3989*bin_width/sigma,round((12./bin_width)-1.),sigma/bin_width] ))
print, toty
for j=0, ndiag-1 do begin
   y = diag[j]*gaussian( findgen(ndiag), [0.3989*bin_width/sigma,j,sigma/bin_width] )/toty
   nondiag[*,j] = y
endfor

print, nondiag
surface, nondiag
mwrfits, nondiag, outfilename, /create

end
