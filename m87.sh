ctioga2 --name m87 --line-style no \
	--legend-inside lt:0.05,0.35 \
	/scale 1.2 \
	--xlog \
	--xrange 0:2.5 --yrange 27:19 \
	-l 'Global Hernquist Model' \
	--yoffset 7 \
	--line-style auto \
	--marker none \
	hern87.dat \
	-l 'Central Hernquist Model' \
	--yoffset 6.5 \
	--xscale 0.4 \
	hern87.dat \
	--line-style no \
	--marker auto \
	-l 'Observed WFC3' \
	--yoffset 0 \
	--xscale 1 \
	photo.dat@'$2+0.5:-2.5*log10(($8+0.0011)/192.6)+10:yerr=1/sqrt(($8+0.0011)/1.5*5599*$11)' \
	-x 'Radius [Arcseconds]' \
	-y '$\mu_\mathrm{F275W}$'
