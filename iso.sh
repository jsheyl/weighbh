ctioga2 --name iso --xlog  \
	--legend-inside lt:0.05,0.45 \
	/scale 1.2 \
	-l 'Global Hernquist Model' \
	--yoffset 0.9 \
	hern.dat \
	-l 'Central Hernquist Model' \
	--xscale 0.1 \
	--yoffset -2.2 \
	hern.dat \
	--xscale 1 \
	--yoffset 0 \
	-l 'Lauer et al.\ Observed WPFC2' \
	wfpc2o.txt \
	-l 'Lauer et al.\ Deconvolved WPFC2' \
	wfpc2d.txt \
	-l 'Michard \& Nieto $B-$Band' \
	--line-style no --marker auto \
	--yoffset -1 \
        nieto.dat \
	--yoffset 0 \
	-l 'Observed WFC3' \
 'iso.dat@$8*0.04:-2.5*log($6)/log(10)+18.3:yerr=2*1.08/sqrt($6*3.14*(3*($7+$8)-sqrt((3*$7+$8)*($7+3*$8))))' \
	--yrange 18:10.5  --xrange -2:1.2 \
	-x 'Semimajor Axis [arcseconds]' -y '$\mu_\mathrm{F555W}$' 
