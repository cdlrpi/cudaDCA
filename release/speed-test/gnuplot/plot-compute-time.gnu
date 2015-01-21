# NOTE: thes ize is given in inches. It can also be given in cm.
set terminal cairolatex pdf size 4,2
set output "~/Dropbox/Research/Publications/conferences/idetc-2015/cudaDCA/figures/speed.tex"

# Line styles
set border linewidth 1
set style line 1 lt rgb "red"   lw 3
set style line 2 lt rgb "green" lw 3
set style line 3 lt rgb "blue"  lw 3
set style line 4 lt rgb "black" lw 3
set style line 5 lt rgb "violet" lw 3

# Set new scaling 
set   autoscale                        # scale axes automatically

# Set key and title
set notitle
set key vert
# set key outside
set key top left
# set key at graph 1.005, graph 1.175
# set key samplen .5
# set key width -1

# Set axis scaling
set autoscale x
set autoscale y
# set format x "$10^{%L}$"
# unset xtic
# set xtic 2000                          # set xtics automatically
# set ytic 0.1                          # set xtics automatically
set xtic auto                          # set xtics automatically
set ytic auto                          # set ytics automatically


set xlabel "Number of Bodies ($N/10^{3}$)"
set ylabel "Compute Time ($ms$)"

plot    "../data.mat" using ($1/1000):($2) title "All OpenMP" with lines ls 1, \
 	"../data.mat" using ($1/1000):($3) title "$1$ Level GPU" with lines ls 2, \
 	"../data.mat" using ($1/1000):($4) title "$3$ Level GPU" with lines ls 3, \
 	"../data.mat" using ($1/1000):($5) title "$6$ Level GPU" with lines ls 4, \
 	"../data.mat" using ($1/1000):($6) title "All GPU" with lines ls 5 
