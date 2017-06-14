set terminal png
set output 'all.png'
set title 'All Postprocessors'
set xlabel 'time'
set ylabel 'values'
plot 'ShearModeIIDisp.dat' using 1:2 title 'resid_x' with linespoints, \
 'ShearModeIIDisp.dat' using 1:3 title 'resid_y' with linespoints

set output 'resid_x.png'
set ylabel 'resid_x'
plot 'ShearModeIIDisp.dat' using 1:2 title 'resid_x' with linespoints

set output 'resid_y.png'
set ylabel 'resid_y'
plot 'ShearModeIIDisp.dat' using 1:3 title 'resid_y' with linespoints

