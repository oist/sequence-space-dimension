# Set the output to a png file
set terminal png
set xrange [0:1]
set yrange [0:600]
set nokey
# The file we'll write to
set output 'Drawing_curve_b0_simulate_evol/SimuOG8862_drawing.png'
# The graphic title
set title 'Betti 0 for SimuOG8862   n = 295'
#plot the graphic
plot "Result_b0curve_simulate_evol/SimuOG8862.b0curve" with steps lw 2