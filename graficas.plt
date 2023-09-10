set datafile separator ' '
set title "Ecuación del Calor en 2D"
set xlabel 'x'
set ylabel 'y'
set zlabel 'U(x,y)'
splot 'Solution.csv' title ''
set terminal png
set output 'sol.png'
replot

