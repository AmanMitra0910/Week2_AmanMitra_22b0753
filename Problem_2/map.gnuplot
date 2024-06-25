set pm3d map
set palette defined ( 0 "blue", 1 "cyan", 2 "green", 3 "yellow", 4 "red" )
set title "Potential Distribution"
set xlabel "x"
set ylabel "y"
set dgrid3d 50,50
splot 'potential_distribution.txt' matrix with pm3d
