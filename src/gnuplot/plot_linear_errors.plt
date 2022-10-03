set term png
set output "graphics/errors_linear_solution.png"
set border 3
set xtics nomirror
set ytics nomirror
plot "data/errors_linearSolution_Nnodes.txt" with lines