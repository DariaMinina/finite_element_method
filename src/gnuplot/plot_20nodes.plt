set term png
set output "solution_20nodes_anal.png"
set border 3
set xtics nomirror
set ytics nomirror
plot "../data/analyticSolution_20nodes.txt" with lines
# , "../data/linearSolution_20nodes.txt" with lines, "../data/cubicSolution_20nodes.txt" with lines