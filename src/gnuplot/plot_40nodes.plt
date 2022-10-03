set term png
set output "solution_40nodes_anal.png"
set border 3
set xtics nomirror
set ytics nomirror
plot "../data/analyticSolution_40nodes.txt" with lines
# , "../data/linearSolution_40nodes.txt" with lines, "../data/cubicSolution_40nodes.txt" with lines