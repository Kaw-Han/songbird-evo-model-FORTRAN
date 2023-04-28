#!/usr/bin/env gnuplot

set datafile separator ","

set term svg font "Calibri, 22"

set xrange [1:100]
set xlabel 'Generation' font "Calibri, 22"

#---------------------------------------------------------
# COLUMNS in data file:
#
# ["GENERATION                 ",  &  !1
#  "NUM-ALIVE                  ",  &  !2
#  "AVERAGE-MASS-ALL           ",  &  !3
#  "AVERAGE-MASS-ALIVE         ",  &  !4
#  "BEST-MASS                  ",  &  !5
#  "AVERAGE-MASS-DOMINANTS     ",  &  !6
#  "AVERAGE-GENE-ALIVE         ",  &  !7
#  "STD-DEV-MASS-ALIVE         "]     !8
#---------------------------------------------------------

# set yrange [1:100000]

set output "plot_gens-1.svg"
set ylabel 'Number of birds alive' font "Calibri, 22"
#plot 'my_model_output.csv' using (column("NUM-ALIVE")) with lines
plot 'my_model_output.csv' using 2 notitle with lines lw 2 lc "blue"

# set yrange [1:45]

set output "plot_gens-2.svg"
set ylabel 'Average mass, all birds' font "Calibri, 22"
#plot 'my_model_output.csv' using (column("AVERAGE-MASS-ALL")) with lines
plot 'my_model_output.csv' using 3 notitle with lines lw 3 lc "blue"


set output "plot_gens-3.svg"
set ylabel 'Average mass, alive birds'font "Calibri, 22"
#plot 'my_model_output.csv' using (column("AVERAGE-MASS-ALIVE")) with lines
plot 'my_model_output.csv' using 4 notitle with lines lw 3 lc "blue" 

set output "plot_gens-4.svg"
set ylabel 'Average mass, best 25% of birdsbirds' font "Calibri, 22"
#plot 'my_model_output.csv' using (column("BEST-MASS")) with lines
plot 'my_model_output.csv' using 5 notitle with lines lw 3 lc "blue"

set output "plot_gens-5.svg"
set ylabel 'Average gene value, alive birds' font "Calibri, 22"
#plot 'my_model_output.csv' using (column("AVERAGE-GENE-ALIVE")) with lines
plot 'my_model_output.csv' using 7 notitle with lines lw 4 lc "blue"

set output "plot_gens-6.svg"
set ylabel 'Standard deviaion of mass, alive birds' font "Calibri, 22"
#plot 'my_model_output.csv' using (column("STD-DEV-MASS-ALIVE")) with lines
plot 'my_model_output.csv' using 8 notitle with lines lw 3 lc "blue"

