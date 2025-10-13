set title "MOS CV Results"
set xlabel "Column 1"
set ylabel "Column 6"
set grid
set key outside

plot for [i=1:10] sprintf("output%d.txt", i) using 1:6 with lines title sprintf("Run %d", i)
