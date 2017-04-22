#!/bin/bash

for i in `seq 1 9` ; do
gnuplot <<EOF
set hidden3d
set term png
set output "frame$i.png"
splot "out.asc" u 1:2:4 i $i w l
EOF
done
