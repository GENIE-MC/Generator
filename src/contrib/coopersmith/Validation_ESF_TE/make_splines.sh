#!/bin/sh
spline_output_suffix=_QE_splines.xml

if [ "$1" = "EFF" ]; then 
  gmkspl -p 14,-14 -t 1000010010,1000010020,1000020040,1000060120,1000100200,1000130270,1000180400,1000260560,1000822080 \
  -n 100 -e 100 -o ./splines/$1$spline_output_suffix --event-generator-list CCQE
else
  gmkspl -p 14,-14 -t 1000060120 -n 100 -e 100 -o ./splines/$1$spline_output_suffix --event-generator-list CCQE
fi

root_output_suffix=_QE_splines.root
gspl2root -f ./splines/$1$spline_output_suffix -p 14 -t 1000060120 -e 100 -o ./splines/$1$root_output_suffix
gspl2root -f ./splines/$1$spline_output_suffix -p -14 -t 1000060120 -e 100 -o ./splines/$1$root_output_suffix
