#!/bin/sh

for s in 1 2 3 4;
do
python Movie_scatter_color.py -s $s -r $s -o ./videoscatter$s 
done
