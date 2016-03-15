#!/bin/bash


x=96
y=96
z=96

for (( j=$x; j>=2; ))
do
    for (( k=$y;k>=2; ))
    do
        for (( i=$z;i>=2; ))
        do 
	        echo  "$i $k $j "
            
            sed -i "s#block_size_x=.*#block_size_x=${i}#g" clover.in
            sed -i "s#block_size_y=.*#block_size_y=${k}#g" clover.in
            sed -i "s#block_size_z=.*#block_size_z=${j}#g" clover.in
            ./clover_leaf 
            cp clover.out clover.out.${i}.${k}.${j}
            i=`expr $i / 2`
        done
        k=`expr $k / 2`
    done
    j=`expr $j / 2`
done
 
echo "Boom!"
