#!/bin/bash


x=128
y=128
z=128

for (( j=$x; j>=2; ))
do
    for (( k=$y;k>=2; ))
    do
        for (( i=$z;i>=2; ))
        do 
	        echo  "$i $k $j "
            i=`expr $i / 2`
            
            sed -i "s#BLOCK_SIZE_x=.*#BLOCK_SIZE_x=${i}#g" clover.in
            sed -i "s#BLOCK_SIZE_y=.*#BLOCK_SIZE_y=${k}#g" clover.in
            sed -i "s#BLOCK_SIZE_z=.*#BLOCK_SIZE_z=${i}#g" clover.in
            ./clover_leaf 
            cp clover.out clover.out.${i}.${k}.${j}
            
        done
        k=`expr $k / 2`
    done
    j=`expr $j / 2`
done
 
echo "Boom!"
