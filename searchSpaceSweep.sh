#!/bin/bash


x=512
y=512
z=512

l3=20971520


for (( j=$x; j>=2; ))
do
    for (( k=$y;k>=2; ))
    do
        for (( i=$z;i>=2; ))
        do 
	        echo  "$i $k $j "
            i=`expr $i / 2`
        done
        k=`expr $k / 2`
    done
    j=`expr $j / 2`
done
 
echo
echo "Boom!"
