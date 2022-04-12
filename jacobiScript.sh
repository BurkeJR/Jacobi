#!/bin/bash
# John Burke ~ Comp 233 A ~ Jacobi With Picture Script
echo -e "John Burke ~ Comp 233 A ~ Jacobi With Picture Script\n"
i=4
./jacobiPictureSerial .01 1000000
while [ $i -le 16 ]
do
    mpirun -np $i -machinefile ~/machines-openmpi ./jacobiPictureMPI .01 1000000
    ((i*=2))
done
echo -e "\nAll Done\n\n"