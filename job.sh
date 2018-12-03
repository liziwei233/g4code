#!/bin/bash

#cat backup/output.C
source /data2/R710liziwei/.bashrc

for thrd in 1 3 5 10 20 30
do
filename="testthrd$thrd.condo"
echo $filename
cp test.condo $filename
#sed -n "/Arguments/p" test.condo 
sed -i "/Arguments/s/thrd/$thrd/" $filename
cs $filename
#sed -i 's//../g' test.condo > test$.condo
done



