#!/bin/bash
#NAME=/baby/one/more/time
#echo $NAME>name.log
FULLNAME=$(cat name.log)
echo $FULLNAME
POS=$(dirname $FULLNAME)
echo $POS
NAME=$(basename $POS)
echo $NAME
for thrd in 1 3 5 10 20 30
do
echo $thrd
hadd $POS/$NAME.root $POS/{0..39}data_$thrd.root
root -b -q "MultiTiers_TACor_sep.C(\"$POS/$NAME\")"
done
