#!/bin/bash


echo "OH static test"
START_TIME=$SECONDS
./mcend ./calc_testing/OH.static.inp > oh.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"

mv expec_spinorbital.t expec_spinorbital.t-oh
python compare_expec_v2.py expec_spinorbital.t-oh expec_spinorbital.t-oh-ref


#echo ELAPSED_TIME=$(($SECONDS - $START_TIME))"s"
#./mcend H2.static.inp > h2p.txt
#mv finalpsi_spinorbital startpsi_spinorbital 
#./mcend LiH.acf.inp-os    > lih-os-acf.txt


