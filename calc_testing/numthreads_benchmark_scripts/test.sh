#!/bin/bash
echo "LiH static test"
START_TIME=$SECONDS
./mcend ./calc_testing/LiH.static.inp-cs > lih-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"

mv expec.t expec.t-lih-cs

START_TIME=$SECONDS
./mcend ./calc_testing/LiH.static.inp-os > lih-os.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"

mv expec_spinorbital.t expec_spinorbital.t-lih-os
python compare_expec_v2.py expec.t-lih-cs expec_spinorbital.t-lih-os

#echo ELAPSED_TIME=$(($SECONDS - $START_TIME))"s"
START_TIME=$SECONDS
echo "LiH acf test"
./mcend ./calc_testing/LiH.static.inp-cs > lih-cs.txt
mv finalpsi startpsi
./mcend ./calc_testing/LiH.acf.inp-cs > lih-acf-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"


mv expec.t expec.t-lih-acf-cs
START_TIME=$SECONDS
./mcend ./calc_testing/LiH.static.inp-os > lih-os.txt
mv finalpsi_spinorbital startpsi_spinorbital
./mcend ./calc_testing/LiH.acf.inp-os > lih-acf-os.txt
mv expec_spinorbital.t expec_spinorbital.t-lih-acf-os
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
python compare_expec_v2.py expec.t-lih-acf-cs expec_spinorbital.t-lih-acf-os

START_TIME=$SECONDS
echo "LiH+ static test"
./mcend ./calc_testing/LiH.static.inp-p > lih-p.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"

mv expec_spinorbital.t expec_spinorbital.t-lih-p-test
python compare_expec_v2.py expec_spinorbital.t-lih-p-test expec_spinorbital.t-lih-p


#mv finalpsi startpsi
#./mcend LiH.acf.inp-cs    > lih-cs-acf.txt
echo "H2 static test"
START_TIME=$SECONDS
./mcend ./calc_testing/H2.static.inp-cs > h2-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"

mv expec.t expec.t-h2-cs

START_TIME=$SECONDS
./mcend ./calc_testing/H2.static.inp-os > h2-os.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"

mv expec_spinorbital.t expec_spinorbital.t-h2-os
python compare_expec_v2.py expec.t-h2-cs expec_spinorbital.t-h2-os



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


