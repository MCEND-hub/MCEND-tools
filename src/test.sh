#!/bin/bash
BINDIR=./bin
pycmd=python3

echo "LiH static test no frozen core real time propagation to check energy conservation"
START_TIME=$SECONDS
${BINDIR}/mcend ./calc_testing/LiH.static.inp-R > lih-nofc-R.txt
echo "closed-shell nofc time" $(($SECONDS - $START_TIME))"s"
${pycmd} check_energy_conserved.py expec.t


echo "LiH static test no frozen core"
START_TIME=$SECONDS
./bin/mcend ./calc_testing/LiH.static.inp > lih-nofc.txt
echo "closed-shell nofc time" $(($SECONDS - $START_TIME))"s"
mv expec.t expec.t-lih-cs-nofc
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec.t-lih-cs-nofc ./calc_testing/expec.t-lih-cs-ref


#mv finalpsi startpsi
START_TIME=$SECONDS
echo "LiH acf test no frozen core"
./bin/mcend ./calc_testing/LiH.acf.inp > lih-nofc.txt
echo "closed-shell nofc time" $(($SECONDS - $START_TIME))"s"
mv expec.t expec.t-lih-cs-acf-nofc
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec.t-lih-cs-acf-nofc ./calc_testing/expec.t-lih-cs-acf-ref


echo "LiH static test"
START_TIME=$SECONDS
./bin/mcend ./calc_testing/LiH.static.inp-cs > lih-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"
mv expec.t expec.t-lih-cs
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec.t-lih-cs ./calc_testing/expec.t-lih-cs-fc-ref


START_TIME=$SECONDS
./bin/mcend ./calc_testing/LiH.static.inp-os > lih-os.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
mv expec_spinorbital.t expec_spinorbital.t-lih-os
echo "Comparing with the closed-shell values"
${pycmd} compare_expec_v2.py expec.t-lih-cs expec_spinorbital.t-lih-os


#echo ELAPSED_TIME=$(($SECONDS - $START_TIME))"s"
START_TIME=$SECONDS
echo "LiH acf test"
./bin/mcend ./calc_testing/LiH.static.inp-cs > lih-cs.txt
#mv finalpsi startpsi
./bin/mcend ./calc_testing/LiH.acf.inp-cs > lih-acf-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"
mv expec.t expec.t-lih-acf-cs
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec.t-lih-acf-cs ./calc_testing/expec.t-lih-acf-cs-fc-ref


START_TIME=$SECONDS
./bin/mcend ./calc_testing/LiH.static.inp-os > lih-os.txt
#mv finalpsi_spinorbital startpsi_spinorbital
./bin/mcend ./calc_testing/LiH.acf.inp-os > lih-acf-os.txt
mv expec_spinorbital.t expec_spinorbital.t-lih-acf-os
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
${pycmd} compare_expec_v2.py expec.t-lih-acf-cs expec_spinorbital.t-lih-acf-os


START_TIME=$SECONDS
echo "LiH+ static test"
./bin/mcend ./calc_testing/LiH.static.inp-p > lih-p.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"

mv expec_spinorbital.t expec_spinorbital.t-lih-p
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec_spinorbital.t-lih-p ./calc_testing/expec_spinorbital.t-lih-p-ref


echo "H2 static test"
START_TIME=$SECONDS
./bin/mcend ./calc_testing/H2.static.inp-cs > h2-cs.txt
echo "closed-shell time" $(($SECONDS - $START_TIME))"s"
mv expec.t expec.t-h2-cs
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec.t-h2-cs  ./calc_testing/expec.t-h2-cs-ref


START_TIME=$SECONDS
./bin/mcend ./calc_testing/H2.static.inp-os > h2-os.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"
mv expec_spinorbital.t expec_spinorbital.t-h2-os
${pycmd} compare_expec_v2.py expec.t-h2-cs expec_spinorbital.t-h2-os



echo "OH static test"
START_TIME=$SECONDS
./bin/mcend ./calc_testing/OH.static.inp > oh.txt
echo "open-shell time" $(($SECONDS - $START_TIME))"s"

mv expec_spinorbital.t expec_spinorbital.t-oh
echo "Comparing with the reference values"
${pycmd} compare_expec_v2.py expec_spinorbital.t-oh ./calc_testing/expec_spinorbital.t-oh-ref



