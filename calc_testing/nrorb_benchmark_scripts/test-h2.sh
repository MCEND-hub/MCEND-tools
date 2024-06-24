./mcend H2.static.inp-cs > h2-cs.txt
mv expec.t expec.t-h2-cs
#mv expec.t expec.t-lih-cs
#mv finalpsi startpsi
#./mcend LiH.acf.inp-cs    > lih-cs-acf.txt
./mcend H2.static.inp-os > h2-os.txt
mv expec_spinorbital.t > expec_spinorbital.t-h2-os
#mv expec_spinorbital.t expec_spinorbital.t-lih-os
#./mcend H2.static.inp > h2p.txt
#mv finalpsi_spinorbital startpsi_spinorbital 
#./mcend LiH.acf.inp-os    > lih-os-acf.txt


