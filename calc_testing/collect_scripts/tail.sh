for i in {10..80..2}
do
  cd freq${i}
  tail -n +2 expec_spinorbital.t >  expec.t
  cd ..
done


