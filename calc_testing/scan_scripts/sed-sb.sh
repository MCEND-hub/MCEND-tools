for i in {10..18..2}
do
 # sed -i 's/numthreads = 1/numthreads = '${i}'/' ./freq${i}/OH.acf.inp
 # sed -i 's/0.22/0.'$i'/' ./freq${i}/OH.acf.inp
  sed -i 's/OHvdz42/OHvdz'$i'/' ./freq${i}/oh-acf.sb
done    

