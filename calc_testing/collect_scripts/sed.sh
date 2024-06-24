for i in {10..20..2}
do
 # sed -i 's/numthreads = 1/numthreads = '${i}'/' ./freq${i}/OH.acf.inp
  sed -i 's/0.20/0.'$i'/' ./freq${i}/B2.acf.inp
#  sed -i 's/  /'                 ./freq${i}$/OH.acf.inp
done    

