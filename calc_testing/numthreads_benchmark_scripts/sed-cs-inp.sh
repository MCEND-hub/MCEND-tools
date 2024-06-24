for i in 1 2 4 6 8 12 16
do
 # sed -i 's/numthreads = 1/numthreads = '${i}'/' ./core${i}/LiH.static.inp-cs
  sed -i 's/nrorb = 2/nrorb = '$i'/' ./4${i}1/LiH.static.inp-os
done    

