for i in 1 2 4 6 8 12 16
do
  sed -i 's/cpus-per-task=1/cpus-per-task='${i}'/' ./core${i}/lih-static.sb
 # sed -i 's/nrorb = 2/nrorb = '$i'/' ./4${i}1/LiH.static.inp-os
done    

