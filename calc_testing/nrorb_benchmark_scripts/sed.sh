for i in 2 3 4 5 6 7 8 9
do
  sed -i 's/nrorb = 2/nrorb = '$i'/' ./4${i}1/LiH.static.inp-cs
 # sed -i 's/nrorb = 2/nrorb = '$i'/' ./4${i}1/LiH.static.inp-os
done    

