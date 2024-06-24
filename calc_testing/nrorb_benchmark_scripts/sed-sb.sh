for i in 2 3 4 5 6 7 8 9
do
  sed -i 's/-test/-4'$i'1/' ./4${i}1/lih-static.sb
 # sed -i 's/nrorb = 2/nrorb = '$i'/' ./4${i}1/LiH.static.inp-os
done    

