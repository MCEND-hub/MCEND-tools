for i in 1 2 4 6 8 12 16
do
  sed -i 's/-test/-'$i'/' ./core${i}/lih-static.sb
done
