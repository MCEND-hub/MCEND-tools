for i in 2 3 4 5 6 7 8 9
do
  cd 4${i}1
  sbatch lih-static.sb
  cd ..
done


