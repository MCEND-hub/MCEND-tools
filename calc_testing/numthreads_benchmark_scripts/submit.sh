for i in 1 2 4 6 8 12 16
do
  cd core${i}
  sbatch lih-static.sb
  cd ..
done


