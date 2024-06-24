for i in {10..18..2}
do
  cd freq${i}
  sbatch oh-acf.sb
  cd ..
done


