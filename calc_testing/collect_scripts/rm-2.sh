for i in {10..80..2}
do
  cd freq${i} 
  rm mcend.out
  mv mcend.out-2 mcend.out
  cd ..
done


