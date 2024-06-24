for i in {10..80..2}
do
  sed -n -e '1,2000p' ./freq${i}/mcend.out  > ./freq${i}/mcend.out-2
done    

