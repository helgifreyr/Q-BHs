# for dir in m=1.0/c2="$1"/w="$2"/"$3"/*; do
for dir in "$@"; do
  echo $dir
  cp extract-data-Q-extremal.m $dir/
  cd $dir
  rm tmp.txt
  math -script extract-data-Q-extremal.m > math-output.txt
  rm extract-data-Q-extremal.m
  cd ../../../../
done
