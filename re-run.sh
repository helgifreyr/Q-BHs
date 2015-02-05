# for dir in m=1.0/c2="$1"/w="$2"/"$3"/*; do
for dir in "$@"; do
  if [ ! -f $dir/tmp.txt ]; # comment out this line to run on already run folders
  then # and this
    echo $dir
    cp extract-data-Q-HBHs.m $dir/
    cd $dir
    rm tmp.txt
    math -script extract-data-Q-HBHs.m > math-output.txt
    rm extract-data-Q-HBHs.m
    cd ../../../../../
  fi # and this
done
