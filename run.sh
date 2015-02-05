# $1 = nr
# $2 = frequency
# $3 = rh
# $4 = alpha
# $5 = branch
# $6 = # of runs


echo "$1" > parameters
echo "$2" >> parameters
echo "$3" >> parameters
echo "$4" >> parameters
for i in `seq 1 $6`;
do
  ./Q-BHs | tee output.txt
done    
math -script extract-data-Q-HBHs.m | tee math-output.txt
mkdir -p m=$1/alpha=$4/w=$2/$5/rh=$3/
mv *txt m=$1/alpha=$4/w=$2/$5/rh=$3/
cp *dat m=$1/alpha=$4/w=$2/$5/rh=$3/

# useful stuff for later
# folders = (*/)
# last_folder = $folders(*/)[-1]
# test test test
