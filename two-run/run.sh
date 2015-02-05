# $1 = nr
# $2 = frequency
# $3 = rh
# $4 = c2
# $5 = branch
# $6 = # of runs


echo "$1" > parameters
echo "$2" >> parameters
echo "$3" >> parameters
echo "$4" >> parameters
for i in `seq 1 $6`;
do
  ./BHs-phi4-extremal | tee output.txt
done    
mkdir -p ../m=$1/c2=$4/w=$2/$5/rh=$3/
mv *txt ../m=$1/c2=$4/w=$2/$5/rh=$3/
cp *dat ../m=$1/c2=$4/w=$2/$5/rh=$3/

# useful stuff for later
# folders = (*/)
# last_folder = $folders(*/)[-1]
