file=$1

if [ ! -f $file ];
then
  echo "file $file does not exist"
else
  echo "file $file exists"
fi
