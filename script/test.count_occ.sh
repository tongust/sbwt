if [ $# -ne 3 ]; then
        echo "usage: $0 [fa] [period] [length of seed]"
        exit 1
fi


rm -rf count_occ
rwd=$PWD

make clean

make count_occ

#mv build_index /home/tongust/Projects/test/
#mv count_occ /home/tongust/Projects/test/

make clean

#cd $rwd

./count_occ  $2 $1 $3 > stats.$1.$2.$3.csv
