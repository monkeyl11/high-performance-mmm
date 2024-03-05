file='./mygemm.h'

for i in $(seq 48 48 960);
do
    for j in $(seq 48 48 960);
    do
        for k in $(seq 48 48 960);
        do
        printf "KC: $i, MC: $j, NC: $k\n"
        echo "\nKC: $i, MC: $j, NC: $k\n" >> runs.txt

        sed -i 's/^\(#define KC \).*/\1'"$i"'/' $file
        sed -i 's/^\(#define MC \).*/\1'"$j"'/' $file
        sed -i 's/^\(#define NC \).*/\1'"$k"'/' $file
        make run >> runs.txt
        done
    done
done