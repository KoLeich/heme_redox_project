#!bin/bash


cd ./database/logfilessplit
rm time.txt
for d in * ; do
	echo "$d" >> time.txt
    tail "$d" --lines=4 | head --lines=1 >> time.txt
    echo " " >> time.txt
done

