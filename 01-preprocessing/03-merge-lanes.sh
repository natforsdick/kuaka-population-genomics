#!/bin/bash -e

cd /nesi/nobackup/ga03186/kuaka-pop-gen/output/02-trimmed/02b-trimgalore/

for i in $(find ./ -type f -name "EXT*val_1.fq.gz" | while read F; do basename $F | rev | cut -c 18- | rev; done | sort | uniq)
	do 
	echo "Merging ${i}"
	
	cat "$i"_L00*_val_1.fq.gz > ../../03-merged/"$i"_val_1.fq.gz

	cat "$i"_L00*_val_2.fq.gz > ../../03-merged/"$i"_val_2.fq.gz

done
