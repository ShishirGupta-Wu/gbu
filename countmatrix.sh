#!/bin/bash


if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` [Please : run as countmatrix.sh real_orthogroups]"
  exit 0
fi

if [ "$1" == "" ]; then
  echo "Usage: `basename $0` [Please see Help with script.sh -h]"
  echo "Additional Files: `basename $0` [Script do not need any additional file. Please use -h option to see usase]"
  exit 0
fi

# passing variables
#var1="$1"         #This will be the required orthogroup.txt file 

	printf "\n"
	echo Creating countmatrix for text file "$1"
	
	mapfile -t speciesarray < species.txt
	mapfile -t array < "$1".txt
	arraylength=${#array[@]}
	for (( i=1; i<${arraylength}+1; i++ ));
		do
		cluster=$(echo ${array[$i-1]} | cut -d : -f 1 | xargs basename)
		for j in "${speciesarray[@]}"
			do
			count=$(grep -o $j <<< ${array[$i-1]} | wc -l)
			echo -e "$cluster\t$j=$count"
			done
		done > parsed_"$1"_cluster_counts_raw.txt #will be in source folder.
	awk '$1==last {printf ",%s",$2; next} NR>1 {print "";} {last=$1; printf "%s",$0;} END{print "";}' parsed_"$1"_cluster_counts_raw.txt > parsed_"$1"_cluster_counts_modified.txt

	sed -i 's/,/;/g' parsed_"$1"_cluster_counts_modified.txt #will be in source folder.
	sed -i 's/=//g' parsed_"$1"_cluster_counts_modified.txt #will be in source folder.
	while read -r line; do echo $line; sed -i "s/$line//g" parsed_"$1"_cluster_counts_modified.txt; done < species.txt 
	a=$(while read -r line; do printf ";$line"; done < species.txt)
	b=$(printf "Orthogroup$a")
	(echo $b && cat parsed_"$1"_cluster_counts_modified.txt | sed -e "s/\t/;/g") > "$1"_count_matrix.csv
	
	rm -rf parsed_"$1"_cluster*
	
	echo Analysis done!
	printf "\n"

	printf "\nHere is how the count_matrix looks\n"
	head "$1"_count_matrix.csv
	
	perl -p -e 's/;/\t/g' "$1"_count_matrix.csv > "$1"_count_matrix.tsv

	echo result saved as "$1"_count_matrix.csv and "$1"_count_matrix.tsv
	printf "\n"
	



