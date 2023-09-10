#!/bin/bash


if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` [Please Enter the four or six letter abbrivation of both the species: run as script.sh afum dmel]"
  exit 0
fi

if [ "$1" == "" ]; then
  echo "Usage: `basename $0` [Please see Help with script.sh -h]"
  echo "Additional Files: `basename $0` [Script needs two perl files dodo_result_align.pl and one2many.pl]"
  exit 0
fi

# for awk passing variables
#var1="$1"
#var2="$2"


printf "\n"
echo Parsing_seed for Species "$1" and "$2"


mkdir seed_sql_"$1"_"$2"
cp sqltable."$1".fasta-"$2".fasta seed_sql_"$1"_"$2"/
cd seed_sql_"$1"_"$2"

cat sqltable."$1".fasta-"$2".fasta | awk '{if ($4==1.000) print $0}' > temp1 
awk '{print $1, $5}' temp1  > temp2 

dodo_result_align.pl temp2 > temp3 
cut -d' ' -f2- temp3 > temp4

one2many.pl temp4 "$1" > temp5
grep -v 'Line' temp5 > nrd

printf "\n"
echo Total 1 to 1 ortholog
perl -lne 'END { print $. }' nrd

mv nrd FINAL_seed_"$1"_"$2"_inparanoid_sql.txt

printf "\n"
echo Creating 1 to 1 ortholog bitscore file

awk '{print $1, $2, $5}' temp1 > temp6 
dodo_result_align.pl temp6 > temp7 
cut -d' ' -f2- temp7 > temp8
one2many.pl temp8 "$1" > temp9
grep -v 'Line' temp9 > temp10
sort temp10 | uniq > temp11
sort -n -r -k2 temp11 > temp12
grep -v "$2" temp12 > temp13

#awk /"$var1"/&&/"$var2"/ temp12 > temp14
grep "$1" temp12 | grep "$2" > temp14



awk 'NR==FNR{a[$1]=$2}{print $0, a[$1]}' OFS='\t' temp13 temp14 > temp15

#awk /"$var1"/&&/"$var2"/  temp15 > temp16
grep "$1" temp12 | grep "$2" temp15 > temp16

sort -n -r -k3 temp16 > temp17

printf "\n"
echo Total 1 to 1 ortholog bitscore file
perl -lne 'END { print $. }' temp17

mv temp17 FINAL_seed_"$1"_"$2"_inparanoid_sql_BitScore.txt

printf "\n"
echo Analysis done!
printf "\n"
echo result saved as FINAL_seed_"$1"_"$2"_inparanoid_sql.txt in seed_sql_"$1"_"$2" directory
printf "\n"
echo result saved as FINAL_seed_"$1"_"$2"_inparanoid_sql_BitScore.txt in seed_sql_"$1"_"$2" directory
printf "\n"

mkdir intermediates_files
mv temp* intermediates_files/

cd ..



