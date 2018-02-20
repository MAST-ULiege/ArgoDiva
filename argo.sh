#!/bin/bash

#Get a file of what is available in merge profiles
#bash data.sh

HOST=vftp.ifremer.fr
USER='anonymous'
PASSWD=''
ftp -i -n $HOST << EOF
quote USER $USER
quote PASS $PASSWD
bin
cd ifremer/argo
get argo_merge-profile_index.txt
bye
EOF

#index path
index=./argo_merge-profile_index.txt

#tmp path
tmp=./tmp.txt

#remove if present
rm ArgoIdParameters.txt
rm ArgoIdParameters.csv

while IFS=, read file date latitude longitude ocean profiler_type institution parameters date_update
	do
	group=`echo $file | awk -F/ '{print $1}'`
	id=`echo $file | awk -F/ '{print $2}'`
	param=`echo $parameters`
	echo "$group $id $param" >> $tmp
done < $index

#remove unecessary lines
sed -i '/#/d' $tmp
sed -i '/file/d' $tmp


#remove duplicates (note : We work we merged profiles so it could be a could thing to only removes id duplicates but we would then need to retain the most recent form of the id.. (add date_update?), note que ce n'est pas forcement une augmentaiton avec le temps -> ex : 6901865
sort tmp.txt | uniq >> ArgoIdParameters.txt
#rm tmp.txt

#create csv file for R
sed 's/ \+/,/g' ArgoIdParameters.txt > ArgoIdParameters.csv
rm argo_merge-profile_index.txt
rm tmp.txt
