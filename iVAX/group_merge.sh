#!/bin/bash
cwd=$(pwd)
## first process every 20 html file in that folder
i=0 
for f in * 
do 
	d=dir_$(printf %03d $((i/20+1)))
	mkdir -p $d
	mv "$f" $d
	let i++ 
done

## loop through all the sub dir in that dir to merge the html


for d in *
do
	echo "Find html in $d "
	html=$(find $d -type f -name '*.html')
	##echo "$html"
 	for f in $html
	do
		#echo "I am in $d to merge html."
		tail -n +48 $f | sed '$d' | sed '$d'  >>  $d/MergedFinalData.html
		echo "<BR style="page-break-after: always">"  >>  $d/MergedFinalData.html
	done
	
done

