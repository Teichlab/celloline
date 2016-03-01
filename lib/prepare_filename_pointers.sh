#!/bin/sh

# this is helper script for making fastq files ready for celloline. It will not fit all situations. It might usefull for inspiration, tough. 
# run the script in the /raw directory, which you place in your input directory (parameter -i), which is given relative to your root (provided in config.) 
# change the function if you wish to decompress or move etc. on the fly as well



shopt -s globstar # this will make it possible to wildcard paths using **, if your files are in several libarites


# this function will take two arguments, a path/filename and a (cell)number to be added
# it assumes that  
link() {
	number=$1
	file=$2
	cellname="singleCell"
	dest=`basename "$file" | sed "s/.*\_\([12]\).fq/$cellname#$number_\1.fq/"`
	echo -n $file
	echo $dest

	ln -s $file $dest
	
}
export -f link


# this adds cell numbers from 1 to Nfiles to files. Actual regex following 'ls' in the two arguments must be the same
# it will print a list of old and new names. Save this list in order to be able to reidentify your cells.
# seq is the bash command that will make the numbers for the files. Do change the path/regex/glob according to where your files are.
# add --dry-run parameter to parallel (before "link") to see the effect - but without actually doing it.
# parallel is available at 
parallel -j16 -k --xapply link ::: $(seq 1 `ls ./**/raw/*.fq | wc -l`) ::: `ls ./**/raw/*.fq`
