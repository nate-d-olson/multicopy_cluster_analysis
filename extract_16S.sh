#!/usr/bin/bash
## Script to extract 16S gene sequences from refseq genomes

extract_16S(){
    prefix=$(echo $fa | sed 's/.fna//')
    if [ ! -e "$prefix.fasta" ]
    then
        echo $prefix
        rnammer -S bacterial -m ssu -xml $prefix-16S.xml -gff $prefix-16S.gff -f $prefix-16S.fasta < $fa
    fi
}

N=8
(
for fa in */*fna; do 
   ((i=i%N)); ((i++==0)) && wait
   extract_16S "$fa" & 
done
)