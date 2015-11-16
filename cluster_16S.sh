!/usr/bin/bash
## Script to cluster 16S rRNA gene copies


for i in 1.00 0.99 0.97 0.94 0.91 0.88 0.85 0.82 0.79 0.76 0.73 0.70 0.67 0.64 0.61;
do 
    ../dnaclust_repo_release3/dnaclust -d -l -t 8 -e 999 -i all_refseq_16S.fasta -s $i > all_refseq_16S_$i.cluster
done