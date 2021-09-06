bedtools getfasta -fi ~/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa -fo test.fa -bed test.bed -nameOnly -s
bowtie ~/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.bowtie1_index test.fa -f > test.align.tsv
