# bedtools getfasta -fi ~/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa -fo test.fa -bed test.bed -nameOnly -s
# bowtie ~/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa.bowtie1_index test.fa -f > test.align.tsv
#python get_detect-seq_info.py \
#    --bowtie_table test.align.tsv \
#    --reference /Users/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa \
#    --bam /Users/zhaohuanan/zhaohn_HD/3.project/2021_DdCBE_topic/20210224_DetectSeq_all_bams/bam/293T-DdCBE-ND6-All-PD_rep2_hg38.MAPQ20.bam \
#    --out test.out.csv

python get_target-seq_info.py \
    --bmat_folder bmat \
    --out test.bmatout.csv
