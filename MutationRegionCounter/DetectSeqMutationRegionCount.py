import pysam
import pandas as pd
from DetectSeqLib_V2.CheckAndLoadFiles import load_reference_fasta_as_dict
from DetectSeqLib_V2.RegionStatsTest import get_align_mismatch_pairs
from DetectSeqLib_V2.RegionStatsTest import get_No_MD_align_mismatch_pairs


# read params
PATH_OUT_TABLE = './test.out.csv'
PATH_REF_FASTA = "/Users/mac/Nutstore/Coding/github/snakepipes_bioinformatics_hermanzhaozzzz/genome_fa/genome_ucsc_hg38.fa"
PATH_INPUT_BAM = "./bam/293T-bat_VEGFA-All-PD_rep1_hg38.MAPQ20.bam"
PATH_INPUT_BED = './bed/GSE151265_293T-VEGFA-Detect-seq_pRBS.bed'
BL_HAS_BED_HEADER = True
NUM_EXTEND_LENGTH = 100
# MUT_DIRECTION = "CT, GA"


# load file
ref_dict = load_reference_fasta_as_dict(ref_fasta_path=PATH_REF_FASTA)
bam_file = pysam.AlignmentFile(PATH_INPUT_BAM, "rb")
# ls_mut_direction = [i.strip() for i in MUT_DIRECTION.split(",")]


with open(PATH_INPUT_BED, "r") as f:
    ls_bed_info = [i.strip().split('\t') for i in f.readlines()]

if BL_HAS_BED_HEADER:
    ls_bed_info = ls_bed_info[1:]



ls_table = []
for ls_bed_line in ls_bed_info:
    chr_ = ls_bed_line[0]
    start = int(ls_bed_line[1]) - NUM_EXTEND_LENGTH
    end = int(ls_bed_line[2]) + NUM_EXTEND_LENGTH
    region_index = ls_bed_line[3].strip()

    for index, align in enumerate(bam_file.fetch(contig=chr_, start=start, end=end)):
        MD_tag_state = align.has_tag("MD")

        if MD_tag_state:
            mis_align_pair_list = get_align_mismatch_pairs(align)
        else:
            mis_align_pair_list = get_No_MD_align_mismatch_pairs(align, ref_dict)

        # 如果这个reads至少存在一个mutation信息    
        if mis_align_pair_list:
            # checked
            align_strand = "-" if align.is_reverse else "+"
            aligned_reads_id = index + 1
            extend_length = NUM_EXTEND_LENGTH

            # 如果只存在一个mutation，不做判断，直接存入list
            if len(mis_align_pair_list) == 1:
                ls_line = []
                if mis_align_pair_list[0][2] == 'C' and mis_align_pair_list[0][3] == 'T':
                    ls_line = [region_index, aligned_reads_id, align_strand, extend_length, 1, 0, 0, 0]
                elif mis_align_pair_list[0][2] == 'G' and mis_align_pair_list[0][3] == 'A':
                    ls_line = [region_index, aligned_reads_id, align_strand, extend_length, 0, 0, 1, 0]

            # 如果存在两个或两个以上的mutation        
            else:
                ls_line1 = []
                ls_line2 = []
                # 先把CT的全筛出来
                ls_this_reads_info = [site for site in mis_align_pair_list if site[2] == 'C' and site[3] == 'T']
                # 如果筛不到，没有就没有
                # 如果筛到至少一个
                if ls_this_reads_info:
                    # 如果有一个就加到ls_table中
                    if len(ls_this_reads_info) == 1:
                        ls_line = [region_index, aligned_reads_id, align_strand, extend_length, 1, 0, 0, 0]
                    # 如果有更多再计算count和gap
                    else:  # todo
                        ls_first = ls_this_reads_info[0]
                        ls_last = ls_this_reads_info[-1]
                        #                         print(ls_first[1])# todo!
                        ls_line1 = [region_index, aligned_reads_id, align_strand, extend_length,
                                    len(ls_this_reads_info), ls_last[0] - ls_first[0], 0, 0]
                # 再把GA的全筛出来
                ls_this_reads_info = [site for site in mis_align_pair_list if site[2] == 'G' and site[3] == 'A']
                # 如果筛不到，没有就没有
                # 如果筛到至少一个
                if ls_this_reads_info:
                    # 如果有一个就加到ls_table中
                    if len(ls_this_reads_info) == 1:
                        ls_line = [region_index, aligned_reads_id, align_strand, extend_length, 0, 0, 1, 0]
                    # 如果有更多再计算count和gap
                    else:  # todo
                        ls_first = ls_this_reads_info[0]
                        ls_last = ls_this_reads_info[-1]
                        #                         print(ls_first[1])# todo!
                        ls_line2 = [region_index, aligned_reads_id, align_strand, extend_length, 0, 0,
                                    len(ls_this_reads_info), ls_last[0] - ls_first[0]]
                if ls_line1 and ls_line2:
                    ls_line = [region_index, aligned_reads_id, align_strand, extend_length, ls_line1[4] + ls_line2[4],
                               ls_line1[5] + ls_line2[5], ls_line1[6] + ls_line2[6], ls_line1[7] + ls_line2[7]]
                elif ls_line1:
                    ls_line = ls_line1
                elif ls_line2:
                    ls_line = ls_line2
                else:
                    pass

                if ls_line:
                    ls_table.append(ls_line)

                # 既没有CT也没有GA的直接pass

        if index > 1000:
            break
df = pd.DataFrame(ls_table, columns=["region_index", "reads_id", "align_strand", "extend_length", "CT_count",
                                     "Ct_max_gap", "GA_count", "GA_max_gap"]
                  )
df.to_csv(PATH_OUT_TABLE, index=False, sep='\t')
