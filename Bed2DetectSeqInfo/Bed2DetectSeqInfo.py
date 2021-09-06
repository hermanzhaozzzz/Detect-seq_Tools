import logging

import pandas as pd
import pysam

from Bed2DetectSeqInfo.loader import (
    load_reference_fasta_as_dict,
    load_region_absolute_position,
)
from Bed2DetectSeqInfo.parse_pysam import (
    get_align_mismatch_pairs,
    get_query_site_mpileup_info,
)

"""initialize"""
pd.set_option("max_colwidth", 300)  # column最大宽度
pd.set_option("display.width", 1000)  # dataframe宽度
pd.set_option("display.max_columns", None)  # column最大显示数
pd.set_option("display.max_rows", 10)  # row最大显示数


pt_aligninfo = "Bed2DetectSeqInfo/test.align.tsv"
# pt_ref = "/Users/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
pt_ref = "/Users/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/chr20.fa"
pt_bam = "/Users/zhaohuanan/zhaohn_HD/3.project/2021_DdCBE_topic/20210224_DetectSeq_all_bams/bam/293T-DdCBE-ND6-All-PD_rep2_hg38.MAPQ20.bam"
bam_file = pysam.AlignmentFile(pt_bam, "rb")


if __name__ == "__main__":
    # load ref
    dt_ref = load_reference_fasta_as_dict(
        ref_fasta_path=pt_ref, ref_name=["chr20"], log_verbose=logging.DEBUG
    )
    # load bam
    bam_file = pysam.AlignmentFile(pt_bam, "rb")
    # load region info
    df_region = load_region_absolute_position(bowtie_table=pt_aligninfo)

    for row in df_region.iterrows():
        # walk through region rows
        # parse region info
        idx_row = row[0]
        info_row = row[1]
        chrom = info_row["chrom"]
        start = info_row["start"]
        stop = info_row["stop"]
        strand = info_row["strand"]
        region_name = info_row["read_id"].split("(")[0]
        seq_ref_strand = info_row["sequence"]
        # for debug, only one region tested
        if chrom != "chr20":
            continue

        bam_file
        for idx_abs_base, pileup_base in enumerate(
            bam_file.pileup(
                contig=chrom,
                start=start,
                stop=stop,
            )
        ):
            # restrict the pileup_base within targeted region
            if pileup_base.pos < start:
                continue
            elif pileup_base.pos > stop:
                break

            # pileup_base.pos  # this base abs position
            # pileup_base.n  # this base depth

            for pileup_read in pileup_base.pileups:
                ls_read_mismatch_pairs = get_align_mismatch_pairs(
                    align=pileup_read.alignment, ref_genome_dict=dt_ref
                )

                if ls_read_mismatch_pairs is None:
                    # None
                    continue
                else:
                    # Not None
                    # [] -> False
                    ls_this_base_mut = [
                        ls for ls in ls_read_mismatch_pairs if ls[0] == pileup_base.pos
                    ]
                    if ls_this_base_mut:
                        # True [[53025969, 137, 'T', 'G']]
                        ls_this_base_mut = ls_this_base_mut[0]
                        print(ls_this_base_mut)
                    # [][0]
                # print(pileup_read.alignment.__dir__())
        # ['__str__', '__init__', '__new__', 'alignment',
        # 'query_position', 'query_position_or_next', 'indel', 'level',
        # 'is_del', 'is_head', 'is_tail', 'is_refskip',
        # '__doc__', '__reduce__', '__setstate__', '__repr__', '__hash__',
        # '__getattribute__', '__setattr__', '__delattr__', '__lt__',
        # '__le__', '__eq__', '__ne__', '__gt__', '__ge__', '__reduce_ex__',
        # '__subclasshook__', '__init_subclass__', '__format__',
        # '__sizeof__', '__dir__', '__class__']
        # print(pileupread.__dir__())
        # print(
        #     get_align_mismatch_pairs(align=pileupread, ref_genome_dict=dt_ref)
        # )
        # if not pileupread.is_del and not pileupread.is_refskip:
        #     # query position is None if is_del or is_refskip is set.
        #     # 取出当前read的该位置的base是哪个碱基
        #     base_info = pileupread.alignment.query_sequence[
        #         pileupread.query_position
        #     ]

    # 提取本site所有的reads，返回迭代器
    # it_reads = samfile.fetch(
    #     contig=dt_line["chr_name"],  # chrome
    #     start=dt_line["chr_index"],  # start
    #     stop=dt_line["chr_index"] + 1,  # stop
    # )
    # # 检查是否位于list中
    # for index, read in enumerate(it_reads):
    #     if not read.is_paired:
    #         # 除去未配对reads
    #         continue

    #     if read.query_name in ls_read_id_maternal:
    #         # 如果匹配ls_read_id_maternal，就写入samfile_out_maternal
    #         samfile_out_maternal.write(read)
    #     elif read.query_name in ls_read_id_fetal:
    #         # 如果匹配ls_read_id_fetal，就写入samfile_out_fetal
    #         samfile_out_fetal.write(read)
    #     else:
    #         # 都不匹配continue扔掉这条reads
    #         continue
    # process_bar(index / totalsites, start_str="", end_str="100%", total_length=15)
    # break
