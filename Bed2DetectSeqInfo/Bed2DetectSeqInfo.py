import logging

import pysam

from Bed2DetectSeqInfo.loader import (
    get_align_mismatch_pairs,
    load_reference_fasta_as_dict,
)

# pt_ref = "/Users/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/genome_ucsc_hg38.fa"
pt_ref = "/Users/zhaohuanan/zhaohn_HD/2.database/db_genomes/genome_fa/genome_ucsc_hg38/chr20.fa"
pt_bam = "/Users/zhaohuanan/zhaohn_HD/3.project/2021_DdCBE_topic/20210224_DetectSeq_all_bams/bam/293T-DdCBE-ND6-All-PD_rep1_hg38.MAPQ20.bam"
bam_file = pysam.AlignmentFile(pt_bam, "rb")
region_start = 53025819
region_stop = 53025866 + 1


if __name__ == "__main__":
    # load ref
    dt_ref = load_reference_fasta_as_dict(
        ref_fasta_path=pt_ref, ref_name=["chr20"], log_verbose=logging.DEBUG
    )
    # load bam
    bam_file = pysam.AlignmentFile(pt_bam, 'rb')
    for index_pileupcolumn, pileupcolumn in enumerate(
        bam_file.pileup(
            contig="chr20",  # chrome
            start=region_start,  # start
            stop=region_stop,  # stop
        )
    ):
        # restrict the pileupcolumn within targeted region
        if pileupcolumn.pos < region_start:
            continue
        elif pileupcolumn.pos > region_stop:
            break

        print(pileupcolumn.pileups.__dir__())
        # print(index_pileupcolumn, pileupcolumn.pos)
        if pileupcolumn.pos == 53025849:
            # 当pileup的当前columns就是这个base的时候，再进行下面的操作
            # 打印 该base下的depth信息
            print("\ncoverage at base %s = %s" % (pileupcolumn.pos, pileupcolumn.n))
            for pileupread in pileupcolumn.pileups:
                # ['__str__', '__init__', '__new__', 'alignment',
                # 'query_position', 'query_position_or_next', 'indel', 'level',
                # 'is_del', 'is_head', 'is_tail', 'is_refskip',
                # '__doc__', '__reduce__', '__setstate__', '__repr__', '__hash__',
                # '__getattribute__', '__setattr__', '__delattr__', '__lt__',
                # '__le__', '__eq__', '__ne__', '__gt__', '__ge__', '__reduce_ex__',
                # '__subclasshook__', '__init_subclass__', '__format__',
                # '__sizeof__', '__dir__', '__class__']
                print(pileupread.__dir__())
                print(
                    get_align_mismatch_pairs(align=pileupread, ref_genome_dict=dt_ref)
                )
                if not pileupread.is_del and not pileupread.is_refskip:
                    # query position is None if is_del or is_refskip is set.
                    # 取出当前read的该位置的base是哪个碱基
                    base_info = pileupread.alignment.query_sequence[
                        pileupread.query_position
                    ]


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
