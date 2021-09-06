import logging
import sys

import pandas as pd
from Bio import SeqIO


def load_reference_fasta_as_dict(
    ref_fasta_path, ref_name="All", log_verbose=logging.DEBUG
) -> dict:
    """Load the reference fasta file and return a dict of reference.

    Args:
        ref_fasta_path (str):
            Reference fasta file path
        ref_name (str/list, optional):
            If set All, load all "seq ids" in reference,
            else only try to load specific chromosome name, such as "chr1". Defaults to "All".
        log_verbose (object, optional): a logging states object, which can be logging.DEBUG,
            logging.INFO, logging.WARNING, logging.ERROR, logging.CRITICAL. Defaults to logging.DEBUG.

    Raises:
        IOError:
            Load file error

    Returns:
        dict/None:
            Normally, this function returns a dict, the key is the right "seq id" in the ref_name list and the value is sequences with .upper();
            it returns a empty_dict when there is no right "seq id" in the ref_name list;
            the loss of the key_value_pair can coursed by the wrong "seq id";
            it returns None when error occurs.
    """
    # log setting
    logging.basicConfig(
        level=log_verbose,
        format="%(levelname)-5s @ %(asctime)s: %(message)s ",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
        filemode="w",
        force=True,
    )
    # load genome as dict
    try:
        genome_fa = SeqIO.parse(handle=ref_fasta_path, format="fasta")
    except OSError:
        raise OSError("Wrong file path: %s" % ref_fasta_path)

    # init var
    ref_seq_dict = {}
    ref_name_set = set(ref_name)

    logging.info("Starting to load the reference genome...")

    for ref in genome_fa:
        if ref_name == "All":
            ref_seq_dict[ref.id] = ref.seq.upper()
            logging.debug("Loading genome...\t" + ref.id)

        elif ref.id in ref_name:
            ref_seq_dict[ref.id] = ref.seq.upper()
            logging.debug("Loading genome...\t" + ref.id)

            # remove already loaded seq
            ref_name_set.remove(ref.id)

            # load all info
            if len(ref_name_set) == 0:
                break

    logging.info("Loading genome done!")

    return ref_seq_dict


def back_indel_shift(info_index_list, cur_index) -> int:
    """return acc_shift(back_indel_shift)

    Args:
        info_index_list (list/tuple): list or tuples generated from align.cigar tuples
        cur_index (int): index related to MD tag in BAM file

    Returns:
        int: acc_shift index
    """
    # parse soft clip and insertion
    if len(info_index_list) == 0:
        return 0

    acc_shift = 0
    for info_start_index, info_len in info_index_list:
        if info_start_index >= cur_index:
            return acc_shift
        else:
            acc_shift += info_len
    return acc_shift


def get_align_mismatch_pairs(align, ref_genome_dict=None) -> list:
    """input a pysam AlignedSegment object

    Args:
        align (pysam.AlignedSeqment object): pysam.AlignedSeqment object
        ref_genome_dict (dict, optional): returned dict from load_reference_fasta_as_dict(). Defaults to None.

    Returns:
        list/None:
            it returns mismatch_pair_list, just like [ref_index, align_index, ref_base, align_base];
            and the "ref_index" is the same coordinate with UCSC genome browser;
            When NM == 0, it returns None.
    """
    # No mismatch
    try:
        if align.get_tag("NM") == 0:
            return None
    except:
        return None

    MD_tag_state = align.has_tag("MD")

    if MD_tag_state:
        # parse softclip, insertion and deletion
        info_index_list = []
        accu_index = 0

        for cigar_type, cigar_len in align.cigartuples:
            if cigar_type == 1 or cigar_type == 4:
                info_index_list.append((accu_index + 1, cigar_len))

            elif cigar_type == 2:
                info_index_list.append((accu_index + 1, -cigar_len))

            accu_index += cigar_len

        # parse MD tag
        mismatch_pair_list = []
        cur_base = ""
        cur_index = 0
        bases = align.get_tag("MD")

        i = 0
        while i < len(bases):
            base = bases[i]

            if base.isdigit():
                cur_base += base
                i += 1

            else:
                cur_index += int(cur_base)
                cur_base = ""

                if base == "^":
                    i += 1
                    del_str = ""

                    while (bases[i].isalpha()) and (i < len(bases)):
                        del_str += bases[i]
                        i += 1

                    cur_index += len(del_str)
                    del_str = ""

                elif base.isalpha():
                    cur_index += 1
                    ref_base = base
                    i += 1

                    # add into list
                    fix_index = cur_index + back_indel_shift(info_index_list, cur_index)

                    if fix_index < len(align.query_sequence):
                        mismatch_pair_list.append(
                            [
                                cur_index + align.reference_start,
                                cur_index - 1,
                                ref_base,
                                align.query_sequence[fix_index - 1],
                            ]
                        )
                    else:
                        return None

        return mismatch_pair_list
    else:
        mismatch_pair_list = []
        for align_idx, ref_idx in align.get_aligned_pairs():
            if (align_idx is not None) and (ref_idx is not None):
                align_base = align.query_sequence[align_idx]
                ref_base = ref_genome_dict[align.reference_name][ref_idx]

                if align_base != ref_base:
                    mismatch_pair_list.append(
                        [ref_idx + 1, align_idx, ref_base, align_base]
                    )

        return mismatch_pair_list


def aligned_segment_get_reference_sequence(align, ref_genome_dict):
    """[summary]

    Args:
        align (pysam.AlignedSeqment object): pysam.AlignedSeqment object
        ref_genome_dict (dict, optional): returned dict from load_reference_fasta_as_dict().
    Returns:
        Bio.Seq object: Bio Seq object
    """
    MD_tag_state = align.has_tag("MD")

    # align_start_idx = align.get_aligned_pairs()[0][0]
    ref_start_idx = align.get_aligned_pairs()[0][1]
    # align_end_idx = align.get_aligned_pairs()[-1][0]
    ref_end_idx = align.get_aligned_pairs()[-1][1]

    if MD_tag_state:
        return align.get_reference_sequence().__str__()
    else:
        return ref_genome_dict[align.reference_name][ref_start_idx : ref_end_idx + 1]


def get_absolute_position(bowtie_table: str) -> pd.DataFrame:
    """get region position from bowtie1 mapping standard output file

    Args:
        bowtie_table (str): path of standard output file (must use the same reference with load_reference_fasta_as_dict())

    Returns:
        pd.DataFrame: the absolution positions of aimed regions.
        cols: [
            "chrom",  # chromosome
            "start",  # start (contain)
            "stop",  # stop (contain)
            "read_id",  # region name
            "score",  # fake score
            "strand",  # strand
            "sequence",  # sequence (ref strand +)
            ]
    """
    df = pd.read_csv(
        bowtie_table,
        sep="\t",
        names=[
            "read_id",
            "strand",
            "chrom",
            "start",
            "sequence",
            "align_info",
            "score",
        ],
        index_col=False,
    )
    df.insert(4, "stop", df["start"] + df.sequence.map(len), allow_duplicates=False)
    return df[["chrom", "start", "stop", "read_id", "score", "strand", "sequence"]]


# unit test
# if __name__ == "__main__":
