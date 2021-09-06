import logging
import sys

import pandas as pd


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
        log_verbose (logging object, optional): a logging states object, which can be logging.DEBUG,
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


def get_reference_sequence_from_aligned_segment(align, ref_genome_dict) -> SeqIO:
    """get reference sequence of pysam.AlignedSeqment object

    Args:
        align (pysam.AlignedSeqment object):
            pysam.AlignedSeqment object
        ref_genome_dict (dict, optional):
            returned dict from load_reference_fasta_as_dict().
    Returns:
        Bio.Seq object:
            it returns reference sequence (Bio Seq object)
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


def load_region_absolute_position(bowtie_table: str) -> pd.DataFrame:
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
