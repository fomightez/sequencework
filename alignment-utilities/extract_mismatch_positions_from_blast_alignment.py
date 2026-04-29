#!/usr/bin/env python
# extract_mismatch_positions_from_blast_alignment.py
__author__ = "based on fomightez/Wayne Decatur style (fomightez on GitHub)"
__license__ = "MIT"
__version__ = "0.1.0"

# extract_mismatch_positions_from_blast_alignment.py
#
# PURPOSE: Takes a BLAST pairwise alignment (text output, default format 0)
# of two sequences and extracts the positions of mismatches relative to the
# QUERY sequence coordinates. This is useful for identifying positions to
# avoid when designing primers and probes (e.g., via IDT PrimerQuest) so
# that they do not span sites that differ between two strains/variants.
#
# Assumes a single pairwise BLAST alignment in standard text format (the
# kind you get from NCBI BLAST or command-line blastn with default output).
# No gaps/insertions/deletions are assumed (sequences are 98%+ identical
# over the region of interest with no indels).
#
# Very useful when you are trying to design qPCR primer and probes for two 
# strains using IDT's Primer Quest Tool, and so want to target where there is 
# identity.
#
# OUTPUT: Prints mismatch positions to stdout (1-based, relative to query).
#         Optionally, saves them to a text file.
#
# USAGE:
#   python extract_mismatch_positions_from_blast_alignment.py <blast_output.txt>
#   python extract_mismatch_positions_from_blast_alignment.py <blast_output.txt> --output mismatches.txt
#   python extract_mismatch_positions_from_blast_alignment.py <blast_output.txt> --show_context
#
# EXAMPLE blast output block this script expects (standard BLAST text format):
#
# Query  1    ATGCGTACGTAGCTAGCTAGC  21
#             ||||||||||||||||| |||
# Sbjct  1    ATGCGTACGTAGCTAGCAAGC  21
#
# The middle line uses '|' for matches and ' ' (space) for mismatches.
#
# DEPENDENCIES: none beyond Python standard library (re, argparse, sys, os)
#
# RELATED SCRIPTS in fomightez/sequencework alignment-utilities:
#   score_differences_between_sequences_by_pairwise_alignment.py
#   roughly_score_relationships_to_subject_seq_pairwise_premsa.py

# ADDITIONAL RESOURCE, a one-liner to take the individual positions and make 
# them ranges for  use in the exluded region in the Primer Quest tool, with 
# example positions = [57, 93, 149, 174, 226, 295, 323, 331, 352, 463, 467, 531, 553, 575, 576, 607, 631, 718, 862, 875, 889, 909, 922, 944, 997, 1108]
#print(",".join(f"{p}-{p}" for p in positions))

import re
import sys
import os
import argparse


###---------------------------------------------------------------------------###
#                           USER-ADJUSTABLE SETTINGS                           #
###---------------------------------------------------------------------------###

default_output_file_name = "mismatch_positions.txt"

###---------------------------------------------------------------------------###


def extract_mismatches(blast_text, show_context=False):
    """
    Parse a BLAST pairwise alignment text and return a list of mismatch
    positions (1-based) in query coordinates.

    Parameters
    ----------
    blast_text : str
        Full text of a BLAST output file (format 0, default pairwise).
    show_context : bool
        If True, also print the local nucleotide context around each mismatch.

    Returns
    -------
    list of int
        1-based positions in query sequence where mismatches occur.
    """
    # Collect all alignment blocks. Each block has three lines:
    #   Query  <start>  <seq>  <end>
    #   <spaces><match_string>
    #   Sbjct  <start>  <seq>  <end>
    #
    # We use a regex that captures:
    #   group 1 = query start
    #   group 2 = query sequence
    #   group 3 = match line (same length as sequences)
    #   group 4 = subject sequence

    block_pattern = re.compile(
        r'Query\s+(\d+)\s+([A-Za-z\-]+)\s+\d+\s*\n'   # Query line
        r'[ \t]+([|\s]+)\s*\n'                          # match/mismatch line
        r'Sbjct\s+\d+\s+([A-Za-z\-]+)\s+\d+',          # Sbjct line
        re.MULTILINE
    )

    mismatches = []

    for match in block_pattern.finditer(blast_text):
        query_start = int(match.group(1))
        query_seq   = match.group(2)
        match_line  = match.group(3)
        sbjct_seq   = match.group(4)

        # The match line may be shorter than seq lines if trailing spaces
        # were stripped; pad it out to be safe.
        match_line = match_line.ljust(len(query_seq))

        for i, symbol in enumerate(match_line):
            if i >= len(query_seq):
                break
            if symbol != '|':
                # This is a mismatch (space) or gap (-); report if not a gap
                q_base = query_seq[i]
                s_base = sbjct_seq[i] if i < len(sbjct_seq) else '?'
                if q_base != '-' and s_base != '-':
                    pos = query_start + i  # 1-based because query_start is 1-based
                    mismatches.append((pos, q_base, s_base))

    return mismatches


def main():
    parser = argparse.ArgumentParser(
        prog='extract_mismatch_positions_from_blast_alignment.py',
        description=(
            "extract_mismatch_positions_from_blast_alignment.py\n"
            "Takes a BLAST pairwise alignment text file (standard format 0)\n"
            "and reports the positions of mismatches in query coordinates\n"
            "(1-based). Useful for identifying sites to exclude when designing\n"
            "qPCR primers and probes across two similar sequences (e.g., viral\n"
            "strains that are ~98%% identical with no indels).\n\n"
            "Written to complement scripts in:\n"
            "  https://github.com/fomightez/sequencework/tree/master/alignment-utilities"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument(
        "blast_file",
        help="Path to the BLAST output text file (default format 0, pairwise).",
        metavar="BLAST_OUTPUT_FILE"
    )

    parser.add_argument(
        '-o', '--output',
        action='store',
        type=str,
        default=None,
        help=(
            "OPTIONAL: Name of file to save mismatch positions to. "
            "If not provided, results are printed to stdout only. "
            "Example: `--output mismatches.txt`"
        )
    )

    parser.add_argument(
        '--show_context',
        action='store_true',
        default=False,
        help=(
            "OPTIONAL: Show the query and subject bases at each mismatch "
            "position in addition to the position number."
        )
    )

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    # Verify input file exists
    if not os.path.isfile(args.blast_file):
        sys.stderr.write(
            "\nERROR: Cannot find the specified BLAST output file: "
            f"'{args.blast_file}'\nPlease check the path and try again.\n\n"
        )
        sys.exit(1)

    # Read input
    with open(args.blast_file, 'r') as f:
        blast_text = f.read()

    # Extract mismatches
    mismatches = extract_mismatches(blast_text, show_context=args.show_context)

    if not mismatches:
        sys.stderr.write(
            "\nNo mismatches found. Check that the file is standard BLAST "
            "pairwise format (format 0) and contains an alignment.\n\n"
        )
        sys.exit(0)

    # Build output text
    lines = []
    lines.append(
        f"\n{len(mismatches)} mismatch(es) found "
        f"(positions are 1-based, in query coordinates):\n"
    )

    if args.show_context:
        lines.append(f"{'Position':>10}  {'Query_base':>10}  {'Sbjct_base':>10}")
        lines.append(f"{'--------':>10}  {'----------':>10}  {'----------':>10}")
        for pos, q_base, s_base in mismatches:
            lines.append(f"{pos:>10}  {q_base:>10}  {s_base:>10}")
    else:
        lines.append("Position (1-based, query):")
        for pos, _, _ in mismatches:
            lines.append(f"  {pos}")

    # Also produce a compact comma-separated list useful for copy-paste
    pos_list = [str(p) for p, _, _ in mismatches]
    lines.append(f"\nComma-separated list of mismatch positions:\n{', '.join(pos_list)}\n")

    output_text = "\n".join(lines)

    # Print to stdout
    print(output_text)

    # Optionally save to file
    if args.output:
        with open(args.output, 'w') as out_f:
            out_f.write(output_text)
        sys.stderr.write(f"\nMismatch positions also saved to: '{args.output}'\n\n")


if __name__ == '__main__':
    main()
