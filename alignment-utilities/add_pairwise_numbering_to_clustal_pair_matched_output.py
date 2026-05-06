#!/usr/bin/env python
# add_pairwise_numbering_to_clustal_pair_matched_output.py
__author__ = "Wayne Decatur" #fomightez on GitHub
__license__ = "MIT"
__version__ = "0.1.0"

# add_pairwise_numbering_to_clustal_pair_matched_output.py
# Based on add_actual_position_numbering_to_sequence_in_top_line_of_clustal_output.py
# by Wayne Decatur (fomightez on GitHub)
#
# PURPOSE: Takes a text document of a PAIRWISE alignment in CLUSTAL format and:
#   1. Adds position numbering above the top (first) sequence, as in the 
#      original script.
#   2. Moves the conservation symbol line (****) to between the two sequences.
#   3. Adds contiguous position numbering below the bottom (second) sequence.
#
# Input format expected (standard CLUSTAL pairwise, 2 sequences per block):
#
# SEQ1_ID   CAGTGGAAATGAAAGTGATGGGG...
# SEQ2_ID   CAGTGGAAATGAAAGTGATGGGG...
#            ************************************************************
#
# Output format:
#
#            1                                                         60
# SEQ1_ID   CAGTGGAAATGAAAGTGATGGGG...
#            ************************************************************
# SEQ2_ID   CAGTGGAAATGAAAGTGATGGGG...
#            1                                                         60
#
# TO RUN:
# python add_pairwise_numbering_to_clustal_pair_matched_output.py ALIGNMENT_TEXT_FILE
#
# Or in a Jupyter notebook after importing:
# add_pairwise_numbering_to_clustal_pair_matched_output("alignment.clustal")

import sys
import os

output_name = 'alignment_pairwise_numbered.clustal'
suffix_for_saving = "_PAIRWISE_NUMBERED"


###---------------------------HELPER FUNCTIONS---------------------------------###

def generate_output_file_name(file_name):
    main_part_of_name, file_extension = os.path.splitext(file_name)
    if '.' in file_name:
        return main_part_of_name + suffix_for_saving + file_extension
    else:
        return file_name + suffix_for_saving + ".clustal"


def make_number_line(seq, start_pos, id_width):
    """
    Given a sequence string (with possible gaps) and the starting residue
    position, return a formatted number line string like:
        '           1                                                         60'
    id_width is the number of characters used by the sequence identifier + 
    whitespace prefix so the numbers align with the sequence characters.
    """
    char_num_without_gaps = len(seq.replace("-", ""))
    end_pos = start_pos + char_num_without_gaps
    if char_num_without_gaps:
        current_position = start_pos + 1
    else:
        current_position = start_pos
    if end_pos < 0:
        end_pos = 0
    spacing = len(seq) - len(str(end_pos))
    if len(seq) > (len(str(current_position)) + len(str(end_pos))):
        num_text = '{:<{}}{}'.format(current_position, spacing, end_pos)
    else:
        if spacing < len(str(end_pos)):
            spacing = len(seq) + len(str(end_pos))
        num_text = '{:>{}}'.format(end_pos, spacing)
    num_line = (" " * id_width) + num_text
    return num_line, end_pos


###--------------------------END OF HELPER FUNCTIONS---------------------------###


def add_pairwise_numbering_to_clustal_pair_matched_output(
        alignment, output_name=output_name):
    '''
    Main function of script.
    Takes a pairwise CLUSTAL format alignment file (or string) and:
      - Adds position numbering above the first sequence
      - Moves the conservation symbol line to between the two sequences
      - Adds position numbering below the second sequence
    '''
    # Read input
    try:
        with open(alignment, 'r') as the_file:
            alignment = the_file.read()
    except (TypeError, OSError, IOError):
        pass  # alignment is already a string

    sys.stderr.write("Alignment read...")

    lines = alignment.split("\n")

    # Identify the two sequence identifiers from the alignment blocks
    # (skip the header line and blank lines)
    seq_ids = []
    for line in lines:
        line_stripped = line.strip()
        if not line_stripped:
            continue
        first_word = line_stripped.split()[0]
        # Skip CLUSTAL header line
        if first_word.upper().startswith("CLUSTAL"):
            continue
        # Skip conservation symbol lines (start with space or * or : or .)
        if line[0] == ' ' or first_word in ('*', ':', '.'):
            continue
        if first_word not in seq_ids:
            seq_ids.append(first_word)
        if len(seq_ids) == 2:
            break

    if len(seq_ids) < 2:
        sys.stderr.write(
            "\nERROR: Could not identify two sequence identifiers. "
            "Is this a pairwise alignment?\n")
        return

    id1, id2 = seq_ids
    sys.stderr.write(
        "\nSequence 1 (top): '{}'\nSequence 2 (bottom): '{}'...".format(
            id1, id2))

    # Process blocks
    # Standard CLUSTAL block order: seq1, seq2, conservation symbol line
    # We want to output:
    #   number line for seq1
    #   seq1 line
    #   conservation symbol line  (moved here, between sequences)
    #   seq2 line
    #   number line for seq2

    growing_output = []
    pos1 = 0  # running position counter for seq1
    pos2 = 0  # running position counter for seq2

    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()

        # Pass through empty lines and header
        if not stripped:
            growing_output.append(line)
            i += 1
            continue

        first_word = stripped.split()[0]

        if first_word.upper().startswith("CLUSTAL"):
            growing_output.append(line)
            i += 1
            continue

        # Detect a pre-existing number/annotation line above seq1 (starts with
        # spaces and contains digits) — skip it, we'll regenerate it
        if line[0] == ' ' and any(c.isdigit() for c in stripped) and \
                not all(c in '*.:-| ' for c in stripped):
            i += 1
            continue

        # Conservation symbol line in its original position (after seq2) —
        # skip here, we'll reinsert it between sequences during block processing
        if line[0] == ' ' and stripped and all(
                c in '*.:-| ' for c in stripped):
            i += 1
            continue

        # Seq1 line — triggers processing of a full block
        if first_word == id1:
            seq1_line = line
            seq1 = seq1_line.split(id1, 1)[1].strip()
            id1_width = seq1_line.index(seq1)  # width of id + whitespace

            # Standard CLUSTAL block order: seq1, seq2, conservation_line
            # Look ahead for seq2 first, then conservation line after it.
            conservation_line = None
            seq2_line = None
            j = i + 1
            # Find seq2
            while j < len(lines):
                next_line = lines[j]
                next_stripped = next_line.strip()
                if not next_stripped:
                    j += 1
                    continue
                next_first = next_stripped.split()[0]
                if next_first == id2:
                    seq2_line = next_line
                    j += 1
                    break
                break
            # Find conservation line after seq2
            while j < len(lines):
                next_line = lines[j]
                next_stripped = next_line.strip()
                if not next_stripped:
                    j += 1
                    continue
                if next_line[0] == ' ' and all(
                        c in '*.:-| ' for c in next_stripped):
                    conservation_line = next_line
                    j += 1
                    break
                break

            seq2 = seq2_line.split(id2, 1)[1].strip() if seq2_line else ""
            id2_width = seq2_line.index(seq2) if seq2_line else id1_width

            # Generate number line for seq1 (above)
            num_line1, new_pos1 = make_number_line(seq1, pos1, id1_width)
            growing_output.append(num_line1)
            growing_output.append(seq1_line.rstrip())

            # Conservation line moved to between sequences
            if conservation_line is not None:
                growing_output.append(conservation_line.rstrip())
            else:
                growing_output.append("")

            # Seq2 line followed by its number line (below)
            if seq2_line:
                growing_output.append(seq2_line.rstrip())
                num_line2, new_pos2 = make_number_line(seq2, pos2, id2_width)
                growing_output.append(num_line2)

            pos1 = new_pos1
            if seq2_line:
                pos2 = new_pos2

            i = j  # advance past all consumed lines
            continue

        # Any other line — pass through
        growing_output.append(line)
        i += 1

    output = "\n".join(growing_output)

    with open(output_name, 'w') as output_file:
        output_file.write(output)
    sys.stderr.write(
        "\nAnnotated alignment saved to '{}'.".format(output_name))


###--------------------------END OF MAIN FUNCTION----------------------------###


def main():
    kwargs = {}
    kwargs['output_name'] = output_name
    add_pairwise_numbering_to_clustal_pair_matched_output(alignment, **kwargs)


if __name__ == "__main__" and '__file__' in globals():
    import argparse
    parser = argparse.ArgumentParser(
        prog='add_pairwise_numbering_to_clustal_pair_matched_output.py',
        description="add_pairwise_numbering_to_clustal_pair_matched_output.py "
                    "takes a pairwise CLUSTAL format alignment and: "
                    "(1) adds position numbering above the top sequence, "
                    "(2) moves the conservation symbol line to between the "
                    "two sequences, and (3) adds position numbering below "
                    "the bottom sequence. "
                    "**** Based on script by Wayne Decatur (fomightez @ github) ***")

    parser.add_argument("align_file",
                        help="CLUSTAL format pairwise alignment file.",
                        metavar="ALIGNMENT_FILE")
    parser.add_argument('-o', '--output', action='store', type=str,
                        default=output_name,
                        help="OPTIONAL: output file name. Default inserts "
                             "'{}' before the file extension.".format(
                            suffix_for_saving))

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()
    alignment = args.align_file
    if args.output == output_name:
        output_name = generate_output_file_name(alignment)
    else:
        output_name = args.output

    main()
