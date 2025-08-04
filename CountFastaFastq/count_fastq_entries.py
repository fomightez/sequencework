#!/usr/bin/env python3
# This is meant to use with `uv` to run. 
# First install `uv` with `pip install uv` then run `!uv run {script_url} my_fasta.fa`, substituting your FASTA-formatted file for `my_fasta.fa`, and `script_url` already defined prior.
#-------------------------------------------------------------#
# Wayne's script `count_fastq_entries.py`.
#-------------------------------------------------------------#
# /// script
# requires-python = ">=3.12"
# dependencies = [
#   "rich",
# ]
# ///

import sys
try:
    fastq_file = sys.argv[1]
except IndexError:
    import rich
    rich.print("\n[bold red]I suspect you forgot to specify the file to read?[/bold red]\n Perhaps try something like the following:\n [bold black]!uv run count_fastq_entries.py my_fasta.fa[/bold black]\n[bold red]**EXITING !!**[/bold red]\n"); sys.exit(1)
def count_fastq_entries(fastq_file_path):
    """
    Count FASTQ entries.
    Blank lines shouldn't cause an issue or miscount, so is 
    insensitive to blank lines
    """
    reads_counted = 0
    with open(fastq_file_path, 'r') as infile:
        while True:
            # Read FASTQ record (4 lines) - exactly like your code
            header = infile.readline().strip()
            if not header:  # End of file
                break
            
            sequence = infile.readline().strip()
            plus_line = infile.readline().strip()
            quality = infile.readline().strip()
            
            # Just count instead of filtering
            reads_counted += 1
    
    print(f"Total reads in {fastq_file_path}: {reads_counted:,}")
    return reads_counted
count = count_fastq_entries(fastq_file)
#print(count)


'''
Other count fastq entries resources and code:

BASH and AWK commands:
bash_commands = """
# Bash methods to count FASTQ entries:

# Method 1: Count lines and divide by 4 (fastest)
wc -l file.fastq | awk '{print $1/4}'

# Method 2: Count header lines (more robust)
grep -c "^@" file.fastq

# Method 3: Count using awk every 4th line
awk 'NR%4==1{count++} END{print count}' file.fastq

# Method 4: Using bioawk (if available)
bioawk -c fastx 'END{print NR}' file.fastq

# Note shown here but biopython iterates on fastq, just as well as fasta, so could
# use biopython, too.


def count_fastq_entries_method1(fastq_file_path):
    """
    Count FASTQ entries by reading in 4-line blocks (similar to your existing code).
    This is the most robust method as it follows FASTQ structure exactly.
    """
    count = 0
    with open(fastq_file_path, 'r') as infile:
        while True:
            # Read FASTQ record (4 lines)
            header = infile.readline()
            if not header:  # End of file
                break
            
            sequence = infile.readline()
            plus_line = infile.readline()
            quality = infile.readline()
            
            # Verify this is a valid FASTQ record
            if header.startswith('@') and plus_line.startswith('+'):
                count += 1
            else:
                print(f"Warning: Invalid FASTQ record at entry {count + 1}")
                break
    
    return count

def count_fastq_entries_method2(fastq_file_path):
    """
    Count FASTQ entries by counting header lines that start with '@'.
    Faster but assumes well-formed FASTQ files.
    """
    count = 0
    with open(fastq_file_path, 'r') as infile:
        for line_num, line in enumerate(infile, 1):
            # Only count every 4th line starting from line 1, and verify it's a header
            if (line_num - 1) % 4 == 0:  # Header lines are at positions 1, 5, 9, etc.
                if line.startswith('@'):
                    count += 1
                else:
                    print(f"Warning: Expected header at line {line_num}, got: {line.strip()[:50]}")
    
    return count

def count_fastq_entries_method3(fastq_file_path):
    """
    Count FASTQ entries using itertools grouper pattern.
    Memory efficient and handles the 4-line structure explicitly.
    """
    from itertools import islice
    
    count = 0
    with open(fastq_file_path, 'r') as infile:
        while True:
            # Read next 4 lines
            four_lines = list(islice(infile, 4))
            if len(four_lines) < 4:
                break
            
            header, sequence, plus_line, quality = four_lines
            
            # Verify FASTQ format
            if header.startswith('@') and plus_line.startswith('+'):
                count += 1
            else:
                print(f"Warning: Invalid FASTQ record at entry {count + 1}")
                break
    
    return count

def count_fastq_entries_bash_style(fastq_file_path):
    """
    Python implementation mimicking bash approach: count lines and divide by 4.
    Fastest but assumes perfect 4-line structure with no blank lines.
    """
    with open(fastq_file_path, 'r') as infile:
        line_count = sum(1 for line in infile)
    
    if line_count % 4 != 0:
        print(f"Warning: Line count ({line_count}) is not divisible by 4")
        print("This might indicate formatting issues in the FASTQ file")
    
    return line_count // 4

# Test all methods and compare results
def test_all_methods(fastq_file_path):
    """
    Test all counting methods and compare results for validation.
    """
    print(f"Counting entries in: {fastq_file_path}")
    print("-" * 50)
    
    methods = [
        ("4-line block method (most robust)", count_fastq_entries_method1),
        ("Header counting method", count_fastq_entries_method2),
        ("Itertools grouper method", count_fastq_entries_method3),
        ("Line count / 4 method (fastest)", count_fastq_entries_bash_style)
    ]
    
    results = {}
    for name, method in methods:
        try:
            count = method(fastq_file_path)
            results[name] = count
            print(f"{name}: {count:,} entries")
        except Exception as e:
            print(f"{name}: Error - {e}")
    
    # Check if all methods agree
    counts = list(results.values())
    if len(set(counts)) == 1:
        print(f"\n✓ All methods agree: {counts[0]:,} entries")
    else:
        print(f"\n⚠ Methods disagree: {dict(zip(results.keys(), counts))}")
    
    return results

# Example usage matching your existing code style
def count_entries_like_your_code(fastq_file_path):
    """
    Count FASTQ entries using the same approach as your filtering code.
    This matches your existing code structure most closely.
    """
    reads_counted = 0
    with open(fastq_file_path, 'r') as infile:
        while True:
            # Read FASTQ record (4 lines) - exactly like your code
            header = infile.readline().strip()
            if not header:  # End of file
                break
            
            sequence = infile.readline().strip()
            plus_line = infile.readline().strip()
            quality = infile.readline().strip()
            
            # Just count instead of filtering
            reads_counted += 1
    
    print(f"Total reads in {fastq_file_path}: {reads_counted:,}")
    return reads_counted



'''
