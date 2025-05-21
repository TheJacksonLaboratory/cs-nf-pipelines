"""
Script to generate random genomic positions while excluding specified regions.

This script reads a genome file (e.g., FASTA.fai) and two BED files containing 
excluded regions. It generates a specified number of random genomic positions 
that do not overlap with the excluded regions and writes the positions to an 
output BED file.

This random file is used in NanoSV for coverage calculation.

Functions:
    read_bed(bed_file):
        Reads a BED file and returns a dictionary of excluded regions per chromosome.

    get_chromosome_lengths(genome_file):
        Reads a genome file and returns a dictionary of chromosome lengths.

    is_excluded(chrom, pos, excluded_regions):
        Checks if a given genomic position is within any excluded regions.

    generate_random_positions(genome_file, bed1_file, bed2_file, num_positions=1000000):
        Generates random genomic positions excluding regions from two BED files.

    write_positions_to_bed(positions, output_file):

Usage:
    python random_position.py --genome_file <genome_file> --bed1_file <bed1_file> 
                              --bed2_file <bed2_file> --output_file <output_file> 
                              [--num_positions <num_positions>]

Example:
    python random_position.py --genome_file Homo_sapiens_assembly38.fasta.fai 
                              --bed1_file simpleRepeat.bed --bed2_file gap.bed 
                              --output_file GRCh38_random_1000000.bed 
                              --num_positions 1000000
"""

import random
import bisect
from collections import defaultdict
import argparse

def read_bed(bed_file):
    """
    Reads a BED file and returns a dictionary where keys are chromosome names
    and values are lists of (start, end) tuples representing the excluded regions.
    Assumes the BED file is sorted.

    Args:
        bed_file (str): Path to the BED file.

    Returns:
        dict: A dictionary of excluded regions per chromosome.
    """
    excluded_regions = defaultdict(list)
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) > 2:  # Ensure the line has enough columns
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    excluded_regions[chrom].append((start, end))
                # Optionally, add a warning for malformed BED lines:
                elif len(parts) > 0:
                    print(f"Warning: Skipping malformed BED line: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: BED file not found: {bed_file}")
        sys.exit(1)  # Exit the program with an error code

    return excluded_regions

def get_chromosome_lengths(genome_file):
    """
    Reads a genome file (e.g., FASTA.fai) and returns a dictionary of
    chromosome lengths.

    Args:
        genome_file (str): Path to the genome file.

    Returns:
        dict: A dictionary of chromosome lengths.  Returns empty dict on error.
    """
    chrom_lengths = {}
    try:
        with open(genome_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    chrom = parts[0]
                    length = int(parts[1])
                    chrom_lengths[chrom] = length
                elif len(parts) > 0:
                    print(f"Warning: Skipping malformed genome file line: {line.strip()}")
    except FileNotFoundError:
        print(f"Error: Genome file not found: {genome_file}")
        sys.exit(1)  # Exit the program with an error code

    return chrom_lengths

def is_excluded(chrom, pos, excluded_regions):
    """
    Checks if a given genomic position is within any of the excluded regions
    for the specified chromosome.

    Args:
        chrom (str): The chromosome name.
        pos (int): The genomic position.
        excluded_regions (dict): A dictionary of excluded regions.

    Returns:
        bool: True if the position is excluded, False otherwise.
    """
    if chrom not in excluded_regions:
        return False  # Chromosome not in excluded regions, so not excluded

    regions = excluded_regions[chrom]
    # Use binary search to efficiently check if the position is in any excluded region
    i = bisect.bisect_left(regions, (pos, pos))  # Find the insertion point
    if i > 0 and regions[i-1][1] >= pos:
        return True #check overlap with the previous region
    if i < len(regions) and regions[i][0] <= pos:
        return True
    return False

def generate_random_positions(genome_file, bed1_file, bed2_file, num_positions=1000000):
    """
    Generates a specified number of random genomic positions, excluding
    regions from two BED files, and distributes them evenly across chromosomes.

    Args:
        genome_file (str): Path to the genome file (e.g., FASTA.fai).
        bed1_file (str): Path to the first BED file.
        bed2_file (str): Path to the second BED file.
        num_positions (int, optional): The number of random positions to generate.
            Defaults to 1,000,000.

    Returns:
        list: A list of (chrom, pos) tuples representing the random genomic positions.
              Returns an empty list on error.
    """

    chrom_lengths = get_chromosome_lengths(genome_file)
    if not chrom_lengths:
        print("Error: Could not retrieve chromosome lengths.  Check genome file.")
        return []

    excluded_regions1 = read_bed(bed1_file)
    excluded_regions2 = read_bed(bed2_file)
    # Combine the excluded regions from both BED files into a single dictionary
    all_excluded_regions = defaultdict(list)
    for chrom, regions in excluded_regions1.items():
        all_excluded_regions[chrom].extend(regions)
    for chrom, regions in excluded_regions2.items():
        all_excluded_regions[chrom].extend(regions)
    # Sort the regions for each chromosome to enable efficient searching
    for chrom in all_excluded_regions:
        all_excluded_regions[chrom].sort()

    random_positions = []
    positions_per_chrom = {chrom: 0 for chrom in chrom_lengths}
    target_positions_per_chrom = {
        chrom: num_positions // len(chrom_lengths) for chrom in chrom_lengths
    }
    # Distribute the remaining positions
    remaining_positions = num_positions % len(chrom_lengths)
    chroms = list(chrom_lengths.keys())
    for i in range(remaining_positions):
        target_positions_per_chrom[chroms[i]] += 1

    attempts = 0
    max_attempts = num_positions * 10  # Avoid infinite loop, allow 10x attempts.
    while len(random_positions) < num_positions and attempts < max_attempts:
        chrom = random.choice(chroms)
        chrom_length = chrom_lengths[chrom]
        pos = random.randint(1, chrom_length)
        attempts += 1

        if not is_excluded(chrom, pos, all_excluded_regions):
            random_positions.append((chrom, pos))
            positions_per_chrom[chrom] += 1
            if positions_per_chrom[chrom] >= target_positions_per_chrom[chrom]:
                del target_positions_per_chrom[chrom] #remove the chrom from the dict
                chroms = list(target_positions_per_chrom.keys())

    if attempts >= max_attempts:
        print(f"Warning: Reached maximum attempts ({max_attempts}).  Could not generate all requested positions.")

    return random_positions
def write_positions_to_bed(positions, output_file):
    """
    Writes the generated random positions to a BED file.

    Args:
        positions (list): A list of (chrom, pos) tuples.
        output_file (str): The path to the output BED file.
    """
    try:
        with open(output_file, 'w') as f:
            for chrom, pos in positions:
                # BED format is 0-based, so subtract 1 from pos.  However, we generated 1-based positions.
                f.write(f"{chrom}\t{pos-1}\t{pos}\n")
    except Exception as e:
        print(f"Error writing to BED file: {e}")

if __name__ == "__main__":
    # Example usage:
    parser = argparse.ArgumentParser(description="Generate random genomic positions excluding specified regions.")
    parser.add_argument("--genome_file", required=True, help="Path to the genome file (e.g., FASTA.fai).")
    parser.add_argument("--bed1_file", required=True, help="Path to the first BED file with excluded regions.")
    parser.add_argument("--bed2_file", required=True, help="Path to the second BED file with excluded regions.")
    parser.add_argument("--output_file", required=True, help="Path to the output BED file.")
    parser.add_argument("--num_positions", type=int, default=1000000, help="Number of random positions to generate (default: 1,000,000).")



    args = parser.parse_args()

    genome_file = args.genome_file
    bed1_file = args.bed1_file
    bed2_file = args.bed2_file
    output_file = args.output_file
    num_positions = args.num_positions

    random_positions = generate_random_positions(genome_file, bed1_file, bed2_file, num_positions=1000000)



    if random_positions:
        # Sort the random positions by chromosome and starting position
        random_positions.sort(key=lambda x: (x[0], x[1]))

        # Count the number of positions on each chromosome
        chrom_counts = defaultdict(int)
        for chrom, _ in random_positions:
            chrom_counts[chrom] += 1

        # Print the counts to the screen
        for chrom, count in sorted(chrom_counts.items()):
            print(f"{chrom}: {count}")
        print(f"Generated {len(random_positions)} random positions.")
        write_positions_to_bed(random_positions, output_file)
        print(f"Wrote positions to {output_file}")
    else:
        print("Error: No random positions generated.")


"""
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/simpleRepeat.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/database/gap.txt.gz
gunzip *.gz
cut -f2,3,4 simpleRepeat.txt > simpleRepeat.bed
cut -f2,3,4 gap.txt > gap.bed
python random_position.py --genome_file Homo_sapiens_assembly38.fasta.fai --bed1_file simpleRepeat.bed --bed2_file gap.bed --output_file GRCh38_random_1000000.bed
rm simpleRepeat*
rm gap*

wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/simpleRepeat.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm39/database/gap.txt.gz
gunzip *.gz
cut -f2,3,4 simpleRepeat.txt > simpleRepeat.bed
cut -f2,3,4 gap.txt > gap.bed
python random_position.py --genome_file Mus_musculus.GRCm38.dna.primary_assembly.fa.fai --bed1_file simpleRepeat.bed --bed2_file gap.bed --output_file GRCm38_random_1000000.bed
rm simpleRepeat*
rm gap*


wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/simpleRepeat.txt.gz
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz
gunzip *.gz
cut -f2,3,4 simpleRepeat.txt > simpleRepeat.bed
cut -f2,3,4 gap.txt > gap.bed
python /random_position.py --genome_file Mus_musculus.GRCm39.dna.primary_assembly.fa.fai --bed1_file simpleRepeat.bed --bed2_file gap.bed --output_file GRCm39_random_1000000.bed
rm simpleRepeat*
rm gap*
"""