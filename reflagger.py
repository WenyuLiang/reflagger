import sys
import pysam
import json
from multiprocessing import Pool
import numpy as np
from collections import defaultdict, namedtuple, Counter
import re
import time
import argparse
parser = argparse.ArgumentParser()
# Create a subparser object for the three modes
subparsers = parser.add_subparsers(dest='mode', required=True)

# Define the "reflag" mode as a subparser
reflag_parser = subparsers.add_parser("reflag", help="Flag to query")
reflag_parser.add_argument("-b", "--bam_file", help="Path to the BAM file", type=str, required=True)
reflag_parser.add_argument("-o", "--output", help="Path to the output file", type=str, required=True)
reflag_parser.add_argument("-v", "--vcf", help="Genomic position to query", type=str, required=True)

# reflag_parser.add_argument("-r", "--ref", help="Path to the reference FASTA file", type=str, required=True)

# Define the "alignment" mode as a subparser
alignment_parser = subparsers.add_parser("alignment", help="Get all alignments w.r.t. a read")
alignment_parser.add_argument("-b", "--bam_file", help="Path to the BAM file", type=str, required=True)
alignment_parser.add_argument("-s", "--seq_id", help="Sequence ID to query", type=str, required=True)

# Define the "pileup" mode as a subparser
pileup_parser = subparsers.add_parser("pileup", help="Get pileup w.r.t. a position")
pileup_parser.add_argument("-b", "--bam_file", help="Path to the BAM file", type=str, required=True)
pileup_parser.add_argument("-p", "--pos", help="Genomic position to query", type=int, required=True)
pileup_parser.add_argument("-r", "--ref", help="Path to the reference FASTA file", type=str, required=True)

# Define "filtering" mode as a subparser
filter_parser = subparsers.add_parser("filter", help="Filtering")
filter_parser.add_argument("-b", "--bam_file", help="Path to the BAM file", type=str, required=True)
filter_parser.add_argument("-o", "--output", help="Path to the output file", type=str, required=True)
filter_parser.add_argument("-m", "--mpileup", help="Filter for mpile", action='store_true', required=False)

args = parser.parse_args()

variant = namedtuple('variant', ('type', 'base'))
alignment = namedtuple('alignment', ('mapq', 'chr', 'start', 'end', 'length'))
stats = defaultdict(Counter)


def calculate_alignment_score(cigar, match_score=2, mismatch_penalty=4, gap_open_penalties=(4, 24), gap_extend_penalties=(2, 1)):
    # Regular expression pattern to extract operations and their lengths from the CIGAR string
    pattern = re.compile(r'(\d+)([MIDNSHP=X])')

    alignment_score = 0

    for length, operation in pattern.findall(cigar):
        length = int(length)

        if operation in ('M', '=', 'X'):  # Match or mismatch
            if operation == 'M' or operation == '=':
                alignment_score += match_score * length
            else:  # operation == 'X'
                alignment_score -= mismatch_penalty * length

        elif operation in ('I', 'D'):  # Insertion or deletion (gap)
                # Calculate gap open penalty
            gap_open_penalty = min(gap_open_penalties[0] + length * gap_extend_penalties[0],
                                    gap_open_penalties[1] + length * gap_extend_penalties[1])
            alignment_score -= gap_open_penalty

        # No need to account for other operations (S, H, N, P) in the score calculation

    return alignment_score

def get_flag(flag):
    if flag & 0x800:
        return 'supplementary'
    elif flag & 0x100:
        return 'secondary'
    elif flag & 0x4: #unmapped
        return 'unmapped'
    else:
        return 'primary'

def get_alignments(bam_file, seq_id):
    """
    Yields set of split reads for each contig separately. Only reads primary alignments
    and infers the split reads from SA alignment tag
    """

    if bam_file == sys.stdin:
        aln_file = pysam.AlignmentFile('-', "r")
    else:
        aln_file = pysam.AlignmentFile(bam_file, "rb")
    for aln in aln_file.fetch(multiple_iterators=True):        
        if aln.query_name==seq_id:
            if aln.has_tag('AS'):
                # Get the alignment score using the 'AS' tag
                score = aln.get_tag('AS')
                print(aln.query_name, get_flag(aln.flag), score, aln.mapping_quality, aln.reference_name, aln.reference_start+1, aln.reference_end + 1, sep="\t")
            else:
                exit("No alignment score found")

def get_info(aln):
    #return alignment(get_flag(aln.flag), aln.get_tag('AS'), aln.mapping_quality, aln.reference_name, aln.reference_start+1, aln.reference_end + 1, aln.reference_end - aln.reference_start, sep="\t")
    return alignment(aln.mapping_quality, aln.reference_name, aln.reference_start+1, aln.reference_end + 1, aln.query_alignment_length)

def get_pileup(bam_file, pos, ref): # bam_file is acrually a pysam.AlignmentFile object; ref is a pysam.FastaFile object
    chr, start_end = pos.split(":")
    start, end = int(start_end.split("-")[0]), int(start_end.split("-")[1])
    print(f'start: {start} end: {end}')
    for pileupcolumn in bam_file.pileup(chr, start, end, truncate=True, stepper='nofilter'):
        print(f'chr: {pileupcolumn.reference_name} pos: {pileupcolumn.reference_pos} base: {pileupcolumn.n}')
        for pileupread in pileupcolumn.pileups:
            if not pileupread.is_refskip:
                if pileupread.indel == 0:
                    try: # aligned base
                        #print('\t' + f'read id: {pileupread.alignment.query_name} flag: {get_flag(pileupread.alignment.flag)} base: {pileupread.alignment.query_sequence[pileupread.query_position]}')
                        #print('\t' + f'read id: {pileupread.alignment.query_name} flag: {get_flag(pileupread.alignment.flag)} base: {variant("A", pileupread.alignment.query_sequence[pileupread.query_position])}')
                        stats[get_flag(pileupread.alignment.flag)][variant("A", pileupread.alignment.query_sequence[pileupread.query_position])] += 1
                    except: # deleted base
                        continue
                        # ref_base = reference.fetch(pileupcolumn.reference_name, pileupcolumn.reference_pos, pileupcolumn.reference_pos + 1).upper()
                        # print('\t' + f'read id: {pileupread.alignment.query_name} flag: {get_flag(pileupread.alignment.flag)} base: {variant("D", ref_base)}')

                elif pileupread.indel > 0:  # Insertion
                    inserted_bases = pileupread.alignment.query_sequence[pileupread.query_position + 1:pileupread.query_position + 1 + pileupread.indel]
                    #print('\t' + f'read id: {pileupread.alignment.query_name} flag: {get_flag(pileupread.alignment.flag)} base: {variant("I", inserted_bases)}')
                    stats[get_flag(pileupread.alignment.flag)][variant("I", inserted_bases)] += 1
                else:  # Deletion
                    deleted_bases = ref.fetch(pileupcolumn.reference_name, pileupcolumn.reference_pos + 1, pileupcolumn.reference_pos + 1 - pileupread.indel).upper()
                    #print('\t' + f'read id: {pileupread.alignment.query_name} flag: {get_flag(pileupread.alignment.flag)} base: {variant("D", deleted_bases)}')
                    stats[get_flag(pileupread.alignment.flag)][variant("D", deleted_bases)] += 1
    for flag in stats:
        print(f'{flag}',stats[flag])

def get_indel_positions(alignment): 
    indel_positions = defaultdict(list)
    ref_pos = alignment.reference_start
    for op, length in alignment.cigartuples:
        if op in (1, 2):
            indel_positions[alignment.reference_name].append(ref_pos)
        if op != 1:  # Move reference position if not an insertion
            ref_pos += length
    return indel_positions

def score(allele_frequency):
    if allele_frequency > 0.5:
        return (allele_frequency-0.5)*2
    else:
        return 8*(allele_frequency-0.5)**3
def indicator(x):
    return 1 if x > 0 else -1

def process_vcf(vcf_file):
    all_indel_positions = defaultdict(lambda: defaultdict(int))
    with open(vcf_file, 'r') as f:
        for line in f:
            ref= line.split()[0]
            pos = int(line.split()[1])
            all_indel_positions[ref][pos] = 0
    return all_indel_positions   
    
def get_all_consensus_indel_positions_parallel(pile, threshould=8):
    consensus_dict = defaultdict(lambda: defaultdict(dict))
    with open(pile) as p:
        for line in p:
            dep = int(line.split()[3])
            allele = line.split()[4].upper()
            if dep <= threshould:
                continue
            del_n = allele.count('-') / dep
            ins_n = allele.count('+') / dep
            #res_n = 1 - del_n - ins_n
            max_n = max(del_n, ins_n) 
            # if max_n == 1: # all the reads are the same
            #     consensus_dict[ref][pos-1] = 1
            #if 0.3 < max_n < 0.6: # limit the depth and max allele frequency  
            if max_n > 0.2:
                pos = int(line.split()[1])
                ref = line.split()[0]
                consensus_dict[ref][pos-1]['-'] = score(allele.count('-') / dep)
                consensus_dict[ref][pos-1]['+'] = score(allele.count('+') / dep)
                if (1 - allele.count('-') / dep - allele.count('+') / dep)<0.3:
                    consensus_dict[ref][pos-1]['o'] = score(1 - allele.count('-') / dep - allele.count('+')/dep) # other should get penalty if not main allele
                else:
                    consensus_dict[ref][pos-1]['o'] = 0
    return consensus_dict


def filter_secondary(bam, out, mpile): # filter out secondary alignments whose mapq is less than 10
    with pysam.AlignmentFile(bam, 'rb') as f, pysam.AlignmentFile(out, 'wb', template=f) as o:
        for aln in f:
            #if aln.is_supplementary and aln.mapping_quality <= mapq:
            if aln.is_supplementary: 
                continue
            elif aln.is_secondary and mpile: #reflag to primary
                aln.is_secondary = False
                aln.flag &= ~0x800
                o.write(aln)    
            elif aln.is_secondary and not mpile:
                o.write(aln)
            else:
                o.write(aln)

def calculate_consistency(aln, consensus_dict):
    # if aln.query_name != '9390f672-ec39-4a92-9f58-47fdd6e4207a':
    #     return 0
    ref_dict = consensus_dict.get(aln.reference_name, {})
    score = 0
    pos = aln.pos
    c = 0 # normalizer
    for op, length in aln.cigar:
        if op in {0, 7, 8}:
            if pos-1 in ref_dict:
                c += 1
                score += ref_dict[pos-1]['o']
                # to_d = tuple([get_flag(aln.flag),pos, ref_dict[pos-1]['o']])
                # if to_d[2] != 0: # only record the position with penalty
                #     my_dict[get_info(aln)].append(to_d)
            pos += length       
        if op in {1, 2}:
            indel_type = '+' if op == 1 else '-' # insertion or deletion                    
            score += ref_dict.get(pos-1, {}).get(indel_type, 0)
            c += 1
            # to_d = tuple([get_flag(aln.flag),pos, ref_dict.get(pos-1, {}).get(indel_type, 0)])
            # if to_d[2] != 0:
            #     my_dict[get_info(aln)].append(to_d)
            # my_dict[get_info(aln)].append()
            pos += length if op == 2 else 0
    #write the dictionary to a file
    # with open('test.txt', 'w') as f:
    #     json.dump(my_dict, f)
    #print(my_dict)
    # if aln.query_name == 'da195d4c-6666-4f9e-a336-d2c61650d643':
    #     print(f"\n {get_info(aln)}, consistency: {score/c}")
    return score/c if c != 0 else 0


def reflag_alignments_parallel(input_bam, vcf, output_bam):
    aln_file = pysam.AlignmentFile(input_bam, "rb", threads= 10)
    out_file = pysam.AlignmentFile(output_bam, "wb", template=aln_file, threads= 10)
    alignments = defaultdict(lambda: {'primary': None, 'secondary': []})
    print(f'start')
    switched = []
    time_read = time.time()
    #all_indel_position = process_vcf(vcf)
    consensus_indel_positions = get_all_consensus_indel_positions_parallel(vcf)
        
    # with open('consensus_indel_positions.json', 'r') as fp:
    #     consensus_indel_positions = json.load(fp)
    # consensus_indel_positions = get_all_consensus_indel_positions_parallel(vcf)
    #print(consensus_indel_positions)
    # Writing the dictionary to a file
    print('done getting consensus indel positions')
    with open('consensus_indel_positions.json', 'w') as f:
        json.dump(consensus_indel_positions, f)

    print(f"Time to read: {time.time() - time_read}")

    print(f'start getting alignments')
    for aln in aln_file.fetch():  
        if aln.is_secondary:
            alignments[aln.query_name]['secondary'].append(aln)
        elif aln.is_supplementary and aln.mapping_quality > 10: # filter out low quality supplementary alignments
            #alignments[aln.query_name]['supplementary'].append(aln)
            out_file.write(aln)
        else:
            alignments[aln.query_name]['primary'] = aln 
    print(f'time to get alignments: {time.time() - time_read}')

    print(f'length of alignments: {len(alignments)}')
    remove = 0
    for values in alignments.values():   
        if len(values['secondary']) == 0:
            out_file.write(values['primary'])
        else:        
            #time_primary = time.time()
            primary = values['primary']
            # protect primary alignment from being switched if its mapping quality is greater than 30

            primary_consistency = calculate_consistency(primary, consensus_indel_positions)
            #print(f"\n {get_info(primary)}, primary consistency: {primary_consistency}")
            #print(f"\n {get_info(primary)}, primary consistency: {primary_consistency}")
            highest_consistency = primary_consistency 
            best_alignment = primary
            for s in values['secondary']:
                if s.query_alignment_length > 0.85 * primary.query_alignment_length:
                    secondary_consistency = calculate_consistency(s, consensus_indel_positions)            
                    if secondary_consistency > highest_consistency:
                        highest_consistency = secondary_consistency
                        best_alignment = s

            for s in values['secondary']:
                if s == best_alignment:
                    s.is_secondary = False
                    primary.is_secondary = True
                    s.flag &= ~0x800
                    primary.flag &= ~0x800
                    # print(f"{get_info(s)}, secondary consistency: {highest_consistency}")
                    # print(f"\t Swapped primary and secondary alignments for {primary.query_name}.")
                    switched.append(primary.query_name)
                out_file.write(s)
            out_file.write(primary)
            if max(primary_consistency, highest_consistency) < 0:
                remove += 1
    print(f"Number of reads removed: {remove}")
    with open('switched.read', 'w') as f:
        f.write('\n'.join(switched))

def main_pipeup(bam, pos, ref):
    get_pileup(bam, pos, ref)

def main_get_alignments(bam, seq_id):
    get_alignments(bam, seq_id)

# def main_reflag(bam, output):
#     reflag_alignments_parallel(bam, output)
def main_reflag(bam, vcf, out):
    # Pass the BAM file path instead of the AlignmentFile object
    reflag_alignments_parallel(bam, vcf, out)
def main_filter(bam, out, m):
    filter_secondary(bam, out, m)


if __name__ == "__main__":
    if args.mode == "pileup":
        main_pipeup(args.bam_file, args.pos, args.ref)
    elif args.mode == "alignment":
        main_get_alignments(args.bam_file, args.seq_id)
    elif args.mode == "reflag":
        main_reflag(args.bam_file, args.vcf, args.output)
    elif args.mode == "filter":
        main_filter(args.bam_file, args.output, args.mpileup)
    else:
        exit('No arguments provided. Use -h for help.')