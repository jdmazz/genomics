"""
Efficiently performs various genomics functionality.
@author jdmazz
"""

from operator import itemgetter

START_CODON = {"ATG", "atg"}
STOP_CODONS = {"TGA", "TAG", "TAA", "tga", "tag", "taa"}

BASE_COMPLEMENTS = {"A":"T", "T":"A", "C":"G", "G":"C", "N":"N", "a":"t",
                  "t":"a", "c":"g", "g":"c", "n":"n"}

def read_fasta_file(fasta_file):
    "Reads a multi-FASTA file and returns a sequence dictionary."
    with open(fasta_file) as f:
        seqs = {}
        for line in f:
            line = line.rstrip()
            if line[0] == '>':
                words = line.split()
                name = words[0][1:]
                seqs[name] = ""
            else:
                seqs[name] = seqs[name] + line
        return seqs

def seq_lengths(seqs):
    "Returns a sorted list of tuples: (seq_name, seq_length)."
    seq_lens = [(seq, len(seqs[seq])) for seq in seqs]
    return sorted(seq_lens, key=itemgetter(1))

def reverse_complement(seq):
    "Returns the reverse complement of the seq."
    return "".join([BASE_COMPLEMENTS[base] for base in seq[::-1]])

def start_codon_pos(seq, frame=0):
    "Returns the positions of the start codons in the reading frame."
    assert (frame <= 2 and frame >= 0)
    codon_pos = []
    for i in range(frame, len(seq), 3):
        codon = seq[i:i+3]
        if codon in START_CODON:
            codon_pos.append(i)
    return codon_pos

def get_orfs(seq, frame=0):
    "Get all orfs in specified frame. Returns dictionary {start_codon:orf}."
    assert (frame <= 2 and frame >= 0)
    orfs = {}
    start_codons = start_codon_pos(seq, frame)
    for start_codon in start_codons:
        orf = [] # orf is deleted on each iteration
        for i in range(start_codon, len(seq), 3):
            codon = seq[i:i+3]
            if codon in STOP_CODONS:
                orf.append(codon)
                orfs[start_codon] = "".join(orf) # Add to orfs only if stop codon found
                break
            else:
                orf.append(codon)
    return orfs

def longest_frame_orf(seq, frame=0):
    "Returns longest orf tuple (start_codon, orf_length)"
    orfs = get_orfs(seq, frame)
    if len(orfs) > 0:
        return seq_lengths(orfs)[-1]
    else:
        return (-1, 0)

def longest_orf_in_seq(seq, frames=[0,1,2], use_rc=False):
    "Returns the longest orf tuple in sequence for all reading frames, 5' and 3'."
    rc_seq = reverse_complement(seq)
    max_orf = (-1,0)
    for i in frames: #Fix this
        long_orf = longest_frame_orf(seq, i)
        if long_orf[1] > max_orf[1]:
            max_orf = long_orf
        if use_rc:
            rc_long_orf = longest_frame_orf(rc_seq, i)
            if rc_long_orf[1] > max_orf[1]:
                max_orf = rc_long_orf
    return max_orf

def longest_orfs(seqs, frames=[0,1,2], use_rc=False):
    "Returns sorted list of [name, (start_codon, max_orf_length)]."
    orfs = [(name, longest_orf_in_seq(seqs[name], frames, use_rc)) for name in seqs]
    return sorted(orfs, key=lambda tup: tup[1][1])

def repeats(seqs, n):
    sub_seqs = {}
    for seq in seqs.values():
        for i in range(len(seq)):
            sub_seq = seq[i:i+n]
            if len(sub_seq) != n:
                break
            if sub_seq in sub_seqs:
                sub_seqs[sub_seq] = sub_seqs[sub_seq] + 1
            else:
                sub_seqs[sub_seq] = 1
    sub_tups = [(seq, sub_seqs[seq]) for seq in sub_seqs]
    return sorted(sub_tups, key=itemgetter(1))

def test():
    seqs = read_fasta_file("dna2.fasta")
    
    print("Number of sequences: ", len(seqs), "\n")
    
    print("Sequence lengths: ")
    for name, lengths in seq_lengths(seqs):
        print(name, "\t", lengths)
    print()
        
    print("Longest orf (start_codon_pos, length): ")
    max_orfs = longest_orfs(seqs, [0,1,2])
    for orf in max_orfs:
        print(orf[0], "\t", orf[1])
    print()
        
    print("Repeats: ")
    repeat_sub_seqs = repeats(seqs, 12)
    print(repeat_sub_seqs)
    print()
    
    print("Repeats: ")
    repeat_sub_seqs = repeats(seqs, 7)
    print(repeat_sub_seqs)
