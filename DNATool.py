from bio_info import NUCLEOTIDES, NUCL_COMPLEMENT_MAP
from main_functions import run_dna_rna_tools
      
def is_dna(seq):    
    for nuc in seq:
        if nuc not in NUCLEOTIDES:
            raise ValueError("It's RNA, not DNA")
    return seq 

def transcribe(seq: str):
    seq = is_dna(seq)
    return seq.replace("T", "U").replace("t","u")

def reverse(seq: str):
    return seq[0][::-1]
    
def complement(seq: str):        
    seq = is_dna(seq)
    return ''.join([NUCL_COMPLEMENT_MAP[nuc] for nuc in seq])          

def reverse_complement(seq: str):
    return complement(reverse(seq)) 
    
def process_multiple_sequences(seqs, operation):       
    list_seq = []
    for seq in seqs:
        list_seq.append(run_dna_rna_tools(operation, [seq]))
    return list_seq

    