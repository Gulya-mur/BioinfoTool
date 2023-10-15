from bio_info import NUCLEOTIDES, NUCL_COMPLEMENT_MAP
from typing import List

      
def is_dna(seq: str) -> str:    
    for nuc in seq:
        if nuc not in NUCLEOTIDES:
            raise ValueError("It's RNA, not DNA")
    return seq 

def transcribe(seq: str) -> str:
    seq = is_dna(seq)
    return seq.replace("T", "U").replace("t","u")

def reverse(seq: str) -> str:
    return seq[0][::-1]
    
def complement(seq: str) -> str:        
    seq = is_dna(seq)
    return ''.join([NUCL_COMPLEMENT_MAP[nuc] for nuc in seq])          

def reverse_complement(seq: str) -> str:
    return complement(reverse(seq)) 
    
def process_multiple_sequences(seqs: List, operation: str) -> list:
    func = {"transcribe": transcribe,
            "reverse": reverse,
            "complement": complement,
            "reverse_complement": reverse_complement}      
    list_seq = []
    for seq in seqs:
        list_seq.append(func[operation](seq))
    return list_seq

