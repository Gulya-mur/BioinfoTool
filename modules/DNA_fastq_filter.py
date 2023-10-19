
from typing import Tuple, Union

def cal_gc_content(dna: str) -> int: 
    """Return percent of GC content"""   
    return dna.upper().count('G') + dna.upper().count('C')/ len(dna) * 100

def is_pass_by_gc(dna: str, parameters: Tuple[Union[int, float]]) -> bool:
    """
    Define filter boundatry and return dna if it passed the filtering using cal_gc_content
    """   
    if isinstance(parameters, tuple):
        return parameters[0] <= cal_gc_content(dna) <= parameters[1]
    else:
        return cal_gc_content(dna) <= parameters  

def is_pass_by_length(dna: str, parameters: Tuple[int]) -> bool:
    """
    Define filter boundatry and return dna if it passed the filtering by set range
    """
    if isinstance(parameters, tuple):
        return parameters[0] <= len(dna) <= parameters[1]
    else:
        return len(dna) <= parameters  

def is_pass_by_quality(quality: str, value: int) -> bool:
    """
    Define quality of each nucleotide using ASCII table and Phred-33 scale
    """
    for qual in quality:
        nuc_qual = ord(qual)-33
        return nuc_qual > value

    

