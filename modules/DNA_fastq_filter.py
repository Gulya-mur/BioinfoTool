
from typing import Tuple, Union

def cal_gc_content(dna: str) -> int: 
    """Return percent of GC content"""   
    return dna.upper().count('G') + dna.upper().count('C')/ len(dna) * 100

def is_pass_by_gc(dna: str, parameters: Tuple[Union[int, float]]) -> str:
    """
    Define filter boundatry and return dna if it passed the filtering using cal_gc_content
    """   
    if isinstance(parameters, tuple):
        if parameters[0] <= cal_gc_content(dna) <= parameters[1]:
            return dna
    else:
        if cal_gc_content(dna) <= parameters:
            return dna   

def is_pass_by_length(dna: str, parameters: Tuple[int]) -> str:
    """
    Define filter boundatry and return dna if it passed the filtering by set range
    """
    if isinstance(parameters, tuple):
        if parameters[0] <= len(dna) <= parameters[1]:
            return dna
    else:
        if len(dna) <= parameters:
            return dna
   

def is_pass_by_quality(quality: str, value: int) -> bool:
    """
    Define quality of each nucleotide using ASCII table and Phred-33 scale
    """
    for qual in quality:
        nuc_qual = ord(qual)-33
        if nuc_qual < value:
            break
    

