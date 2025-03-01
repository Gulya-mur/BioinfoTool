
from typing import Tuple, Union

def cal_gc_content(dna: str) -> int: 
    """Return percent of GC content"""   
    return dna.count('G') + dna.count('C')/ len(dna) * 100

def gc_bounds(dna: str, value: Tuple[Union[int, float]]) -> str:
    """
    Define filter boundatry and return dna if it passed the filtering using cal_gc_content
    """   
    if isinstance(value, tuple):
        if value[0] <= cal_gc_content(dna) <= value[1]:
            return dna
    else:
        if cal_gc_content(dna) <= value:
            return dna   

def length_bounds(dna: str, value: Tuple[int]) -> str:
    """
    Define filter boundatry and return dna if it passed the filtering by set range
    """
    if isinstance(value, tuple):
        if value[0] <= len(dna) <= value[1]:
            return dna
    else:
        if len(dna) <= value:
            return dna
   

def quality_threshold(quality: str, value: int):
    """
    Define quality of each nucleotide using ASCII table and Phred-33 scale
    """
    for qual in quality:
        nuc_qual = ord(qual)-33
        if nuc_qual < value:
            break
    return True

