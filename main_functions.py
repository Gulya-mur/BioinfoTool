from DNATool import *
from typing import List
from ProteinTool import *

def run_dna_rna_tools(operation: str, seqs: List[str]):
    if len(seqs) >= 2 and operation == "process_multiple_sequences":
        return process_multiple_sequences(seqs)
    else:
        if operation == "transcribe":
            return transcribe(seqs)
        if operation == "reverse":
            return reverse(seqs)
        if operation == "complement":
            return complement(seqs)
        if operation == "reverse_complement":
            return reverse_complement(seqs)
        
FUNC_STR_INPUT = {
    'gravy': calc_gravy,
    'iso': calc_iso_point,
    'rename': transform_to_three_letters,
    'lengths': sequence_length,
    'molw': calc_protein_mass}

FUNC_LIST_INPUT = {
    'heavy': find_heaviest_proteins,
    'light': find_lightest_proteins}        

        
def process_seqs(option: str, seqs: List[str]):
    """
    Perform some simple operations on amino acids sequences.
    """
    check_sequences(seqs)
    if option in FUNC_STR_INPUT.keys():
        results = []
        for seq in seqs:
            result_tmp = FUNC_STR_INPUT[option](seq.upper())
            results.append(result_tmp)
        return results
    elif option in FUNC_LIST_INPUT.keys():
        return FUNC_LIST_INPUT[option](seqs)
    else:
        raise ValueError("Enter valid operation")