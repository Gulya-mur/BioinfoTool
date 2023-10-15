from DNATool import transcribe, reverse, complement, reverse_complement, process_multiple_sequences
from typing import Dict, List, Tuple, Union
from ProteinTool import calc_gravy, calc_iso_point, transform_to_three_letters, sequence_length, calc_protein_mass, find_heaviest_proteins, find_lightest_proteins, check_sequences
from DNA_fastq_filter import length_bounds, gc_bounds, quality_threshold

def run_dna_rna_tools(operation: str, seqs: List[str]):
    if len(seqs) >= 2:
        return process_multiple_sequences(seqs, operation)
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
    
def read_fastq(seqs: Dict[str, tuple] , gc_bounds: Tuple[Union[int, float]], length_bounds: Tuple[int], quality_threshold: int) -> Dict:

    """
    The main function takes 4 arguments. It checks DNA for compliance with conditions and return new Dict with filtered DNA. 
    """
    filtered_reads = {}
    for seq_name, values in seqs.items():
        seq, seq_quality = values 
        if (is_pass_by_length(seq[0], length_bounds) and is_pass_by_gc(seq[0], gc_bounds) and is_pass_by_quality(seq_quality[1], quality_threshold)):                   
                   filtered_reads[seq_name] = values
    return filtered_reads
