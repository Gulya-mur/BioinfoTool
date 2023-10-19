from typing import Dict, List, Tuple, Union
import os.path
from modules.ProteinTool import calc_gravy, calc_iso_point, transform_to_three_letters, sequence_length, calc_protein_mass, find_heaviest_proteins, find_lightest_proteins, check_sequences
from modules.DNA_fastq_filter import is_pass_by_gc, is_pass_by_length, is_pass_by_quality
from modules.DNATool import transcribe, reverse, complement, reverse_complement, process_multiple_sequences

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
    
def read_fastq(fastq_file: str) -> Dict:
    with open(fastq_file, 'r') as fastq:
        fastq_dict = {}
        while True:
            name = fastq.readline()
            seq = fastq.readline().rstrip()
            comment = fastq.readline()
            qual = fastq.readline().rstrip()
            if len(seq) == 0:
                break
            fastq_dict[name] = (seq,comment, qual)
    fastq.closed
    return fastq_dict


def filter_fastq(input_path: str , gc_bounds: Tuple[Union[int, float]], length_bounds: Tuple[int], quality_threshold: int, output_filename: Union[str, None]) -> Dict:

    """
    The main function takes 4 arguments. It checks DNA for compliance with conditions and return new Dict with filtered DNA. 
    """
    def final_result(filtered_reads: Dict, output_filename: Union[str, None]):
        save_path = 'fastq_filtrator_resuls/'
        if output_filename is None:
            output_filename = input_path 
        with open (output_filename, 'w') as outfile:
            for name, values in filtered_reads.items():
                    seq, comment, quality = values
                    outfile.write(name)
                    outfile.write(seq)
                    outfile.write('\n')
                    outfile.write(comment)
                    outfile.write(quality)
                    outfile.write('\n')
#        os.path.join(save_path, output_filename)        
                
    fastq_dict = read_fastq(input_path)
    filtered_reads = {}
    for seq_name, values in fastq_dict.items():
        seq, _, seq_quality = values 
        if (is_pass_by_length(seq, length_bounds) and is_pass_by_gc(seq, gc_bounds) and is_pass_by_quality(seq_quality, quality_threshold)):                   
                   filtered_reads[seq_name] = values
    print(filtered_reads)
#    return final_result(filtered_reads, output_filename)

filter_fastq('example_fastq.fastq', (20, 80), (0, 2**32), 0, 'output_filename.fastq')





            

         
     