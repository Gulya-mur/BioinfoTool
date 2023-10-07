from DNATool import transcribe, reverse, complement, reverse_complement, process_multiple_sequences

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