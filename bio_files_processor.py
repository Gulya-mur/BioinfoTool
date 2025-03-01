from typing import Union, TextIO

def convert_multiline_fasta_to_oneline(input_fasta: str, output_file: str) -> TextIO:
    """Function takes fasta file and convert internal data to oneline deleting unnecessary inforamtion"""
    if output_filename == None:
         output_filename = input_fasta
    oneline = []
    with open(input_fasta, 'r') as lines:
        for line in lines:
            if line[0] != '>':
                oneline.append(line.strip())
        full_line = ''.join(oneline)
    with open(output_file, 'w') as outfile:
        outfile.write(full_line)
    

def select_genes_from_gbk_to_fasta(input_gbk, genes: str, output_fasta: Union[str, None], n_before = 1, n_after = 1):
    """Не хватило сил на эту функцию. Хотелось и доп. функции сделать, но больше времени займет. Если возможно позже досдать на проверку
    (даже не за быллы, а для себя), то буду только рада сделать)"""
    pass
