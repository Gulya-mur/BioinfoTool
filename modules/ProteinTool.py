from bio_info import GRAVY_AA_VALUES, AMINO_ACIDS_NAMES
from typing import List


VALID_SYMBOLS = set(AMINO_ACIDS_NAMES)


def calc_gravy(seq: str) -> float:
    """
    Calculate GRAVY (grand average of hydropathy) value
    of given amino acids sequence
    """
    gravy_aa_sum = 0
    for amino_ac in seq:
        gravy_aa_sum += GRAVY_AA_VALUES[amino_ac]
    return round(gravy_aa_sum / len(seq), 3)


def calc_total_charge(charged_amino_ac_numbers_list: list,
                      ph_value: float) -> float:
    """
    Calculate the approximate total charge of some amino acid sequence
    for given pH value
    based only on a list of the number of key charged amino acids.
    """
    n_terminal_charge = 1 / (1 + 10 ** (ph_value - 8.2))
    c_terminal_charge = -1 / (1 + 10 ** (3.65 - ph_value))
    cys_charge = -charged_amino_ac_numbers_list[0] / (1 + 10 ** (8.18 - ph_value))
    asp_charge = -charged_amino_ac_numbers_list[1] / (1 + 10 ** (3.9 - ph_value))
    glu_charge = -charged_amino_ac_numbers_list[2] / (1 + 10 ** (4.07 - ph_value))
    tyr_charge = -charged_amino_ac_numbers_list[3] / (1 + 10 ** (10.46 - ph_value))
    his_charge = charged_amino_ac_numbers_list[4] / (1 + 10 ** (ph_value - 6.04))
    lys_charge = charged_amino_ac_numbers_list[5] / (1 + 10 ** (ph_value - 10.54))
    arg_charge = charged_amino_ac_numbers_list[6] / (1 + 10 ** (ph_value - 12.48))
    total_charge = (n_terminal_charge +
                    c_terminal_charge +
                    cys_charge +
                    asp_charge +
                    glu_charge +
                    tyr_charge +
                    his_charge +
                    lys_charge +
                    arg_charge)
    return total_charge


def calc_iso_point(seq: str):
    """
    Calculate approximate isoelectric point of given amino acids sequence
    """
    charged_amino_ac_numbers = []
    for amino_ac in ("C", "D", "E", "Y", "H", "K", "R"):
        charged_amino_ac_numbers.append(seq.count(amino_ac))
    total_charge_tmp = 1
    ph_iso_point = -0.1
    while total_charge_tmp > 0:
        ph_iso_point += 0.1
        total_charge_tmp = calc_total_charge(
            charged_amino_ac_numbers,
            ph_iso_point)
    return round(ph_iso_point, 1)


def transform_to_three_letters(seq: str) -> str:
    """
    Transform 1-letter aminoacid symbols in
    sequence to 3-letter symbols separated by
    hyphens.
    """
    new_name = ''
    for amino_acid in seq:
        new_name += AMINO_ACIDS_NAMES[amino_acid] + '-'
    return new_name[:-1]


def sequence_length(seq: str) -> int:
    """
    Function counts number of aminoacids in
    given sequence
    """
    return len(seq)


def calc_protein_mass(seq: str) -> int:
    """
    Calculate protein molecular weight using the average
    molecular weight of amino acid - 110 Da
    """
    return sequence_length(seq) * 110


def find_heaviest_proteins(proteins: List[str]) -> List[str]:
    """
    Return the sequence of the heaviest protein from list
    """
    protein_mass = {}
    for protein in proteins:
        protein_mass[protein] = calc_protein_mass(protein)
    return max_mass(protein_mass)


def max_mass(protein_mass):
    """
    Count amount of proteins with the same maximum mass and return them
    """
    max_weight = max(protein_mass.values())
    proteins = []
    for protein, weight in protein_mass.items():
        if weight == max_weight:
                proteins.append(protein)
    return proteins


def find_lightest_proteins(proteins: List[str]) -> List[str]:
    """
    Return the sequence of the lightest protein from list
    """
    protein_mass = {}
    for protein in proteins:
        protein_mass[protein] = calc_protein_mass(protein)
    return min_mass(protein_mass)


def min_mass(protein_mass):
    """
    Count amount of proteins with the same minimum mass and return it
    """
    min_weight = min(protein_mass.values())
    proteins = []
    for protein, weight in protein_mass.items():
        if weight == min_weight:
                proteins.append(protein)
    return proteins


def check_sequences(seqs: list):
    """
    Raise ValueError if at least one sequence
    contains non valid symbols
    """
    if not (isinstance(seqs, list)):
        raise ValueError("Enter valid protein sequence")
    for seq in seqs:
        if (not (isinstance(seq, str))) or (not (set(seq.upper()).issubset(VALID_SYMBOLS))):
            raise ValueError("Enter valid protein sequence")


# Didn't place at the beginning because the functions are defined above



