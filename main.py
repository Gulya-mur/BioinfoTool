from typing import Dict, List, Tuple, Union, TextIO, Optional
import os
from abc import ABC, abstractmethod
from Bio import SeqIO, SeqUtils
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


def save_file_at_dir(
    file_content: Dict[str, Tuple[str, str, str]],
    output_filename: Optional[str],
    input_path: str,
    mode: str = "w",
) -> TextIO:
    if not output_filename:
        output_filename = os.path.basename(input_path)

    output_dir = os.path.join(os.path.dirname(input_path), "fastq_filtrator_results")
    os.makedirs(output_dir, exist_ok=True)
    full_path = os.path.join(output_dir, output_filename)

    records = []
    for header, (seq, comment, qual) in file_content.items():
        record = SeqRecord(
            seq=Seq(seq),
            id=header.strip("@").split()[0],
            description=comment.strip("+").strip(),
            letter_annotations={"phred_quality": [ord(q) - 33 for q in qual]},
        )
        records.append(record)

    with open(full_path, mode) as file:
        SeqIO.write(records, file, "fastq")

    return file


def filter_fastq(
    input_path: str,
    gc_bounds: Union[Tuple[Union[int, float], Union[int, float]]] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
    output_filename: Optional[str] = None,
) -> Dict[str, Tuple[str, str, str]]:

    filtered_seqs = {}

    total_reads = 0
    passed_reads = 0

    for record in SeqIO.parse(input_path, "fastq"):
        total_reads += 1
        gc = SeqUtils.gc_fraction(record.seq) * 100
        if not (gc_bounds[0] <= gc <= gc_bounds[1]):
            continue

        len_seq = len(record.seq)
        if not (length_bounds[0] <= len_seq <= length_bounds[1]):
            continue

        qualities = record.letter_annotations.get("phred_quality", [])
        if not qualities:
            continue
        avg_quality = sum(qualities) / len(qualities)
        if avg_quality < quality_threshold:
            continue

        passed_reads += 1
        filtered_seqs[record.id] = (
            str(record.seq),
            "+",
            "".join([chr(q + 33) for q in qualities]),
        )

    save_file_at_dir(
        file_content=filtered_seqs,
        output_filename=output_filename,
        input_path=input_path,
    )

    return total_reads, passed_reads


class BiologicalSequence(ABC):
    def __init__(self, sequence: str):
        self.sequence = sequence.upper()
        if not self.is_valid():
            raise ValueError("Invalid sequence characters")

    def __len__(self) -> int:
        return len(self.sequence)

    def __getitem__(self, index: int | slice):
        if isinstance(index, slice):
            return self.__class__(self.sequence[index])
        return self.sequence[index]

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}({self.sequence})"

    def __str__(self) -> str:
        return self.sequence

    @abstractmethod
    def is_valid(self) -> bool:
        pass


class NucleicAcidSequence(BiologicalSequence):
    _alphabet = set()
    _nucl_complement_map = {}

    def is_valid(self) -> bool:
        return all(nuc in self._alphabet for nuc in self.sequence)

    def complement(self) -> str:
        return self.__class__(
            "".join(self._nucl_complement_map[nuc] for nuc in self.sequence)
        )

    def reverse(self) -> str:
        return self.__class__(self.sequence[::-1])

    def reverse_complement(self) -> str:
        return self.complement().reverse()


class DNASequence(NucleicAcidSequence):
    _alphabet = {"A", "T", "G", "C"}
    _nucl_complement_map = {"A": "T", "T": "A", "G": "C", "C": "G"}

    def transcribe(self):
        return RNASequence(self.sequence.replace("T", "U"))


class RNASequence(NucleicAcidSequence):
    _alphabet = {"A", "U", "G", "C"}
    _nucl_complement_map = {"A": "U", "U": "A", "C": "G", "G": "C"}


class AminoAcidSequence(BiologicalSequence):
    _alphabet = {
        "A",
        "R",
        "N",
        "D",
        "V",
        "H",
        "G",
        "Q",
        "E",
        "I",
        "L",
        "K",
        "M",
        "P",
        "S",
        "Y",
        "T",
        "W",
        "F",
        "C",
    }
    amino_acid_names = {
        "A": "Ala",
        "R": "Arg",
        "N": "Asn",
        "D": "Asp",
        "V": "Val",
        "H": "His",
        "G": "Gly",
        "Q": "Gln",
        "E": "Glu",
        "I": "Ile",
        "L": "Leu",
        "K": "Lys",
        "M": "Met",
        "P": "Pro",
        "S": "Ser",
        "Y": "Tyr",
        "T": "Thr",
        "W": "Trp",
        "F": "Phe",
        "C": "Cys",
    }

    def is_valid(self) -> bool:
        return all(aa in self._alphabet for aa in self.sequence)

    def transform_to_three_letters(self) -> str:
        return "".join(self.amino_acid_names[aa] for aa in self.sequence)


if __name__ == "__main__":
    # Пример работы с ДНК
    dna = DNASequence("ATGC")
    print(f"Исходная ДНК: {dna}")
    print(f"Обратная комплементарная ДНК: {dna.reverse_complement()}")
    print(f"Транскрибированная РНК: {dna.transcribe()}")

    # Пример фильтрации FASTQ
    filtered_1 = filter_fastq(
        input_path="example.fastq",
        gc_bounds=(20, 80),
        length_bounds=(50, 150),
        quality_threshold=20,
        output_filename="filtered.fastq",
    )
    print(f"После фильтрации c условиями: {filtered_1}")
    filtered_2 = filter_fastq(input_path="example.fastq")
    print(f"После фильтрации со значениями по умолчанию: {filtered_2}")
