from typing import Tuple, Union, Optional
import os
from abc import ABC, abstractmethod
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from Bio.SeqRecord import SeqRecord
import argparse
from loguru import logger

logger.add("fastq_filtrator.log", format="{time} - {level} - {message}", level="INFO")

def write(
    file_content: list[SeqRecord],
    input_path: str,
    output_filename: Optional[str],
    output_dir: Optional[str] = None
) -> None:
    
    try:
        if not output_filename:
            output_filename = os.path.basename(input_path)


        if not output_dir:
            input_dir = os.path.dirname(os.path.abspath(input_path))
            output_dir = os.path.join(input_dir, "../fastq_filtrator_results")

        os.makedirs(output_dir, exist_ok=True)
        full_path = os.path.join(output_dir, output_filename)

        with open(full_path,'w') as file:
            SeqIO.write(file_content, file, "fastq")
        
        logger.info(f"Сохранено {len(file_content)} ридов в {full_path}")

    except Exception as e:
        logger.error(f"Ошибка сохранения: {e}")
        raise


def filter_fastq(
    input_path: str,
    gc_bounds: Union[Tuple[Union[int, float], Union[int, float]]] = (0, 100),
    length_bounds: Union[Tuple[int, int], int] = (0, 2**32),
    quality_threshold: int = 0,
    output_filename: Optional[str] = None,
    output_dir: Optional[str] = None
) -> Tuple[int, int]:

    filtered_seqs = []

    total_reads = 0
    passed_reads = 0

    try:
        for record in SeqIO.parse(input_path, "fastq"):
            total_reads += 1
            gc = gc_fraction(record.seq) * 100
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
            filtered_seqs.append(record)

        if filtered_seqs:
            write(filtered_seqs, input_path, output_filename, output_dir)
    except Exception as e:
        logger.error(f"Ошибка при обработке файла {input_path}: {e}")
        raise

    # write(
    #     file_content=filtered_seqs,
    #     output_filename=output_filename,
    #     input_path=input_path,
    #     output_dir=output_dir
    # )
    logger.info(f"Прочитано ридов: {total_reads}, прошло фильтрацию: {passed_reads}")
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

def main():
    parser = argparse.ArgumentParser(description="Фильтрация FASTQ файлов")
    parser.add_argument("input_path", type=str, help="Путь к исходному FASTQ файлу")
    parser.add_argument('--gc_bounds', type=float, nargs=2, default=(0, 100), help="Диапазон GC содержимого (по умолчанию: 0-100)")
    parser.add_argument('--length_bounds', type=float, nargs=2, default=(0, 2**32), help="Диапазон длины последовательностей (по умолчанию: 0-2^32)")
    parser.add_argument('--quality_threshold', type=float, default=0, help="Минимальное среднее качество для сохранения (по умолчанию: 0)")
    parser.add_argument('--output_filename', type=str, nargs=1, default=None, help="Имя выходного файла (по умолчанию будет использовано имя исходного файла)")

    args = parser.parse_args()

    filtered = filter_fastq(
        input_path=args.input_path,
        gc_bounds=(args.gc_bounds[0], args.gc_bounds[1]),
        length_bounds=(args.length_bounds[0], args.length_bounds[1]),
        quality_threshold=args.quality_threshold,
        output_filename=args.output_filename,
    )

    print(f"После фильтрации: {filtered}")


if __name__ == "__main__":
    main()

