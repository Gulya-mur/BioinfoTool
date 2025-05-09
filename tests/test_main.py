from BioInfo_Tool.main import filter_fastq, DNASequence, RNASequence
import pytest


def test_example_fastq_exists(example_fastq_path):
    assert example_fastq_path.exists()

def test_filter_by_length(example_fastq_path):
    total, passed = filter_fastq(
        input_path=example_fastq_path,
        length_bounds=(5, 10),
        gc_bounds=(0, 100),
        quality_threshold=0,
        output_filename="test_out.fastq"
    )
    assert total == 89
    assert passed == 6

def test_filter_by_gc(example_fastq_path):
    total, passed = filter_fastq(
        input_path=example_fastq_path,
        gc_bounds=(50, 100), 
        quality_threshold=0,
        length_bounds=(0, 100),
        output_filename="test_out.fastq"
    )
    assert total == 89
    assert passed <= 26 

def test_invalid_dna_sequence():
    with pytest.raises(ValueError):
        DNASequence("ATBXZ")

def test_transcription():
    dna = DNASequence("ATGC")
    rna = dna.transcribe()
    assert str(rna) == "AUGC"

def test_reverse_complement():
    dna = DNASequence("ATGC")
    rc = dna.reverse_complement()
    assert str(rc) == "GCAT"

def test_dna_sequence_indexing():
    dna = DNASequence("ATGC")
    assert len(dna) == 4
    assert dna[0] == "A"
    assert str(dna[1:3]) == "TG"

def test_output_file_creation(example_path, filtered_output_dir):
    input_file = example_path / "example_1.fastq"
    input_file.write_text("@read1\nATGC\n+\nIIII\n")

    filter_fastq(str(input_file), output_dir=str(filtered_output_dir))

    output_file = filtered_output_dir / "example_1.fastq"
    assert output_file.exists()