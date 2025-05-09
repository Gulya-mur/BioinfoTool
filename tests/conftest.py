import pytest
from pathlib import Path
from Bio import SeqIO

@pytest.fixture(scope="session")
def example_fastq_path():
    return Path(__file__).parent / "data" / "example.fastq"

@pytest.fixture(scope="session")
def example_path():
    return Path(__file__).parent / "data" 

@pytest.fixture
def filtered_output_dir(tmp_path):

    output_dir = tmp_path.parent / "fastq_filtrator_results"
    output_dir.mkdir(exist_ok=True)
    return output_dir

