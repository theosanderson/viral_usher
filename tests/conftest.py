import pytest
from pathlib import Path


@pytest.fixture(scope="session")
def config_zika_no_workdir_path():
    """Return path to test config file."""
    return Path(__file__).parent / "fixtures" / "zika_no_workdir.toml"


@pytest.fixture(scope="session")
def config_zika_no_workdir_data():
    """Return expected config data for testing."""
    return {
        "viral_usher_version": "0.2.0",
        "refseq_acc": "NC_035889.1",
        "refseq_assembly": "GCF_002366285.1",
        "taxonomy_id": "64320",
        "extra_fasta": ""
    }
