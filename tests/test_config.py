import pytest
import tempfile
from pathlib import Path
import viral_usher.config


@pytest.fixture(scope="session")
def written_config_file(config_zika_no_workdir_path, config_zika_no_workdir_data):
    """Create a config file once for the entire test session."""
    # Create temporary directory that persists for the session
    with tempfile.TemporaryDirectory() as temp_dir:
        # Temporary path for written config file
        config_path = Path(temp_dir) / "test_config.toml"

        # Test data to write: add missing workir (in temp_dir) to zika example
        config_data = config_zika_no_workdir_data
        config_data["workdir"] = Path(temp_dir) / "zika"

        # Write the config file
        viral_usher.config.write_config(config_data, config_path)

        # Verify it was written correctly
        assert config_path.exists()

        yield config_path
        # File gets cleaned up automatically when temp_dir exits


def test_check_refseq_acc_good():
    assert viral_usher.config.check_refseq_acc("NC_012345.1") is None


def test_check_refseq_acc_bad():
    with pytest.raises(ValueError):
        viral_usher.config.check_refseq_acc(None)
    with pytest.raises(ValueError):
        viral_usher.config.check_refseq_acc("")
    with pytest.raises(ValueError):
        viral_usher.config.check_refseq_acc("AB034985")


def test_check_taxonomy_id_good():
    assert viral_usher.config.check_taxonomy_id("12345") is None


def test_check_taxonomy_id_bad():
    with pytest.raises(ValueError):
        viral_usher.config.check_taxonomy_id(None)
    with pytest.raises(ValueError):
        viral_usher.config.check_taxonomy_id("")
    with pytest.raises(ValueError):
        viral_usher.config.check_taxonomy_id("measles")


def test_check_refseq_assembly_good():
    assert viral_usher.config.check_refseq_assembly("GCA_123456789.1") is None
    assert viral_usher.config.check_refseq_assembly("GCF_123456789.1") is None


def test_check_refseq_assembly_bad():
    with pytest.raises(ValueError):
        viral_usher.config.check_refseq_assembly(None)
    with pytest.raises(ValueError):
        viral_usher.config.check_refseq_assembly("")
    with pytest.raises(ValueError):
        viral_usher.config.check_refseq_assembly("NC_012345.1")


def test_read_config_no_workdir(config_zika_no_workdir_path, config_zika_no_workdir_data):
    result = viral_usher.config.read_config(config_zika_no_workdir_path)
    assert result == config_zika_no_workdir_data


def test_read_config_file_not_found():
    with pytest.raises(FileNotFoundError):
        viral_usher.config.read_config("nonexistent.toml")


def test_write_config_creates_file(written_config_file):
    """Test that the config file was created and exists."""
    assert written_config_file.exists()
    assert written_config_file.suffix == ".toml"


def test_read_written_config(written_config_file):
    """Test reading the config file that was written."""
    config_data = viral_usher.config.read_config(written_config_file)
    assert config_data["refseq_acc"] == "NC_035889.1"
    assert config_data["taxonomy_id"] == "64320"
