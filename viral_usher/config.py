import sys
import tomllib


def read_config(config_path):
    with open(config_path, 'rb') as f:
        config = tomllib.load(f)
    return config


def check_refseq_acc(refseq_acc):
    if not refseq_acc or not isinstance(refseq_acc, str):
        raise ValueError("RefSeq accession must be a non-empty string.")
    if not refseq_acc.startswith("NC_") and not refseq_acc.startswith("NZ_"):
        raise ValueError("RefSeq accession must start with 'NC_' or 'NZ_'.")


def check_taxonomy_id(taxonomy_id):
    if not taxonomy_id or not isinstance(taxonomy_id, str):
        raise ValueError("Taxonomy ID must be a non-empty string.")
    if not taxonomy_id.isdigit():
        raise ValueError("Taxonomy ID must be a numeric string.")


def check_refseq_assembly(refseq_assembly):
    if not refseq_assembly or not isinstance(refseq_assembly, str):
        raise ValueError("RefSeq assembly must be a non-empty string.")
    if not refseq_assembly.startswith("GCA_") and not refseq_assembly.startswith("GCF_"):
        raise ValueError("RefSeq assembly must start with 'GCA_' or 'GCF_'.")


def parse_config(config_path):
    """Read the config file and validate its contents."""
    config = read_config(config_path)
    print(f"refseq_acc: {config['refseq_acc']}")
    print(f"taxonomy_id: {config['taxonomy_id']}")
    print(f"refseq_assembly: {config['refseq_assembly']}")
    try:
        check_refseq_acc(config['refseq_acc'])
        check_taxonomy_id(config['taxonomy_id'])
        check_refseq_assembly(config['refseq_assembly'])
    except (ValueError, FileNotFoundError) as e:
        print(f"Error in config file {config_path}: {e}", file=sys.stderr)
        raise
    return config


def write_config(config, config_path):
    """Write a config TOML file to the specified path."""
    with open(config_path, 'w') as f:
        print(f"# viral_usher config for RefSeq {config['refseq_acc']}, taxonomy ID {config['taxonomy_id']}\n", file=f)
        for key in config:
            print(f"{key} = '{config[key]}'", file=f)
