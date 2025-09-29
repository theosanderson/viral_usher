import os
import sys
import tomllib


DEFAULT_MIN_LENGTH_PROPORTION = 0.8
DEFAULT_MAX_N_PROPORTION = 0.25
DEFAULT_MAX_PARSIMONY = 1000
DEFAULT_MAX_BRANCH_LENGTH = 10000
DEFAULT_DOCKER_IMAGE = "angiehinrichs/viral_usher"


def read_config(config_path):
    try:
        with open(config_path, 'rb') as f:
            config = tomllib.load(f)
    except FileNotFoundError:
        print(f"Config file '{config_path}' does not exist.", file=sys.stderr)
        raise
    return config


def check_refseq_vs_ref(config, check_paths=True):
    """Ensure that the config specifies either RefSeq info or local reference files, but not both."""
    refseq_acc = config.get('refseq_acc', '')
    refseq_assembly = config.get('refseq_assembly', '')
    ref_fasta = config.get('ref_fasta', '')
    ref_gbff = config.get('ref_gbff', '')
    if (refseq_acc or refseq_assembly) and (ref_fasta or ref_gbff):
        raise ValueError("Config must specify either refseq_acc and refseq_assembly, or ref_fasta and ref_gbff, but not both.")
    if (not refseq_acc or not refseq_assembly) and (not ref_fasta or not ref_gbff):
        raise ValueError("Config must specify either refseq_acc and refseq_assembly, or ref_fasta and ref_gbff.")
    if ref_fasta:
        if check_paths and (not os.path.exists(ref_fasta) or not os.access(ref_fasta, os.R_OK)):
            raise ValueError(f"Reference FASTA file '{ref_fasta}' does not exist or is not readable.")
    if ref_gbff:
        if check_paths and (not os.path.exists(ref_gbff) or not os.access(ref_gbff, os.R_OK)):
            raise ValueError(f"Reference GenBank file '{ref_gbff}' does not exist or is not readable.")


def check_refseq_acc(refseq_acc):
    if not refseq_acc or not isinstance(refseq_acc, str):
        raise ValueError("RefSeq accession must be a non-empty string.")
    if not (refseq_acc.startswith("NC_") or refseq_acc.startswith("NZ_") or refseq_acc.startswith("AC_")):
        raise ValueError("RefSeq accession must start with 'NC_', 'NZ_' or 'AC_'.")


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


def check_config(config, check_paths=True):
    check_taxonomy_id(config['taxonomy_id'])
    check_refseq_vs_ref(config, check_paths)
    if config.get('refseq_acc', ''):
        check_refseq_acc(config['refseq_acc'])
        check_refseq_assembly(config['refseq_assembly'])
    extra_fasta = config.get('extra_fasta', '')
    if extra_fasta and check_paths and not (os.path.exists(extra_fasta) and os.access(extra_fasta, os.R_OK)):
        raise ValueError(f"Fasta file {extra_fasta} does not exist or is not readable.")


def parse_config(config_path):
    """Read the config file and validate its contents."""
    config = read_config(config_path)
    check_config(config)
    return config


def write_config(config, config_path, check_paths=True):
    """Write a config TOML file to the specified path."""
    check_config(config, check_paths)
    with open(config_path, 'w') as f:
        if config.get('refseq_acc'):
            ref_desc = f"RefSeq {config['refseq_acc']}"
        else:
            ref_desc = f"reference genome {config['ref_fasta']}"
        print(f"# viral_usher config for {ref_desc}, taxonomy ID {config['taxonomy_id']}\n", file=f)
        for key in config:
            print(f"{key} = '{config[key]}'", file=f)
