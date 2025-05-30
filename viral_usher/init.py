# 'init' subcommand

import sys
import os
import json
import re
import argparse
import logging
from . import ncbi_helper

def prompt_int_choice(prompt, min_value=1, max_value=None):
    """Get an integer choice from the user, with optional min and max values.
    min_value will be presented as default in case of empty input."""
    while True:
        try:
            choice = input(f"{prompt} [{min_value}]: ")
            choice = int(choice) if choice else min_value
            if min_value <= choice and (max_value is None or choice <= max_value):
                return choice
            else:
                print(f"Please enter a number between {min_value} and {max_value}.")
        except ValueError:
            if (re.match(r"^[Qq](uit)?", choice)):
                print("Exiting...")
                sys.exit(0)
            print("Invalid input. Please enter a number.")

def prompt_taxonomy_id(ncbi, species_name):
    """Prompt the user to select a Taxonomy ID from the list of matches for a given species name."""
    matches = ncbi.get_taxonomy_entries(species_name)
    if matches:
        print(f"Found the following Taxonomy IDs for '{species_name}':")
        for idx, entry in enumerate(matches):
            print(f"{idx+1}. {entry['sci_name']} (Taxonomy ID: {entry['tax_id']})")
        choice = prompt_int_choice("Select the correct species by number", 1, len(matches))
        return matches[choice-1]["tax_id"]
    else:
        print(f"No matches found for '{species_name}'.")
        sys.exit(1)

def prompt_refseq_id(ncbi, taxid):
    """Prompt the user to select a RefSeq ID from the list of available IDs for a given Taxonomy ID."""
    refseq_entries = ncbi.get_refseqs_for_taxid(taxid)
    if refseq_entries:
        print(f"Found the following RefSeq IDs for Taxonomy ID {taxid}:")
        for idx, r in enumerate(refseq_entries):
            label = f"{r['accession']}: {r['title']}"
            if (r["strain"] and r["strain"] != "No strain"):
                label += f" ({r['strain']})"
            print(f"{idx+1}. {label}")
        choice = prompt_int_choice("Select the correct RefSeq ID by number", 1, len(refseq_entries))
        return refseq_entries[choice-1]['accession']
    else:
        print(f"No RefSeq IDs found for Taxonomy ID {taxid}.")
        sys.exit(1)

def write_config(config, config_path):
    """Write a config TOML file to the specified path."""
    with open(config_path, 'w') as f:
        print(f"# usher_viral config for RefSeq {config['refseq_acc']}, taxonomy ID {config['taxonomy_id']}\n", file=f)
        print(f"refseq_acc = '{config['refseq_acc']}'", file=f)
        print(f"taxonomy_id = '{config['taxonomy_id']}'", file=f)
        print(f"refseq_assembly = '{config['refseq_assembly']}'", file=f)
        print(f"workdir = '{config['workdir']}'", file=f)

def handle_init(args):
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    ncbi = ncbi_helper.NcbiHelper()
    if args.refseq:
        refseq_id = args.refseq
        taxid = ncbi.get_taxid_for_refseq(refseq_id)
        if args.taxonomy_id:
            if taxid != args.taxonomy_id:
                print(f"\nWARNING: RefSeq ID {refseq_id} is associated with Taxonomy ID {taxid}, not the provided Taxonomy ID {args.taxonomy_id}.\n")
                taxid = args.taxonomy_id
    elif args.taxonomy_id:
        taxid = args.taxonomy_id
        refseq_id = prompt_refseq_id(ncbi, taxid)
    elif args.species:
        species = args.species
        taxid = prompt_taxonomy_id(ncbi, species)
        refseq_id = prompt_refseq_id(ncbi, taxid)
    else:
        species = input("Enter viral species name: ").strip()
        if not species:
            print("Can't proceed without a species name.")
            sys.exit(1)
        taxid = prompt_taxonomy_id(ncbi, species)
        refseq_id = prompt_refseq_id(ncbi, taxid)
    assembly_id = ncbi.get_assembly_acc_for_refseq_acc(refseq_id)
    if not assembly_id:
        print(f"Could not find assembly ID for RefSeq ID {refseq_id} -- can't download RefSeq.")
        sys.exit(1)

    workdir = input("Enter directory where sequences should be downloaded and trees built [.]: ") or "."
    # Create the workdir if it doesn't exist
    if not os.path.exists(workdir):
        print(f"Creating directory {workdir}...")
        os.makedirs(workdir)

    config_path_default = f"{workdir}/viral_usher_config_{refseq_id}_{taxid}.toml"
    if args.config:
        config_path = args.config
    else:
        config_path = input(f"Enter path for config file [{config_path_default}]: ") or config_path_default
    workdir = os.path.abspath(workdir)
    config = {
        "refseq_acc": refseq_id,
        "refseq_assembly": assembly_id,
        "taxonomy_id": taxid,
        "workdir": workdir,
    }
    print(f"Writing config file to {config_path}")
    write_config(config, config_path)

    print(f"\nReady to roll!  Next, try running 'viral_usher build --config {config_path}'\n")
