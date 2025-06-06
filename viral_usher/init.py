# 'init' subcommand: get user's preferences and write config file

import sys
import os
import importlib.metadata
import json
import re
import logging
from . import ncbi_helper
from .config import write_config

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

def check_optional_file_readable(filename):
    """If filename is empty, that's fine; if non-empty, make sure it's a readable file"""
    if filename:
        if not os.path.exists(filename):
            return False, f"Sorry, file '{filename}' does not exist"
        if not os.path.isfile(filename):
            return False, f"Sorry, file '{filename}' does not appear to be a regular file"
        else:
            try:
                with open(filename, "r") as f:
                    pass
            except IOError:
                return False, f"Sorry, unable to read file '{filename}'."
    return True, None

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

    if args.extra_fasta:
        extra_fasta = args.extra_fasta
        ok, error_message = check_optional_file_readable(extra_fasta)
        if not ok:
            print(f"{error_message}\nPlease try again with a different file for --extra-fasta.", file=sys.stderr)
            sys.exit(1)
    else:
        ok = False
        while not ok:
            extra_fasta = input("If you have your own fasta file, then enter its path []:")
            ok, error_message = check_optional_file_readable(extra_fasta)
            if not ok:
                print(error_message)

    workdir = input("Enter directory where sequences should be downloaded and trees built [.]: ") or "."
    # Create the workdir if it doesn't exist
    if not os.path.exists(workdir):
        print(f"Creating directory {workdir}...")
        try:
            os.makedirs(workdir)
        except:
            print(f"Sorry, can't create directory path '{workdir}', please enter a different one.")
            sys.exit(1)

    config_path_default = f"{workdir}/viral_usher_config_{refseq_id}_{taxid}.toml"
    if args.config:
        config_path = args.config
    else:
        config_path = input(f"Enter path for config file [{config_path_default}]: ") or config_path_default
    workdir = os.path.abspath(workdir)
    viral_usher_version = importlib.metadata.version('viral_usher')
    config = {
        "viral_usher_version": viral_usher_version,
        "refseq_acc": refseq_id,
        "refseq_assembly": assembly_id,
        "taxonomy_id": taxid,
        "extra_fasta": extra_fasta,
        "workdir": workdir,
    }
    print(f"Writing config file to {config_path}")
    try:
        write_config(config, config_path)
    except:
        print(f"Sorry, can't write to file '{config_path}'; please try again with a different path.")
        sys.exit(1)

    print(f"\nReady to roll!  Next, try running 'viral_usher build --config {config_path}'\n")
