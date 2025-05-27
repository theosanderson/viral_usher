import argparse
from .init import handle_init
from .run import handle_run

def main():
    parser = argparse.ArgumentParser(prog="viral-usher")
    subparsers = parser.add_subparsers(dest="command", required=True)

    # Init subcommand
    init_parser = subparsers.add_parser("init", help="Initialize configuration for a new species/taxon and reference")
    init_parser.add_argument("-r", "--refseq", help="RefSeq accession to use as reference/root, if known")
    init_parser.add_argument("-t", "--taxonomy_id", help="NCBI Taxonomy ID of the viral species, if known")
    init_parser.add_argument("-s", "--species", help="Viral species name (if Taxonomy ID is not known)")
    init_parser.add_argument("--config", type=str, help="Path to config file output")

    # Run subcommand
    run_parser = subparsers.add_parser("run", help="Run the pipeline")
    run_parser.add_argument("--config", type=str, required=True, help="Path to config file input")

    args = parser.parse_args()

    if args.command == "init":
        handle_init(args)
    elif args.command == "run":
        handle_run(args)
    else:
        parser.print_help()
        exit(1)

if __name__ == "__main__":
    main()
