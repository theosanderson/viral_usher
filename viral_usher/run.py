# 'run' subcommand
import gzip
import io
import json
import lzma
import tomllib
import os
import shutil
import subprocess
import sys
import zipfile
from Bio import SeqIO

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
    
def check_zip_file(zip_file):
    if not zip_file or not isinstance(zip_file, str):
        raise ValueError("Zip file path must be a non-empty string.")
    if not zip_file.endswith(".zip"):
        raise ValueError("Zip file must have a '.zip' extension.")
    if not os.path.exists(zip_file):
        raise FileNotFoundError(f"Zip file '{zip_file}' does not exist.")

def parse_config(config_path):
    """Read the config file and validate its contents."""
    config = read_config(config_path)
    print(f"workdir: {config['workdir']}")
    print(f"refseq_acc: {config['refseq_acc']}")
    print(f"taxonomy_id: {config['taxonomy_id']}")
    print(f"refseq_assembly: {config['refseq_assembly']}")
    print(f"refseq_zip: {config['refseq_zip']}")
    print(f"genbank_zip: {config['genbank_zip']}")
    try:
        check_refseq_acc(config['refseq_acc'])
        check_taxonomy_id(config['taxonomy_id'])
        check_refseq_assembly(config['refseq_assembly'])
        check_zip_file(config['refseq_zip'])
        check_zip_file(config['genbank_zip'])
    except (ValueError, FileNotFoundError) as e:
        print(f"Configuration error: {e}")
        raise
    return config

def unpack_refseq_zip(refseq_zip, workdir, refseq_acc):
    """Extract and rename the fasta and gbff files from the RefSeq zip."""
    with zipfile.ZipFile(refseq_zip, 'r') as zip_ref:
        fasta_path = os.path.join(workdir, 'refseq.fasta')
        gbff_path = os.path.join(workdir, 'refseq.gbff')
        length = 0
        fasta_found = False
        gbff_found = False
        for name in zip_ref.namelist():
            if name.endswith('.fna'):
                with io.TextIOWrapper(zip_ref.open(name, 'r'), encoding='utf-8') as fasta_in, open(fasta_path, 'w') as fasta_out:
                    for record in SeqIO.parse(fasta_in, 'fasta'):
                        if record.id == refseq_acc:
                            record.description = record.id
                            SeqIO.write(record, fasta_out, 'fasta')
                            length = len(record.seq)
                            fasta_found = True
            elif name.endswith('.gbff'):
                with zip_ref.open(name, 'r') as gbff_in, open(gbff_path, 'wb') as gbff_out:
                    shutil.copyfileobj(gbff_in, gbff_out)
                gbff_found = True
        if not fasta_found:
            raise ValueError(f"Failed to find .fna file with {refseq_acc} in {refseq_zip}.")
        if not gbff_found:
            raise ValueError(f"Failed to find .gbff file in {refseq_zip}.")
    return fasta_path, gbff_path, length

def unpack_genbank_zip(genbank_zip, workdir, min_length, max_N_proportion):
    """Process the fasta and data_report.jsonl files from the GenBank zip,
    filtering the fasta by length and proportion of ambiguous characters and
    keeping only the accession part of sequence names, and extracting basic
    metadata from data_report.jsonl into much more compact TSV.
    Compress fasta with lzma (xz) and TSV with gzip"""
    import json
    fasta_out_path = os.path.join(workdir, 'genbank.fasta.xz')
    tsv_out_path = os.path.join(workdir, 'data_report.tsv.gz')
    fasta_written = False
    report_written = False
    with zipfile.ZipFile(genbank_zip, 'r') as zip_ref:
        for name in zip_ref.namelist():
            if name.endswith('.fna'):
                # Filter fasta file and chop descriptions to get accession as name in nextclade
                with io.TextIOWrapper(zip_ref.open(name, 'r'), encoding='utf-8') as fasta_in, lzma.open(fasta_out_path, 'wt') as fasta_out:
                    for record in SeqIO.parse(fasta_in, 'fasta'):
                        # Filter sequences by length and proportion of ambiguous characters
                        if len(record.seq) < min_length or record.seq.count('N') / len(record.seq) > max_N_proportion:
                            continue
                        # Chop sequence name after first space to limit to accession only
                        record.description = record.id
                        SeqIO.write(record, fasta_out, 'fasta')
                print(f"Wrote fasta file: {fasta_out_path}")
                fasta_written = True
            elif name.endswith('data_report.jsonl'):
                # Extract only the bits we want from data_report.jsonl to data_report.tsv
                with zip_ref.open(name) as jsonl_in, gzip.open(tsv_out_path, 'wt') as tsv_out:
                    # Write header
                    tsv_out.write("\t".join(["accession", "isolate", "length", "date", "biosample", "submitter"]) + "\n")
                    for line in jsonl_in:
                        item = json.loads(line)
                        # Only keep selected fields; adjust as needed
                        accession = item.get('accession', '')
                        isolate = item.get('isolate', {}).get('name', '')
                        length = item.get('length', '')
                        date = item.get(isolate, {}).get('collectionDate', '')
                        location = item.get('location', {}).get('geographicLocation', '')
                        biosample = item.get('biosample', '')
                        submitter = item.get('submitter', {}).get('affiliation', '')
                        submitter_country = item.get('submitter', {}).get('country', '')
                        if submitter and submitter_country:
                            submitter = f"{submitter}, {submitter_country}"
                        tsv_out.write("\t".join([accession, isolate, str(length), date, biosample, submitter]) + "\n")

                print(f"Wrote filtered TSV file: {tsv_out_path}")
                report_written = True
    if not fasta_written:
        raise ValueError(f"Failed to find .fna file in {genbank_zip}.")
    if not report_written:
        raise ValueError(f"Failed to find data_report.jsonl file in {genbank_zip}.")
    return fasta_out_path, tsv_out_path

def align_sequences(refseq_fasta, genbank_fasta, workdir):
    """Run nextclade to align the filtered sequences to the reference."""
    msa_fasta = os.path.join(workdir, 'msa.fasta.xz')
    command = ['nextclade', 'run', '--input-ref', refseq_fasta, '--include-reference', '--output-fasta', msa_fasta, genbank_fasta]
    try:
        result = subprocess.run(command, check=True)
        print(f"Nextclade alignment completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Nextclade alignment failed: {e.stderr}")
        raise
    return msa_fasta

def msa_to_vcf(msa_fasta, workdir):
    """Run faToVcf as a pipeline: lzma-decompress msa_fasta.xz | faToVcf | gzip > vcf_output.gz, using Python's lzma/gzip libs and chunked I/O."""
    import subprocess
    vcf_output = os.path.join(workdir, 'msa.vcf.gz')
    with lzma.open(msa_fasta, 'rb') as fasta_in, gzip.open(vcf_output, 'wb') as vcf_out:
        fatovcf_proc = subprocess.Popen(
            ['faToVcf', '-includeNoAltN', 'stdin', 'stdout'],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE
        )
        # Read and write in chunks for efficiency
        for chunk in iter(lambda: fasta_in.read(8192), b''):
            fatovcf_proc.stdin.write(chunk)
        fatovcf_proc.stdin.close()
        for chunk in iter(lambda: fatovcf_proc.stdout.read(8192), b''):
            vcf_out.write(chunk)
        fatovcf_proc.stdout.close()
        fatovcf_proc.wait()
    print(f"Wrote VCF file: {vcf_output}")
    return vcf_output

def make_empty_tree(workdir):
    """Create an empty tree file to be used as a starting point for usher-sampled."""
    empty_tree_path = os.path.join(workdir, 'empty_tree.nwk')
    with open(empty_tree_path, 'w') as f:
        f.write("()")
    print(f"Created empty tree file: {empty_tree_path}")
    return empty_tree_path

def run_usher_sampled(tree, vcf, workdir):
    pb_out = os.path.join(workdir, 'usher_sampled.pb.gz')
    command = ['usher-sampled', '-t', tree, '-v', vcf, '-o', pb_out]
    try:
        result = subprocess.run(command, check=True)
        print(f"Usher-sampled completed successfully")
    except subprocess.CalledProcessError as e:
        print(f"Usher-sampled failed: {e.stderr}")
        raise
    return pb_out

def handle_run(args):
    print(f"Running pipeline with config file: {args.config}")
    try:
        config = parse_config(args.config)
        refseq_fasta, refseq_gbff, refseq_length = unpack_refseq_zip(config['refseq_zip'], config['workdir'], config['refseq_acc'])
        # TODO: Parameterize these into config
        min_length_proportion = 0.8
        min_length = int(refseq_length * min_length_proportion)
        max_N_proportion = 0.25
        print(f"Using min_length_proportion={min_length_proportion} ({min_length} bases) and max_N_proportion={max_N_proportion}")
        genbank_fasta, data_report = unpack_genbank_zip(config['genbank_zip'], config['workdir'], min_length, max_N_proportion)
        msa_fasta = align_sequences(refseq_fasta, genbank_fasta, config['workdir'])
        vcf = msa_to_vcf(msa_fasta, config['workdir'])
        empty_tree = make_empty_tree(config['workdir'])
        preopt = run_usher_sampled(empty_tree, vcf, config['workdir'])
    # TODO: Run usher-sampled to build a tree
    # TODO: Run matOptimize to clean up after usher-sampled
    # TODO: Rename tree names to include location, isolate name and date; prepare metadata for Taxonium
    # TODO: Run usher_to_taxonium to visualize the tree
    except (ValueError, FileNotFoundError, FileNotFoundError, PermissionError):
        sys.exit(1)
