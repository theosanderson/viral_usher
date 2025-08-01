import argparse
import gzip
import io
import json
import lzma
import re
import shutil
import subprocess
import sys
import time
import zipfile
from Bio import SeqIO
from . import config
from . import ncbi_helper


def run_command(command, stdout_filename=None, stderr_filename=None, fail_ok=False):
    """Run a command, complain and raise if it fails (unless fail_ok, then just return False)"""
    success = False
    stdout = None
    stderr = None
    if stdout_filename:
        stdout = open(stdout_filename, 'w')
    if stderr_filename:
        stderr = open(stderr_filename, 'w')
    try:
        subprocess.run(command, check=True, stdout=stdout, stderr=stderr)
        success = True
    except subprocess.CalledProcessError as e:
        if fail_ok:
            pass
        else:
            if e.stderr:
                print(f"{command[0]} failed: {e.stderr}\n{command}", file=sys.stderr)
            elif stderr_filename:
                print(f"{command[0]} failed, see {stderr_filename}\n{command}", file=sys.stderr)
            else:
                print(f"{command[0]} failed\n{command}", file=sys.stderr)
            raise
    if stdout:
        stdout.close()
    if stderr:
        stderr.close()
    return success


def unpack_refseq_zip(refseq_zip, refseq_acc):
    """Extract and rename the fasta and gbff files from the RefSeq zip."""
    with zipfile.ZipFile(refseq_zip, 'r') as zip_ref:
        fasta_path = 'refseq.fasta'
        gbff_path = 'refseq.gbff'
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


def passes_seq_filter(record, min_length, max_N_proportion):
    """Fasta record passes the filter if it is at least min_length bases and has at most max_N_proportion of Ns"""
    return len(record.seq) >= min_length and record.seq.count('N') / len(record.seq) <= max_N_proportion


def filter_genbank_sequences(fasta_in, fasta_out, min_length, max_N_proportion, fasta_out_path):
    # Filter fasta file and chop descriptions to get accession as name in nextclade
    start_time = time.perf_counter()
    record_count = 0
    passed_count = 0
    for record in SeqIO.parse(fasta_in, 'fasta'):
        record_count += 1
        # Filter sequences by length and proportion of ambiguous characters
        if not passes_seq_filter(record, min_length, max_N_proportion):
            continue
        passed_count += 1
        # Chop sequence name after first space to limit to accession only
        record.description = record.id
        SeqIO.write(record, fasta_out, 'fasta')
    elapsed_time = time.perf_counter() - start_time
    print(f"Processed {record_count} sequences from GenBank and wrote fasta file with {passed_count} sequences passing filters: {fasta_out_path} in {elapsed_time:.1f}s")


def extract_metadata(jsonl_in, tsv_out, tsv_out_path):
    # Extract only the bits we want from data_report.jsonl to data_report.tsv
    start_time = time.perf_counter()
    # Write header
    tsv_out.write("\t".join(["accession", "isolate", "country", "location", "date", "length", "biosample", "submitter"]) + "\n")
    for line in jsonl_in:
        item = json.loads(line)
        # Only keep selected fields; adjust as needed
        accession = item.get('accession', '')
        isolate = item.get('isolate', {}).get('name', '')
        length = item.get('length', '')
        date = item.get('isolate', {}).get('collectionDate', '')
        location = item.get('location', {}).get('geographicLocation', '')
        country = location.split(":")[0]
        location = ":".join(location.split(":")[1:]).strip()
        biosample = item.get('biosample', '')
        submitter = item.get('submitter', {}).get('affiliation', '')
        submitter_country = item.get('submitter', {}).get('country', '')
        if submitter and submitter_country:
            submitter = f"{submitter}, {submitter_country}"
        tsv_out.write("\t".join([accession, isolate, country, location, date, str(length), biosample, submitter]) + "\n")
    elapsed_time = time.perf_counter() - start_time
    print(f"Wrote metadata TSV file: {tsv_out_path} in {elapsed_time:.1f}s")


def unpack_genbank_zip(genbank_zip, min_length, max_N_proportion):
    """Process the fasta and data_report.jsonl files from the GenBank zip,
    filtering the fasta by length and proportion of ambiguous characters and
    keeping only the accession part of sequence names, and extracting basic
    metadata from data_report.jsonl into much more compact TSV.
    Compress fasta with lzma (xz) and TSV with gzip"""
    fasta_out_path = 'genbank.fasta.xz'
    tsv_out_path = 'data_report.tsv.gz'
    fasta_written = False
    report_written = False
    with zipfile.ZipFile(genbank_zip, 'r') as zip_ref:
        for name in zip_ref.namelist():
            if name.endswith('.fna'):
                with io.TextIOWrapper(zip_ref.open(name, 'r'), encoding='utf-8') as fasta_in, lzma.open(fasta_out_path, 'wt') as fasta_out:
                    filter_genbank_sequences(fasta_in, fasta_out, min_length, max_N_proportion, fasta_out_path)
                fasta_written = True
            elif name.endswith('data_report.jsonl'):
                with zip_ref.open(name) as jsonl_in, gzip.open(tsv_out_path, 'wt') as tsv_out:
                    extract_metadata(jsonl_in, tsv_out, tsv_out_path)
                report_written = True
    if not fasta_written:
        raise ValueError(f"Failed to find .fna file in {genbank_zip}.")
    if not report_written:
        raise ValueError(f"Failed to find data_report.jsonl file in {genbank_zip}.")
    return fasta_out_path, tsv_out_path


def open_maybe_decompress(filename, mode='rt', encoding='utf-8'):
    """
    Open a file, automatically using gzip or lzma decompression based on the file extension.
    Supports .gz, .xz, and uncompressed files.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, mode, encoding=encoding)
    elif filename.endswith('.xz'):
        return lzma.open(filename, mode, encoding=encoding)
    else:
        return open(filename, mode, encoding=encoding)


def record_to_fasta_bytes(record):
    """Convert a SeqRecord to FASTA-encoded bytes (utf-8)"""
    buf = io.StringIO()
    SeqIO.write(record, buf, 'fasta')
    return buf.getvalue().encode('utf-8')


def align_sequences(refseq_fasta, extra_fasta, genbank_fasta, refseq_acc, min_length, max_N_proportion):
    """Run nextclade to align the filtered sequences to the reference, and pipe its output to faToVcf and gzip."""
    msa_vcf_gz = 'msa.vcf.gz'
    nextclade_err_txt = 'nextclade.err.txt'
    start_time = time.perf_counter()
    # Set up the pipeline: | nextclade | faToVcf | gzip > msa.vcf.gz
    nextclade_cmd = [
        'nextclade', 'run', '--input-ref', refseq_fasta, '--include-reference', 'true', '--output-fasta', '/dev/stdout'
    ]
    fatovcf_cmd = ['faToVcf', '-includeNoAltN', '-ref=' + refseq_acc, 'stdin', 'stdout']
    with gzip.open(msa_vcf_gz, 'wb') as vcf_out, open(nextclade_err_txt, 'wb') as nextclade_stderr:
        try:
            # Start nextclade
            nextclade_proc = subprocess.Popen(nextclade_cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=nextclade_stderr)
            # Start faToVcf, reading from nextclade's stdout
            fatovcf_proc = subprocess.Popen(fatovcf_cmd, stdin=nextclade_proc.stdout, stdout=subprocess.PIPE)
            nextclade_proc.stdout.close()  # Allow nextclade to receive SIGPIPE if faToVcf exits
            # Send fasta input(s) to nextclade_proc's stdin.  If extra_fasta duplicates sequences
            # also found in genbank_fasta, prefer the ones in extra_fasta.
            extra_ids = dict()
            if extra_fasta:
                with open_maybe_decompress(extra_fasta) as ef:
                    for record in SeqIO.parse(ef, 'fasta'):
                        if record.id in extra_ids:
                            print(f"Fasta ID '{record.id}' has already been used in {extra_fasta}, discarding duplicate.", file=sys.stderr)
                        else:
                            extra_ids[record.id] = True
                            if not passes_seq_filter(record, min_length, max_N_proportion):
                                continue
                            record.description = record.id
                            nextclade_proc.stdin.write(record_to_fasta_bytes(record))
            # genbank_fasta has already been filtered for length & Ns, but discard duplicates
            with open_maybe_decompress(genbank_fasta) as gf:
                for record in SeqIO.parse(gf, 'fasta'):
                    if record.id in extra_ids:
                        print(f"GenBank ID '{record.id}' was already used in {extra_fasta}, discarding duplicate.", file=sys.stderr)
                    else:
                        nextclade_proc.stdin.write(record_to_fasta_bytes(record))
            nextclade_proc.stdin.close()
            # Write faToVcf's output to gzip file
            for chunk in iter(lambda: fatovcf_proc.stdout.read(8192), b''):
                vcf_out.write(chunk)
            fatovcf_proc.stdout.close()
            fatovcf_proc.wait()
            nextclade_proc.wait()
        except subprocess.CalledProcessError as e:
            print(f"nextclade | faToVcf failed: {e.stderr}", file=sys.stderr)
            raise
    elapsed_time = time.perf_counter() - start_time
    print(f"Nextclade alignment and VCF conversion completed successfully. Wrote VCF file: {msa_vcf_gz} in {elapsed_time:.1f}s")
    # gzip nextclade.err.txt (note: using gzip.open above did not work, it was written uncompressed.)
    run_command(['gzip', '-f', nextclade_err_txt])
    return msa_vcf_gz


def make_empty_tree():
    """Create an empty tree file to be used as a starting point for usher-sampled."""
    empty_tree_path = 'empty_tree.nwk'
    with open(empty_tree_path, 'w') as f:
        f.write("()\n")
    print(f"Created empty tree file: {empty_tree_path}")
    return empty_tree_path


def run_usher_sampled(tree, vcf):
    """Build the tree.  Use docker --platform linux/amd64 so this will work even on Mac with ARM CPU. """
    pb_out = 'usher_sampled.pb.gz'
    start_time = time.perf_counter()
    command = ['usher-sampled', '-A', '-e', '5', '-t', tree, '-v', vcf, '-o', pb_out,
               '--optimization_radius', '0', '--batch_size_per_process', '100']
    run_command(command, stdout_filename='usher-sampled.out.log', stderr_filename='usher-sampled.err.log')
    elapsed_time = time.perf_counter() - start_time
    print(f"Ran usher-sampled in {elapsed_time:.1f}s")
    return pb_out


def run_matoptimize(pb_file, vcf_file):
    """Run matOptimize to clean up after usher-sampled"""
    pb_out = 'optimized.unfiltered.pb.gz'
    start_time = time.perf_counter()
    # Try with VCF, which is less tested but should give better results because it includes
    # info about which bases are ambiguous or N.  That info is lost when usher imputes values.
    command = ['matOptimize', '-m', '0.00000001', '-M', '1',
               '-i', pb_file, '-v', vcf_file, '-o', pb_out]
    if not run_command(command, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log', fail_ok=True):
        print("matOptimize with VCF failed, trying again without VCF.")
        command = ['matOptimize', '-m', '0.00000001', '-M', '1',
                   '-i', pb_file, '-o', pb_out]
        run_command(command, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log')
    elapsed_time = time.perf_counter() - start_time
    print(f"Ran matOptimize in {elapsed_time:.1f}s")
    return pb_out


def run_matutils_filter(opt_unfiltered_tree, max_parsimony, max_branch_length):
    """Run matUtils extract to filter sequences/branches and collapse tree post-matOptimize"""
    pb_out = 'optimized.pb.gz'
    start_time = time.perf_counter()
    command = ['matUtils', 'extract', '-i', opt_unfiltered_tree,
               '--max-parsimony', str(max_parsimony),
               '--max-branch-length', str(max_branch_length),
               '--collapse-tree',
               '-o', pb_out]
    run_command(command, stdout_filename='matUtils.filter.out.log', stderr_filename='matUtils.filter.err.log')
    elapsed_time = time.perf_counter() - start_time
    print(f"Ran matUtils to filter (--max-parsimony {max_parsimony} --max-branch-length {max_branch_length}) in {elapsed_time:.1f}s")
    return pb_out


def sanitize_name(name):
    """Replace characters that could cause trouble, such as spaces or Newick special characters,
    with underscores."""
    for c in ['[', ']', '(', ')', ':', ';', ',', "'", ' ']:
        name = name.replace(c, '_')
    return name


def fudge_isolate(isolate, accession, date, country, year):
    """If isolate looks like a full /-separated descriptor then use it.
    Otherwise try to approximate it from available metadata bits and pieces."""
    full = None
    if '/' in isolate:
        full = isolate
    else:
        if country and isolate and year:
            full = '/'.join([country, isolate, year])
        elif country and year:
            full = '/'.join([country, year])
        elif isolate and year:
            full = '/'.join([isolate, year])
        elif isolate:
            full = isolate
    return full


def rename_seqs(pb_in, data_report_tsv, acc_to_strain):
    """Rename sequence names in tree to include country, isolate name and date when available.
    Prepare metadata for taxonium."""
    rename_out = "rename.tsv"
    metadata_out = "metadata.tsv.gz"
    pb_out = "viz.pb.gz"
    start_time = time.perf_counter()
    with gzip.open(data_report_tsv, 'rt') as tsv_in, open(rename_out, 'w') as r_out, gzip.open(metadata_out, 'wt') as m_out:
        header = tsv_in.readline().split('\t')
        accession_idx = header.index('accession')
        isolate_idx = header.index('isolate')
        date_idx = header.index('date')
        country_idx = header.index('country')
        # No header for rename_out; write header for metadata_out
        m_out.write('strain' + '\t' + '\t'.join(header))

        # Create a mapping from accession to metadata
        acc_to_name = {}
        for line in tsv_in:
            fields = line.split('\t')
            accession = fields[accession_idx].strip()
            isolate = sanitize_name(fields[isolate_idx].strip())
            strain = sanitize_name(acc_to_strain.get(accession, ''))
            if not isolate:
                isolate = strain                
            elif strain and len(isolate) < len(strain):
                # Sometimes people use nice strain like "MVi/Marseille.FRA/21.19/20[D8]" but goofy isolate
                # like "9061923678/9051834309" for MZ031240.1.  If strain is longer, use it.
                isolate = strain
            date = fields[date_idx].strip()
            country = sanitize_name(fields[country_idx].strip())
            groups = re.match('^([0-9]{4})(-[0-9]{2})?(-[0-9]{2})?$', date)
            year = None
            if groups:
                year = groups[1]
            if not date:
                date = '?'
            full = fudge_isolate(isolate, accession, date, country, year)
            if full:
                name = '|'.join([full, accession, date])
            else:
                name = '|'.join([accession, date])
            acc_to_name[accession] = name
            r_out.write("\t".join([accession, name]) + '\n')
            m_out.write(name + '\t' + '\t'.join(fields))
    command = ['matUtils', 'mask', '-i', pb_in, '--rename-samples', rename_out, '-o', pb_out]
    run_command(command, stdout_filename='matUtils.rename.out.log', stderr_filename='matUtils.rename.err.log')
    elapsed_time = time.perf_counter() - start_time
    print(f"Ran matUtils to rename sequences for display in {elapsed_time:.1f}s")
    return rename_out, metadata_out, pb_out


def get_header(tsv_in):
    with gzip.open(tsv_in, 'rt') as tsv:
        header = tsv.readline().split('\t')
        for idx, field in enumerate(header):
            header[idx] = field.strip()
    return header


def usher_to_taxonium(pb_in, metadata_in):
    jsonl_out = "tree.jsonl.gz"
    start_time = time.perf_counter()
    columns = ','.join(get_header(metadata_in))
    command = ['usher_to_taxonium', '--input', pb_in, '--metadata', metadata_in,
               '--columns', columns, '--output', jsonl_out]
    run_command(command, stdout_filename='utt.out.log', stderr_filename='utt.err.log')
    elapsed_time = time.perf_counter() - start_time
    print(f"Ran usher_to_taxonium in {elapsed_time:.1f}s")
    return jsonl_out


def main():
    parser = argparse.ArgumentParser(prog="viral_usher_build")
    parser.add_argument("--config", type=str, required=True, help="Path to config file input")

    args = parser.parse_args()

    print(f"Running pipeline with config file: {args.config}")
    config_contents = config.parse_config(args.config)
    refseq_acc = config_contents['refseq_acc']
    assembly_id = config_contents['refseq_assembly']
    taxid = config_contents['taxonomy_id']
    min_length_proportion = float(config_contents.get('min_length_proportion', config.DEFAULT_MIN_LENGTH_PROPORTION))
    max_N_proportion = float(config_contents.get('max_N_proportion', config.DEFAULT_MAX_N_PROPORTION))
    max_parsimony = int(config_contents.get('max_parsimony', str(config.DEFAULT_MAX_PARSIMONY)))
    max_branch_length = int(config_contents.get('max_branch_length', str(config.DEFAULT_MAX_BRANCH_LENGTH)))
    extra_fasta = config_contents.get('extra_fasta', '')
    refseq_zip = f"{refseq_acc}.zip"
    genbank_zip = f"genbank_{taxid}.zip"

    # Download the RefSeq genome and all GenBank genomes for the given Taxonomy ID
    ncbi = ncbi_helper.NcbiHelper()
    print(f"Downloading RefSeq {refseq_acc} (Assembly {assembly_id}) genome to {refseq_zip}...")
    ncbi.download_refseq(assembly_id, refseq_zip)
    print(f"Downloading all GenBank genomes for taxid {taxid} to {genbank_zip}...")
    ncbi.download_genbank(taxid, genbank_zip)
    # Get strain metadata from NCBI Virus API
    print(f"Querying NCBI Virus API for extra metadata for taxid {taxid}...")
    acc_to_strain = ncbi.query_ncbi_virus_metadata(taxid)

    refseq_fasta, refseq_gbff, refseq_length = unpack_refseq_zip(refseq_zip, refseq_acc)

    min_length = int(refseq_length * min_length_proportion)
    print(f"Using min_length_proportion={min_length_proportion} ({min_length} bases) and max_N_proportion={max_N_proportion}")
    genbank_fasta, data_report = unpack_genbank_zip(genbank_zip, min_length, max_N_proportion)
    msa_vcf = align_sequences(refseq_fasta, extra_fasta, genbank_fasta, refseq_acc, min_length, max_N_proportion)
    empty_tree = make_empty_tree()
    preopt_tree = run_usher_sampled(empty_tree, msa_vcf)
    opt_unfiltered_tree = run_matoptimize(preopt_tree, msa_vcf)
    opt_tree = run_matutils_filter(opt_unfiltered_tree, max_parsimony, max_branch_length)
    rename_tsv, metadata_tsv, viz_tree = rename_seqs(opt_tree, data_report, acc_to_strain)
    usher_to_taxonium(viz_tree, metadata_tsv)


if __name__ == "__main__":
    main()
