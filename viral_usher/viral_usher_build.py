import argparse
import csv
import datetime
import gzip
import io
import lzma
import os
import re
import shutil
import subprocess
import sys
import time
import zipfile
from Bio import SeqIO
from . import config
from . import ncbi_helper

update_tree_input = "optimized.pb.gz"
update_nextclade_input = "nextclade.clade.tsv"


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
                print(f"{command[0]} failed: {e.stderr}\nFailed command: {' '.join(command)}", file=sys.stderr)
            elif stderr_filename:
                print(f"{command[0]} failed, see {stderr_filename}\nFailed command: {' '.join(command)}", file=sys.stderr)
            else:
                print(f"{command[0]} failed\nFailed command: {' '.join(command)}", file=sys.stderr)
            sys.exit(1)
    if stdout:
        stdout.close()
    if stderr:
        stderr.close()
    return success


def run_command_get_stdout(command):
    """Run a command and return its stdout output."""
    try:
        result = subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        return result.stdout.decode('utf-8')
    except subprocess.CalledProcessError as e:
        print(f"Command failed: {e.stderr}\nFailed command: {' '.join(command)}", file=sys.stderr)
        sys.exit(1)


def start_timing(message):
    print(message)
    return time.perf_counter()


def finish_timing(start_time):
    elapsed_time = time.perf_counter() - start_time
    print(f"... done in {elapsed_time:.1f}s")


def normalize_segment(segment):
    # Many viruses have segments S (small), M (medium), and L (large).  The vast majority of records use S, M or L
    # but there are a few oddballs that spell it out or use middle instead of medium.  Change those to S/M/L.
    if segment and len(segment) > 1:
        seg_lower = segment.lower()
        if seg_lower == 'small':
            segment = 'S'
        elif seg_lower == 'large':
            segment = 'L'
        elif seg_lower == 'medium' or seg_lower == 'middle':
            segment = 'M'
    return segment


def unpack_refseq_zip(refseq_zip, refseq_acc):
    """Extract and rename the fasta and gbff files from the RefSeq zip."""
    start_time = start_timing(f"Unpacking {refseq_zip}...")
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
    # Search for segment in refseq.gbff's record for refseq_acc
    segment = None
    with open(gbff_path, 'r') as gbff_in:
        for record in SeqIO.parse(gbff_in, 'genbank'):
            if record.id == refseq_acc:
                for feature in record.features:
                    if "segment" in feature.qualifiers:
                        segment = normalize_segment(feature.qualifiers["segment"][0])
                        print(f"Found segment {segment} for {refseq_acc} in {gbff_path}")
                        break
                break
    finish_timing(start_time)
    return fasta_path, gbff_path, length, segment


def scan_ncbi_virus_metadata(ncbi_virus_metadata):
    """Get just a couple attributes that are useful for filtering sequences: length and segment."""
    acc_to_length_segment = {}
    with open(ncbi_virus_metadata, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            accession = row.get("accession")
            length = int(row.get("length", 0))
            segment = normalize_segment(row.get("segment", ""))
            acc_to_length_segment[accession] = (length, segment)
    return acc_to_length_segment


def passes_seq_filter(record, length, max_N_proportion):
    """Fasta record passes the filter if it has at most max_N_proportion of Ns"""
    return record.seq.count('N') / length <= max_N_proportion


def filter_genbank_sequences(fasta_in, fasta_out, min_length, max_N_proportion, acc_to_length_segment, refseq_segment, fasta_out_path):
    # Filter fasta file and chop descriptions to get accession as name in nextclade
    start_time = start_timing(f"Filtering GenBank sequences by length >= {min_length} and proportion of Ns <= {max_N_proportion}...")
    record_count = 0
    passed_count = 0
    for record in SeqIO.parse(fasta_in, 'fasta'):
        record_count += 1
        (length, segment) = acc_to_length_segment.get(record.id, (0, None))
        if length == 0:
            length = len(record.seq)
        # Skip if segment is given but is different from refseq_segment or sequence is shorter than min_length
        if length < min_length or (segment and segment != refseq_segment):
            continue
        # Filter by proportion of ambiguous characters
        if not passes_seq_filter(record, length, max_N_proportion):
            continue
        passed_count += 1
        # Chop sequence name after first space to limit to accession only
        record.description = record.id
        SeqIO.write(record, fasta_out, 'fasta')
    finish_timing(start_time)
    print(f"Processed {record_count} sequences from GenBank and wrote {passed_count} sequences passing filters to {fasta_out_path}.")
    return record_count, passed_count


def unpack_genbank_zip(genbank_zip, acc_to_length_segment, min_length, max_N_proportion, refseq_segment):
    """Process the fasta and data_report.jsonl files from the GenBank zip,
    filtering the fasta by length and proportion of ambiguous characters and
    keeping only the accession part of sequence names, and extracting basic
    metadata from data_report.jsonl into much more compact TSV.
    Compress fasta with lzma (xz) and TSV with gzip"""
    fasta_out_path = 'genbank.fasta.xz'
    with zipfile.ZipFile(genbank_zip, 'r') as zip_ref:
        namelist = zip_ref.namelist()
        genomic_fna_paths = [name for name in namelist if name.endswith('genomic.fna')]
        if len(genomic_fna_paths) != 1:
            raise ValueError(f"Expected exactly one genomic.fna file in {genbank_zip}, found: {genomic_fna_paths}")
        with io.TextIOWrapper(zip_ref.open(genomic_fna_paths[0], 'r'), encoding='utf-8') as fasta_in, lzma.open(fasta_out_path, 'wt') as fasta_out:
            genbank_count, passed_count = filter_genbank_sequences(fasta_in, fasta_out, min_length, max_N_proportion, acc_to_length_segment, refseq_segment, fasta_out_path)
    return fasta_out_path, genbank_count, passed_count


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


def sanitize_name(name):
    """Replace characters that could cause trouble, such as spaces or Newick special characters,
    with underscores."""
    for c in ['[', ']', '(', ')', ':', ';', ',', "'", ' ']:
        name = name.replace(c, '_')
    return name


def align_sequences(refseq_fasta, extra_fasta, genbank_fasta, refseq_acc, min_length, max_N_proportion):
    """Run nextclade to align the filtered sequences to the reference, and pipe its output to faToVcf and gzip."""
    msa_vcf_gz = 'msa.vcf.gz'
    nextclade_err_txt = 'nextclade.align.err.log'
    start_time = start_timing(f"Aligning sequences with nextclade and converting to VCF, writing to {msa_vcf_gz}...")
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
                        record.description = sanitize_name(record.description)
                        record.id = record.description
                        if record.id in extra_ids:
                            print(f"Fasta ID '{record.id}' has already been used in {extra_fasta}, discarding duplicate.", file=sys.stderr)
                        else:
                            extra_ids[record.id] = True
                            length = len(record.seq)
                            if length < min_length or not passes_seq_filter(record, length, max_N_proportion):
                                continue
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
    finish_timing(start_time)
    start_time = start_timing(f"Compressing {nextclade_err_txt}...")
    # gzip nextclade.align.err.log (note: using gzip.open above did not work, it was written uncompressed.)
    run_command(['gzip', '-f', nextclade_err_txt])
    # Count the number of samples in the #CHROM header line of msa.vcf.gz
    aligned_count = 0
    with gzip.open(msa_vcf_gz, 'rt') as vcf_in:
        for line in vcf_in:
            if line.startswith('#CHROM'):
                header = line.strip().split('\t')
                aligned_count = len(header) - 9  # Subtract fixed columns
                break
    finish_timing(start_time)
    return msa_vcf_gz, aligned_count


def make_empty_tree():
    """Create an empty tree file to be used as a starting point for usher-sampled."""
    empty_tree_path = 'empty_tree.nwk'
    with open(empty_tree_path, 'w') as f:
        f.write("()\n")
    return empty_tree_path


def run_usher_sampled(tree, vcf):
    """Build the tree.  Use docker --platform linux/amd64 so this will work even on Mac with ARM CPU. """
    pb_out = 'usher_sampled.pb.gz'
    start_time = start_timing(f"Running usher-sampled on {vcf}...")
    tree_flag = "-t" if tree.endswith('.nwk') else "-i"
    command = ['usher-sampled', '-A', '-e', '5', tree_flag, tree, '-v', vcf, '-o', pb_out,
               '--optimization_radius', '0', '--batch_size_per_process', '100']
    run_command(command, stdout_filename='usher-sampled.out.log', stderr_filename='usher-sampled.err.log')
    finish_timing(start_time)
    return pb_out


def run_matoptimize(pb_file, vcf_file, update):
    """Run matOptimize to clean up after usher-sampled"""
    pb_out = 'optimized.unfiltered.pb.gz'
    start_time = start_timing(f"Running matOptimize on {pb_file}...")
    command_no_vcf = ['matOptimize', '-m', '0.00000001', '-M', '1', '-i', pb_file, '-o', pb_out]
    # Unless we're doing an update, in which case the VCF covers only the new sequences,
    # first try with VCF, which is less tested but should give better results because it includes
    # info about which bases are ambiguous or N.  That info is lost when usher imputes values.
    if update:
        run_command(command_no_vcf, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log')
    else:
        command_vcf = command_no_vcf.copy()
        command_vcf.append("-v")
        command_vcf.append(vcf_file)
        if not run_command(command_vcf, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log', fail_ok=True):
            finish_timing(start_time)
            start_time = start_timing("matOptimize with VCF failed, trying again without VCF...")
            run_command(command_no_vcf, stdout_filename='matOptimize.out.log', stderr_filename='matOptimize.err.log')
    finish_timing(start_time)
    return pb_out


def run_matutils_filter(opt_unfiltered_tree, max_parsimony, max_branch_length):
    """Run matUtils extract to filter sequences/branches and collapse tree post-matOptimize"""
    pb_out = update_tree_input
    sample_names_out = 'tree_samples.txt'
    start_time = start_timing(f"Running matUtils to filter (--max-parsimony {max_parsimony} --max-branch-length {max_branch_length})...")
    command = ['matUtils', 'extract', '-i', opt_unfiltered_tree,
               '--max-parsimony', str(max_parsimony),
               '--max-branch-length', str(max_branch_length),
               '--collapse-tree',
               '--used-samples', sample_names_out,
               '-o', pb_out]
    run_command(command, stdout_filename='matUtils.filter.out.log', stderr_filename='matUtils.filter.err.log')
    finish_timing(start_time)
    # Run matUtils summary to get the final number of samples in the tree
    start_time = start_timing("Counting samples in the final tree...")
    command = ['matUtils', 'summary', '-i', pb_out]
    summary_output = run_command_get_stdout(command)
    for line in summary_output.split('\n'):
        if line.startswith('Total Samples in Tree:'):
            tree_tip_count = int(line.split(':')[1].strip())
            break
    finish_timing(start_time)
    return pb_out, sample_names_out, tree_tip_count


def get_existing_nextclade_assignments(update, starting_tree_accessions, nextclade_path):
    """Get existing Nextclade assignments, and the number of columns to expect, from the old TSV file about to be replaced."""
    if not update or not nextclade_path:
        return {}, None
    existing_assignments = {}
    nextclade_tsv = update_nextclade_input if os.path.exists(update_nextclade_input) else update_nextclade_input + ".gz"
    existing_column_count = None
    if os.path.exists(nextclade_tsv):
        with open_maybe_decompress(nextclade_tsv) as tsv_in:
            header_cols = tsv_in.readline().rstrip('\n').split('\t')
            # 'dynamic' might have added clade and/or phenotype, motif etc. columns; keep all of them
            col_count = len(header_cols) - 1
            if existing_column_count is None:
                existing_column_count = col_count
            elif existing_column_count != col_count:
                raise ValueError(f"Inconsistent Nextclade column count in {nextclade_tsv}: {existing_column_count} vs {col_count}")
            for line in tsv_in:
                fields = line.rstrip('\n').split('\t')
                seq_name = fields[0]
                clades = fields[1:]
                existing_assignments[seq_name] = clades
    else:
        print(f"Nextclade TSV file {nextclade_tsv} not found. Only new items will have clade assignments.", file=sys.stderr)
    return existing_assignments, existing_column_count


def run_nextclade(nextclade_path, nextclade_clade_columns, existing_nextclade_assignments, existing_column_count, genbank_fasta, extra_fasta, refseq_fasta):
    """Run nextclade to assign sequences to clades.  Return dict mapping accession to clade."""
    if not nextclade_path:
        return {}, ""
    nextclade_tsv = update_nextclade_input
    # Nextclade --output-columns-selection does not accept individual customClades names; instead use 'dynamic'
    # to get all of them.  See https://github.com/nextstrain/nextclade/issues/1667 .
    nextclade_col_list = nextclade_clade_columns.split(',')
    if len(nextclade_col_list) > 1:
        output_columns = 'seqName,clade,dynamic'
    else:
        output_columns = 'seqName,clade'
    start_time = start_timing(f"Running nextclade to assign clades (dataset: {nextclade_path})...")
    command = ['nextclade', 'run',
               '--dataset-name', nextclade_path,
               '--output-columns-selection', output_columns,
               '--output-tsv', nextclade_tsv,
               genbank_fasta, refseq_fasta]
    if extra_fasta:
        command.append(extra_fasta)
    run_command(command, stdout_filename='nextclade.clade.out.log', stderr_filename='nextclade.clade.err.log')
    finish_timing(start_time)
    # Read the TSV file and return a dict mapping sequence names to clades
    start_time = start_timing(f"Reading nextclade assignments from {nextclade_tsv} to add to metadata...")
    new_assignments = {}
    with open(nextclade_tsv, 'r') as tsv_in:
        header_cols = tsv_in.readline().rstrip('\n').split('\t')
        # 'dynamic' might have added clade and/or phenotype, motif etc. columns; keep all of them
        nextclade_clade_columns = ','.join(header_cols[1:])
        nextclade_column_count = len(header_cols) - 1
        if existing_column_count is not None and nextclade_column_count != existing_column_count:
            print(f"Existing nextclade annotations had {existing_column_count} columns but new annotations have {nextclade_column_count} columns ({header_cols[1:]})." +
                  f"\nRemove {update_nextclade_input} and run again, to rebuild with the latest nextclade dataset '{nextclade_path}'", file=sys.stderr)
            sys.exit(1)
        for line in tsv_in:
            fields = line.rstrip('\n').split('\t')
            seq_name = fields[0]
            clades = fields[1:]
            new_assignments[seq_name] = clades
    if existing_column_count is not None:
        nextclade_assignments = existing_nextclade_assignments
        nextclade_assignments.update(new_assignments)
    else:
        nextclade_assignments = new_assignments
    finish_timing(start_time)
    if existing_column_count is not None:
        # Append the existing assignments to the new little update nextclade.clade.tsv so they're not lost.
        with open(nextclade_tsv, 'a') as tsv_out:
            for acc, clade_columns in existing_nextclade_assignments.items():
                if acc not in new_assignments:
                    print(acc + "\t" + "\t".join(clade_columns), file=tsv_out)
    return nextclade_assignments, nextclade_clade_columns


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


def make_display_name(isolate, strain, date, accession, country):
    """Make a display name with whatever info is available in isolate or strain, the country and date in addition to accession."""
    if not isolate:
        isolate = strain
    elif strain and len(isolate) < len(strain):
        # Sometimes people use nice strain like "MVi/Marseille.FRA/21.19/20[D8]" but goofy isolate
        # like "9061923678/9051834309" for MZ031240.1.  If strain is longer, use it.
        isolate = strain
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
    return name


def finalize_metadata(ncbi_virus_metadata, nextclade_assignments, nextclade_clade_columns, refseq_segment, sample_names):
    """Make a file for renaming sequence names in tree to include country, isolate name and date when available.
    Prepare metadata for taxonium."""
    rename_out = "rename.tsv"
    metadata_out = "metadata.tsv.gz"
    start_time = start_timing(f"Finalizing sequence names for display and metadata, writing to {metadata_out}...")
    remove_segment = not refseq_segment
    # Use a subprocess to get input from grep -Fwf {sample_names} {ncbi_virus_metadata}
    grep_cmd = ['grep', '-Fwf', sample_names, ncbi_virus_metadata]
    with subprocess.Popen(grep_cmd, stdout=subprocess.PIPE, text=True) as grep_proc, \
            open(rename_out, 'w') as r_out, gzip.open(metadata_out, 'wt') as m_out:
        # Read header from ncbi_virus_metadata (not filtered by grep)
        with open(ncbi_virus_metadata, 'rt') as csv_in:
            header = csv_in.readline().rstrip('\n').split(',')
        accession_idx = header.index('accession')
        isolate_idx = header.index('isolate')
        strain_idx = header.index('strain')
        date_idx = header.index('date')
        country_loc_idx = header.index('country_location')
        # In metadata header, make separate columns for country and location:
        header_mod = header[:country_loc_idx] + ['country', 'location'] + header[country_loc_idx+1:]
        segment_idx = header_mod.index('segment')
        if remove_segment:
            header_mod.pop(segment_idx)
        # No header for rename_out; write header for metadata_out
        metadata_header = 'strain\t' + '\t'.join(header_mod)
        if nextclade_clade_columns:
            nextclade_col_list = nextclade_clade_columns.split(',')
            nextclade_col_list = [col if col.startswith('nextclade_') or col.startswith('Nextclade)_') else 'nextclade_' + col for col in nextclade_col_list]
            metadata_header += '\t' + '\t'.join(nextclade_col_list)
        m_out.write(metadata_header + '\n')

        # Process filtered lines from grep
        reader = csv.reader(grep_proc.stdout, delimiter=',')
        for row in reader:
            accession = row[accession_idx].strip()
            isolate = sanitize_name(row[isolate_idx].strip())
            strain = sanitize_name(row[strain_idx].strip())
            date = row[date_idx].strip()
            country_loc = row[country_loc_idx].strip()
            # Make separate columns for country and location
            country, location = country_loc.split(':', 1) if ':' in country_loc else (country_loc, '')
            country = sanitize_name(country.strip())
            location = location.strip()
            row_mod = row[:country_loc_idx] + [country, location] + row[country_loc_idx+1:]
            if remove_segment:
                row_mod.pop(segment_idx)
            name = make_display_name(isolate, strain, date, accession, country)
            metadata = name + '\t' + '\t'.join(row_mod)
            if nextclade_clade_columns:
                clades = nextclade_assignments.get(accession, [''] * len(nextclade_col_list))
                metadata += '\t' + '\t'.join(clades)
            # Write output to both renaming file and metadata file
            r_out.write("\t".join([accession, name]) + '\n')
            m_out.write(metadata + '\n')
    finish_timing(start_time)
    return metadata_out, rename_out


def rename_seqs(pb_in, rename_tsv):
    pb_out = "viz.pb.gz"
    start_time = start_timing(f"Renaming sequences in {pb_in} to make {pb_out}...")
    command = ['matUtils', 'mask', '-i', pb_in, '--rename-samples', rename_tsv, '-o', pb_out]
    run_command(command, stdout_filename='matUtils.rename.out.log', stderr_filename='matUtils.rename.err.log')
    finish_timing(start_time)
    return pb_out


def dump_newick(pb_in):
    """Dump the phylogenetic tree in Newick format."""
    newick_out = "viz.nwk"
    start_time = start_timing(f"Writing tree in Newick format to {newick_out}.gz...")
    command = ['matUtils', 'extract', '-i', pb_in, '--write-tree', newick_out]
    run_command(command, stdout_filename='/dev/null', stderr_filename='/dev/null')
    command = ['gzip', '-f', newick_out]
    run_command(command, stdout_filename='/dev/null', stderr_filename='/dev/null')
    finish_timing(start_time)


def get_header(tsv_in):
    with gzip.open(tsv_in, 'rt') as tsv:
        header = tsv.readline().split('\t')
        for idx, field in enumerate(header):
            header[idx] = field.strip()
    return header


def usher_to_taxonium(pb_in, metadata_in, refseq_gbff):
    jsonl_out = "tree.jsonl.gz"
    start_time = start_timing(f"Running usher_to_taxonium to make {jsonl_out}...")
    columns = ','.join(get_header(metadata_in))
    command = ['usher_to_taxonium', '--input', pb_in, '--metadata', metadata_in,
               '--columns', columns, '--genbank', refseq_gbff, '--output', jsonl_out]
    if not run_command(command, stdout_filename='utt.out.log', stderr_filename='utt.err.log', fail_ok=True):
        finish_timing(start_time)
        start_time = start_timing("usher_to_taxonium failed with --genbank, trying again without --genbank...")
        command = ['usher_to_taxonium', '--input', pb_in, '--metadata', metadata_in,
                   '--columns', columns, '--output', jsonl_out]
        run_command(command, stdout_filename='utt.out.log', stderr_filename='utt.err.log')
    finish_timing(start_time)
    return jsonl_out


def write_output_stats(refseq_acc, refseq_length, gb_count, filtered_count, aligned_count, tree_tip_count):
    output_stats = [["refseq_acc", "refseq_length", "gb_count", "filtered_count", "aligned_count", "tree_tip_count"],
                    [refseq_acc, refseq_length, gb_count, filtered_count, aligned_count, tree_tip_count]]
    print("Writing output stats to output_stats.tsv.")
    with open("output_stats.tsv", "w") as f:
        for row in output_stats:
            f.write("\t".join(map(str, row)) + "\n")


def get_tree_names(pb_in):
    """Return a set containing the names of the tips in pb_in."""
    tree_names = set()
    start_time = start_timing(f"Extracting tree names from {pb_in}...")
    command = ["matUtils", "extract", "-i", pb_in, "--used-samples", "tree_samples.txt"]
    run_command(command, stdout_filename="/dev/null", stderr_filename="/dev/null")
    with open("tree_samples.txt", "r") as f:
        for line in f:
            tree_names.add(line.strip())
    finish_timing(start_time)
    return tree_names


def find_new_accessions(tree_names, acc_to_length_segment, min_length, refseq_segment):
    """Return a list of GenBank accessions in acc_to_length_segment that are not in the current tree and probably wouldn't be filtered later."""
    all_accessions = set(acc_to_length_segment.keys())
    # Exclude accessions whose metadata has a segment that is different from RefSeq
    exclude_accessions = {acc for acc, (length, seg) in acc_to_length_segment.items() if length < min_length or seg and seg != refseq_segment}
    new_accessions = all_accessions - tree_names - exclude_accessions
    return list(new_accessions)


def get_new_extra_fasta(tree_names, extra_fasta):
    """Get a new extra_fasta file containing only sequences not already in the tree."""
    if not extra_fasta:
        return ""
    now = datetime.datetime.now()
    new_extra_fasta = "new_extra." + now.strftime("%Y-%m-%d_%H-%M-%S") + ".fasta"
    start_time = start_timing(f"Filtering {extra_fasta} to get new sequences not already in tree, writing to {new_extra_fasta}...")
    new_count = 0
    with open_maybe_decompress(extra_fasta) as ef, open(new_extra_fasta, 'w') as nef:
        for record in SeqIO.parse(ef, 'fasta'):
            record.description = sanitize_name(record.description)
            record.id = record.description
            if record.id not in tree_names:
                SeqIO.write(record, nef, 'fasta')
                new_count += 1
    finish_timing(start_time)
    if new_count == 0:
        print(f"No new sequences found in {extra_fasta}, skipping it.")
        return ""
    else:
        print(f"Wrote {new_count} new sequences to {new_extra_fasta}.")
        return new_extra_fasta


def run_nextclade_on_existing(ncbi, starting_tree_accessions, nextclade_path, nextclade_clade_columns,
                              extra_fasta, refseq_fasta):
    """Download GenBank sequences that are already in the tree.  Run nextclade on those plus extra_fasta."""
    old_genbank_zip = "genbank_old.zip"
    start_time = start_timing("Downloading and unpacking GenBank genomes that were already in the tree, for nextclade...")
    ncbi.download_genbank_accessions(list(starting_tree_accessions), old_genbank_zip)
    old_genbank_fasta, _, _ = unpack_genbank_zip(old_genbank_zip, {}, 0, 1.0, "")
    finish_timing(start_time)
    print("Running nextclade on genomes that were already in the tree.")
    run_nextclade(nextclade_path, nextclade_clade_columns, {}, None, old_genbank_fasta, extra_fasta, refseq_fasta)
    os.remove(old_genbank_zip)
    if not os.path.exists(update_nextclade_input):
        print(f"Expected {update_nextclade_input} to be created but it's not found.", sys.stderr)
        sys.exit(1)


def main():
    parser = argparse.ArgumentParser(prog="viral_usher_build")
    parser.add_argument("--config", type=str, required=True, help="Path to config file input")
    parser.add_argument("--update", action="store_true",
                        help="Add only new sequences to existing tree instead of building tree from scratch")

    args = parser.parse_args()

    print(f"Running pipeline with config file: {args.config}")
    config_contents = config.parse_config(args.config)
    refseq_acc = config_contents['refseq_acc']
    assembly_id = config_contents['refseq_assembly']
    taxid = config_contents['taxonomy_id']
    nextclade_path = config_contents.get('nextclade_dataset', '')
    nextclade_clade_columns = config_contents.get('nextclade_clade_columns', 'clade')
    min_length_proportion = float(config_contents.get('min_length_proportion', config.DEFAULT_MIN_LENGTH_PROPORTION))
    max_N_proportion = float(config_contents.get('max_N_proportion', config.DEFAULT_MAX_N_PROPORTION))
    max_parsimony = int(config_contents.get('max_parsimony', str(config.DEFAULT_MAX_PARSIMONY)))
    max_branch_length = int(config_contents.get('max_branch_length', str(config.DEFAULT_MAX_BRANCH_LENGTH)))
    extra_fasta = config_contents.get('extra_fasta', '')
    refseq_zip = f"{refseq_acc}.zip"
    genbank_zip = f"genbank_{taxid}.zip"

    # Download the RefSeq genome and annotations
    ncbi = ncbi_helper.NcbiHelper()
    start_time = start_timing(f"Downloading RefSeq {refseq_acc} (Assembly {assembly_id}) genome to {refseq_zip}...")
    ncbi.download_refseq(assembly_id, refseq_zip)
    finish_timing(start_time)
    refseq_fasta, refseq_gbff, refseq_length, refseq_segment = unpack_refseq_zip(refseq_zip, refseq_acc)
    min_length = int(refseq_length * min_length_proportion)

    starting_tree_accessions = []
    if (args.update):
        # Make sure we have required inputs
        if not os.path.exists(update_tree_input):
            print(f"No existing {update_tree_input} found, cannot do update.  Run without --update to build tree from scratch.", file=sys.stderr)
            sys.exit(1)
        starting_tree_accessions = get_tree_names(update_tree_input)
        if nextclade_path and not os.path.exists(update_nextclade_input) and not os.path.exists(update_nextclade_input + ".gz"):
            print(f"No existing {update_nextclade_input} found; recreating it.")
            run_nextclade_on_existing(ncbi, starting_tree_accessions, nextclade_path, nextclade_clade_columns,
                                      extra_fasta, refseq_fasta)

    # Get metadata from NCBI Virus API
    start_time = start_timing(f"Querying NCBI Virus API for metadata for all GenBank sequences for taxid {taxid}...")
    ncbi_virus_metadata = "ncbi_virus_metadata.csv"
    ncbi.query_ncbi_virus_metadata(taxid, ncbi_virus_metadata)
    acc_to_length_segment = scan_ncbi_virus_metadata(ncbi_virus_metadata)
    finish_timing(start_time)

    # Download GenBank genomes for the given Taxonomy ID: only the new ones if --update, otherwise all of them
    if (args.update):
        new_accessions = find_new_accessions(starting_tree_accessions, acc_to_length_segment, min_length, refseq_segment)
        start_time = start_timing(f"Downloading {len(new_accessions)} new GenBank genomes for taxid {taxid} to {genbank_zip}...")
        ncbi.download_genbank_accessions(new_accessions, genbank_zip)
        finish_timing(start_time)
        extra_fasta = get_new_extra_fasta(starting_tree_accessions, extra_fasta)
    else:
        start_time = start_timing(f"Downloading all GenBank genomes for taxid {taxid} to {genbank_zip}...")
        ncbi.download_genbank(taxid, genbank_zip)
        finish_timing(start_time)
    genbank_fasta, gb_count, filtered_count = unpack_genbank_zip(genbank_zip, acc_to_length_segment, min_length, max_N_proportion, refseq_segment)

    # The core of the pipeline: align sequences, build tree, finalize metadata
    msa_vcf, aligned_count = align_sequences(refseq_fasta, extra_fasta, genbank_fasta, refseq_acc, min_length, max_N_proportion)
    if args.update:
        if aligned_count == 0:
            preopt_tree = "usher_sampled.pb.gz"
            shutil.copy(update_tree_input, preopt_tree)
        else:
            preopt_tree = run_usher_sampled(update_tree_input, msa_vcf)
    else:
        empty_tree = make_empty_tree()
        preopt_tree = run_usher_sampled(empty_tree, msa_vcf)
    opt_unfiltered_tree = run_matoptimize(preopt_tree, msa_vcf, args.update)
    opt_tree, sample_names, tree_tip_count = run_matutils_filter(opt_unfiltered_tree, max_parsimony, max_branch_length)
    existing_nextclade_assignments, existing_column_count = get_existing_nextclade_assignments(args.update, starting_tree_accessions, nextclade_path)
    nextclade_assignments, nextclade_clade_columns = run_nextclade(nextclade_path, nextclade_clade_columns,
                                                                   existing_nextclade_assignments, existing_column_count,
                                                                   genbank_fasta, extra_fasta, refseq_fasta)
    metadata_tsv, rename_tsv = finalize_metadata(ncbi_virus_metadata, nextclade_assignments, nextclade_clade_columns,
                                                 refseq_segment, sample_names)
    viz_tree = rename_seqs(opt_tree, rename_tsv)
    dump_newick(viz_tree)
    write_output_stats(refseq_acc, refseq_length, gb_count, filtered_count, aligned_count, tree_tip_count)
    usher_to_taxonium(viz_tree, metadata_tsv, refseq_gbff)


if __name__ == "__main__":
    main()
