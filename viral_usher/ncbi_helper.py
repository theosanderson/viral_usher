"""Get data from NCBI using APIs (EUtils, Virus, Datasets) or the datasets command-line tool."""

import requests
import logging
import sys
import time
import json
import subprocess
import xml.etree.ElementTree as ET

NCBI_DATASETS_V2_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_DATASETS_TAXONOMY_SUGGEST_URL = f"{NCBI_DATASETS_V2_BASE}/taxonomy/taxon_suggest"
NCBI_EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
NCBI_EUTILS_SEARCH_URL = f"{NCBI_EUTILS_BASE}/esearch.fcgi"
NCBI_EUTILS_SUMMARY_URL = f"{NCBI_EUTILS_BASE}/esummary.fcgi"
NCBI_EUTILS_ELINK_URL = f"{NCBI_EUTILS_BASE}/elink.fcgi"
NCBI_EUTILS_EFETCH_URL = f"{NCBI_EUTILS_BASE}/efetch.fcgi"


class NcbiHelper:
    max_retries = 5
    retry_delay = 60
    retry_delay_factor = 3

    def __init__(self):
        self.logger = logging.getLogger(__name__)

    def set_debug(self, debug):
        if debug:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

    def run_command(self, command):
        """Run a command using subprocess.run, check for errors."""
        self.logger.debug(f"Running command: {' '.join(command)}")
        try:
            result = subprocess.run(command, check=True, capture_output=True, text=True)
            self.logger.debug(f"Command output: {result.stdout}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"Command '{' '.join(command)}' failed with exit code {e.returncode}")
            self.logger.error(f"stderr: {e.stderr}")
            sys.exit(1)

    def query_with_retry(self, func, *args, **kwargs):
        """Call func with args and kwargs, retrying on failure."""
        time.sleep(1)
        attempt = 0
        delay = self.retry_delay
        while attempt < self.max_retries:
            attempt += 1
            try:
                return func(*args, **kwargs)
            except Exception as e:
                self.logger.error(f"Function {func.__name__} failed on attempt {attempt} with error: {e}")
            if attempt < self.max_retries:
                time.sleep(delay)
                delay *= self.retry_delay_factor
            else:
                self.logger.error(f"Function {func.__name__} failed after {self.max_retries} attempts.")
                sys.exit(1)

    def check_response(self, resp):
        if resp.status_code != 200:
            self.logger.error(f"API request failed: {resp.status_code}")
            try:
                self.logger.error(json.dumps(resp.json(), indent=2))
            except json.JSONDecodeError:
                self.logger.error(resp.text)
            return False
        return True

    def get_request(self, url, params):
        """Send a GET request, check the response, and return the JSON data.  Retry on failure.  Sleep to rate-limit."""
        def get_request_core(url, params):
            resp = requests.get(url, params=params)
            if self.check_response(resp):
                try:
                    result = resp.json()
                    return result
                except json.JSONDecodeError:
                    raise Exception(f"Expected JSON from status 200 response but got decode error from url {url} with params:\n" +
                                    json.dumps(params, indent=2) + "\nResponse:\n" + resp.text)
            else:
                raise Exception(f"Bad response from {url} with params:\n{json.dumps(params, indent=2)}")
        return self.query_with_retry(get_request_core, url, params)

    def eutils_request(self, url, params):
        """Send a request, check the response, and return the JSON data."""
        params["retmode"] = "json"
        return self.get_request(url, params)

    def esearch_request(self, db, term):
        """Use ESearch to search for term in db.  Return a list of GI numbers."""
        data = self.eutils_request(NCBI_EUTILS_SEARCH_URL, {"db": db, "term": term})
        return data.get("esearchresult", {}).get("idlist", [])

    def esearch_request_one(self, db, term):
        """Use ESearch to search for term in db.  Warn if there is not exactly one result. Return the first GI number."""
        idlist = self.esearch_request(db, term)
        if not idlist or len(idlist) == 0:
            self.logger.warning(f"No ID found for {db} term '{term}'.")
            return None
        elif len(idlist) > 1:
            self.logger.warning(f"Multiple IDs found for {db} term '{term}'. Using the first one.")
        return idlist[0]

    def elink_request(self, db, dbfrom, id):
        """Use ELink to get links from dbfrom to db for the given ID.  Return a list of linked GI numbers.  Assumes that there is only one linksetdb."""
        data = self.eutils_request(NCBI_EUTILS_ELINK_URL, {"db": db, "dbfrom": dbfrom, "id": id})
        linksets = data.get("linksets", [])
        if not linksets:
            self.logger.warning(f"No linksets found for {dbfrom} GI {id}.  data: {data}")
            return []
        if len(linksets) > 1:
            self.logger.warning(f"Multiple linksets found for {dbfrom} GI {id}. Using the first one.  data: {data}")
        return linksets[0].get("ids")

    def elink_request_one(self, db, dbfrom, id):
        """Use ELink to get links from dbfrom to db for the given ID.  Warn if there is not exactly one result. Return the first linked GI number."""
        links = self.elink_request(db, dbfrom, id)
        if not links or len(links) == 0:
            self.logger.warning(f"No {db} links found for {dbfrom} GI {id}.")
            return None
        elif len(links) > 1:
            self.logger.warning(f"Multiple {db} links found for {dbfrom} GI {id}. Using the first one.")
        return links[0]

    def esummary_request(self, db, ids):
        """Use ESummary to get summaries for the given IDs in the specified database.  Return a dictionary of summaries."""
        data = self.eutils_request(NCBI_EUTILS_SUMMARY_URL, {"db": db, "id": ",".join(ids)})
        return data.get("result", {})

    def esummary_request_one(self, db, id):
        """Use ESummary to get a summary for the given ID in the specified database.  Return the summary dictionary for id."""
        data = self.eutils_request(NCBI_EUTILS_SUMMARY_URL, {"db": db, "id": id})
        if not data or "result" not in data or id not in data["result"]:
            self.logger.warning(f"No summary found for {db} ID {id}.")
            return None
        return data.get("result", {}).get(id, {})

    def get_taxonomy_entries(self, species_name):
        """Use the NCBI Datasets API to look up the Taxonomy ID(s) for a given species name.
        Return [ { "sci_name": ..., "tax_id": ... }, ...] for all matches."""
        url = f"{NCBI_DATASETS_TAXONOMY_SUGGEST_URL}/{species_name}"
        data = self.get_request(url, {})
        return data.get("sci_name_and_ids", [])

    def get_assembly_acc_for_refseq_acc(self, refseq_acc):
        """Use the NCBI EUtils API to look up the Assembly accession for a given RefSeq accession."""
        # First use NCBI Esearch to search for the Assembly GI number for the given RefSeq GI number
        refseq_gi = self.esearch_request_one("assembly", f"{refseq_acc}[Accession]")
        if not refseq_gi:
            return None
        # Next use NCBI ELink to get the assembly GI number for the given RefSeq GI number
        assembly_gi = self.elink_request_one("assembly", "nuccore", refseq_gi)
        if not assembly_gi:
            return None
        # Now use the assembly GI number to get the assembly accession
        return self.esummary_request_one("assembly", assembly_gi).get("assemblyaccession")

    def get_refseqs_for_taxid(self, taxid):
        """Use the NCBI EUtils API to look up RefSeq IDs for a given Taxonomy ID.
        Return [ { "accession": ..., "title": ..., "strain": ..., "uid": ... }, ...] for all matches."""
        # First use NCBI Esearch to search for RefSeq GI numbers for the given Taxonomy ID
        ids = self.esearch_request("nuccore", f"txid{taxid}[Organism:exp] AND srcdb_refseq[PROP]")
        if not ids:
            return None
        # Next use NCBI Esummary to get actual RefSeq accession numbers, titles etc for the GI numbers
        summaries = self.esummary_request("nuccore", ids)
        # Reformat the results into just what we need
        refseqs = []
        for refseq_gi in ids:
            entry = summaries.get(refseq_gi, {})
            accession = entry.get("accessionversion")
            if not accession:
                self.logger.warning(f"No accession for RefSeq GI number {refseq_gi}.")
                continue
            refseqs.append({
                "accession": accession,
                "title": entry.get("title", "No title"),
                "strain": entry.get("strain", "No strain"),
                "uid": refseq_gi
            })
        return refseqs

    def get_species_taxid_for_refseq(self, refseq_id):
        """Look up the Taxonomy ID for the given RefSeq ID using the NCBI EUtils API and Assembly database."""
        refseq_gi = self.esearch_request_one("nuccore", f"{refseq_id}[Accession]")
        if not refseq_gi:
            self.logger.warning(f"No RefSeq GI number for accession {refseq_id}")
            return None, None
        nuc_record = self.esummary_request_one("nuccore", refseq_gi)
        if not nuc_record:
            self.logger.warning(f"No nuccore record for RefSeq GI number {refseq_gi}")
            return None, None
        species = nuc_record.get("organism")
        taxid = nuc_record.get("taxid")
        if not taxid:
            self.logger.warning(f"No nuccore taxid for RefSeq GI number {refseq_gi}, trying assembly database")
            assembly_record = self.esummary_request_one("assembly", refseq_gi)
            species = assembly_record.get("organism", {}).get("scientificname", species)
            taxid = assembly_record.get("taxid")
        return species, taxid

    def get_species_from_taxid(self, taxid):
        """Return the species name for a given Taxonomy ID using the NCBI EUtils API."""
        tax_record = self.esummary_request_one("taxonomy", taxid)
        if not tax_record:
            self.logger.warning(f"No taxonomy record for Taxonomy ID {taxid}")
            return None
        return tax_record.get("scientificname")

    def download_zip_file(self, url, filename):
        """Download a zip file from the given URL and save it to the specified filename."""
        headers = {"accept": "application/zip"}

        def download_zip_file_core(url, filename):
            resp = requests.get(url, headers=headers, stream=True)
            if self.check_response(resp):
                with open(filename, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                return filename
            else:
                raise Exception(f"Bad response from {url}")
        return self.query_with_retry(download_zip_file_core, url, filename)

    def post_download_zip_file(self, url, params, filename):
        """Download a zip file from the given URL and save it to the specified filename."""
        headers = {"accept": "application/zip"}

        def post_download_zip_file_core(url, params, filename):
            resp = requests.post(url, headers=headers, json=params, stream=True)
            if self.check_response(resp):
                with open(filename, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                return filename
            else:
                raise Exception(f"Bad response from {url} with params:\n{json.dumps(params, indent=2)}")
        return self.query_with_retry(post_download_zip_file_core, url, params, filename)

    def download_refseq(self, assembly_id, filename):
        """Download RefSeq using the NCBI Datasets API which requires the GCF_* assembly ID."""
        url = f"{NCBI_DATASETS_V2_BASE}/genome/accession/{assembly_id}/download?include_annotation_type=GENOME_FASTA&include_annotation_type=GENOME_GBFF&filename={filename}"
        return self.download_zip_file(url, filename)

    def download_genbank(self, taxid, filename):
        """Download all GenBank genomes for the Taxonomy ID.  Return the .zip file name."""
        # Use the NCBI datasets command-line tool because the plain API seems to get cut off more often
        command = ["datasets", "download", "virus", "genome", "taxon", str(taxid),
                   "--filename", filename, "--no-progressbar"]
        return self.run_command(command)

    def download_genbank_accessions(self, accessions, filename):
        """Download specific GenBank accessions for the Taxonomy ID.  Return the .zip file name."""
        # If it looks weird that the URL has a GET-style param even though other params are separately POSTed,
        # I'm just following the API doc as of Aug. 2025:
        # https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/rest-api/#post-/virus/genome/download
        url = f"{NCBI_DATASETS_V2_BASE}/virus/genome/download?filename={filename}"
        params = {"include_sequence": ["GENOME"],
                  "aux_report": ["NONE"],
                  "accessions": accessions}
        return self.post_download_zip_file(url, params, filename)

    def query_ncbi_virus_metadata(self, taxid, filename):
        """Query the undocumented NCBI Virus API for metadata for all sequences for taxid.
        Unlike NCBI Datasets' data_report, this includes strain and serotype.  Results written to filename are CSV."""
        url = ("https://www.ncbi.nlm.nih.gov/genomes/VirusVariation/vvsearch2/" +
               "?fq=%7B%21tag%3DSeqType_s%7DSeqType_s%3A%28%22Nucleotide%22%29"
               f"&fq=VirusLineageId_ss%3A%28{taxid}%29&q=%2A%3A%2A&cmd=download&dlfmt=csv" +
               "&fl=accession%3AAccVer_s" +
               "%2Cisolate%3AIsolate_s" +
               "%2Ccountry_location%3ACountryFull_s" +
               "%2Cdate%3ACollectionDate_s" +
               "%2Clength%3ASLen_i" +
               "%2Cbiosample%3ABioSample_s" +
               "%2Chost%3AHost_s" +
               "%2Csubmitter%3ASubmitterAffilFull_s" +
               "%2Cauthors%3AAuthors_csv" +
               "%2Cstrain%3AStrain_s" +
               "%2Cserotype%3ASerotype_s" +
               "%2Csegment%3ASegment_s"
               )
        # Request the URL, write results to filename

        def query_ncbi_virus_metadata_core(url, filename):
            resp = requests.get(url, stream=True)
            if self.check_response(resp):
                with open(filename, "wb") as f:
                    for chunk in resp.iter_content(chunk_size=8192):
                        if chunk:
                            f.write(chunk)
                return filename
            else:
                raise Exception(f"Bad response from {url}")
        return self.query_with_retry(query_ncbi_virus_metadata_core, url, filename)

    def get_species_level_taxid(self, subspecies_taxid):
        """Return the species-level Taxonomy ID for a given subspecies Taxonomy ID."""
        # EFETCH with retmode=json returns only the subspecies ID, but I need the whole record with lineage etc.
        # That requires parsing the XML response.
        def get_species_level_taxid_core(subspecies_taxid):
            resp = requests.get(NCBI_EUTILS_EFETCH_URL, params={"id": subspecies_taxid, "db": "taxonomy"})
            if self.check_response(resp):
                # Parse the XML response to find the species-level Taxonomy ID
                root = ET.fromstring(resp.content)
                for taxon in root.findall(".//Taxon"):
                    if taxon.find("Rank").text == "species":
                        print(f"Found species-level Taxonomy ID: {taxon.find('TaxId').text} for subspecies Taxonomy ID: {subspecies_taxid}")
                        return taxon.find("TaxId").text
                return subspecies_taxid
            else:
                raise Exception(f"Bad response from {NCBI_EUTILS_EFETCH_URL} with params: id={subspecies_taxid}, db=taxonomy")
        return self.query_with_retry(get_species_level_taxid_core, subspecies_taxid)
