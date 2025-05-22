# Usher Helper

This project provides a Python script to interactively query NCBI for viral species, retrieve Taxonomy and RefSeq IDs, and download relevant genome data using the NCBI Datasets API.

## Features
- Interactive user prompts for viral species name
- Automated NCBI Taxonomy and RefSeq queries
- Download of selected RefSeqs and all GenBank genomes for a given Taxonomy ID

## Requirements
- Python 3.8+
- `requests` library

## Usage
Run the driver script:

```
python driver.py
```

Follow the prompts to select a viral species and download genome data.
