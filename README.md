# viral_usher

**Viral Usher** is a command-line tool to set up and run a pipeline to build an UShER tree for a new viral species (or type, subtype, etc.) using genomes downloaded from NCBI.

---

## 🔧 Features

- Subcommands:
  - `init`: Generate a config file (interactive or via command line options)
  - `run`: Execute the pipeline using the config file
- Conda/Bioconda integration for tool installation
- Portable to laptops, servers, or cloud platforms

---

## 📦 Installation

1. Clone the repo:
   ```bash
   git clone https://github.com/yourname/viral_usher.git
   cd viral_usher
   ```
2. Install the package:
    ```bash
    pip install .

---

## Requirements
- conda or mamba
- Python 3.11+
- Python libraries: `requests`, `biopython`

---

## 🚀 Usage

### Create a config file with `viral_usher init`
Interactive setup:

    ```bash
    viral_usher init
    ```
With command line options:

    ```bash
    viral_usher init --refseq NC_004162.2
    ```

### Run pipeline
    ```bash
    viral_usher run --input data/input.fasta
    ```

---

## 🧪 Development

### Install dev dependencies

    ```bash
    pip install -e .
    ```

### Run tests
someday
