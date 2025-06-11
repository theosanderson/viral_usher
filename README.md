# viral_usher

**Viral Usher** is a command-line tool to set up and run a pipeline to build an [UShER](https://usher-wiki.readthedocs.io/en/latest/) tree for a new viral species (or type, subtype, etc.) using genomes downloaded from [NCBI](https://www.ncbi.nlm.nih.gov/).

---

## 🔧 Features

- Subcommands:
  - `init`: Generate a config file (interactive or via command line options)
  - `build`: Download sequences and build a tree, guided by the config file
- Uses [Docker](https://www.docker.com/) for portability to laptops, servers, or cloud platforms

---

## 📦 Installation

0. Install prerequisites (if not already installed)
- [Docker](https://www.docker.com/)
- [Python](https://www.python.org/) version 3.11 or later (we highly recommend using an environment manager such as [venv](https://docs.python.org/3/library/venv.html), [miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main), [mamba](https://mamba.readthedocs.io/en/latest/installation/mamba-installation.html) etc.)

1. Clone the repo:
   ```bash
   git clone https://github.com/AngieHinrichs/viral_usher.git
   cd viral_usher
   ```
2. Install the package (again, we highly recommended using an environment manager):
    ```bash
    pip install .
    ```

---

## 🚀 Usage

### Create a config file with `viral_usher init`
For interactive setup, with the script prompting you to make a series of selections:
   ```bash
   viral_usher init
   ```
Or with command line options:
   ```bash
   viral_usher init --refseq NC_004162.2 --workdir chikungunya --config chikungunya/config.toml
   ```

### Build a tree using config file with `viral_usher build`:
   ```bash
   viral_usher build --config chikungunya/config.toml
   ```

That's all!  viral_usher will create the following files in workdir:
- a tree in UShER protobuf format (optimized.pb.gz)
- a metadata file in TSV format (metadata.tsv.gz)
- a Taxonium tree file that you can view using https://taxonium.org/ (tree.jsonl.gz)

---

## 🧪 Development

### Install dev dependencies

```bash
pip install -e .
```

### Run tests
someday
