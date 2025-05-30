# viral_usher

**Viral Usher** is a command-line tool to set up and run a pipeline to build an UShER tree for a new viral species (or type, subtype, etc.) using genomes downloaded from NCBI.

---

## 🔧 Features

- Subcommands:
  - `init`: Generate a config file (interactive or via command line options)
  - `build`: Download sequences and build a tree, guided by the config file
- Uses Docker for portability to laptops, servers, or cloud platforms

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
    ```

---

## Requirements
- Docker
- Python 3.11+
- Python libraries: `docker`

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

### Build a tree using config file:
    ```bash
    viral_usher build --config chikungunya/config.toml
    ```

---

## 🧪 Development

### Install dev dependencies

    ```bash
    pip install -e .
    ```

### Run tests
someday
