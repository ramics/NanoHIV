# NanoHIV

## Installation

There is no installation for NanoHIV - it can be invoked directly on the command line using Python 3. Dependencies can be installed using `apt` or `yum`, e.g.:

```
sudo apt install python3 samtools nanopolish minimap2  

```

## Usage

```
Usage: nanohiv.py [OPTIONS]

  Create an HIV consensus sequence from Oxford Nanopore data.

Options:
  --reference PATH
  --reads PATH
  --fast5 PATH
  --output PATH
  --standard-gap-penalty INTEGER
  --lower-gap-penalty INTEGER
  --help                          Show this message and exit.
```
