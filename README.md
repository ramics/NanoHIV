# NanoHIV

Build an HIV consensus sequence from Oxford Nanopore reads.

NanoHIV maps reads to a reference, polishes a consensus against the raw signal, then remaps with a lower gap penalty so real indels survive. Built for sequencing the HIV proviral reservoir, and published in *Cells* (2021).

## Requirements

Python 3 with `click`, plus `minimap2`, `samtools` and `nanopolish`:

```
sudo apt install python3 python3-click samtools nanopolish minimap2
```

## Run

A single script, run in place:

```
python nanohiv.py --reference ref.fasta --reads reads.fastq --fast5 fast5_dir --output consensus.fasta
```

| Option | Default | Purpose |
|---|---|---|
| `--standard-gap-penalty` | 4 | gap-open penalty for the first mapping |
| `--lower-gap-penalty` | 3 | gap-open penalty for the remap |
