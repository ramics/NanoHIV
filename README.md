# NanoHIV

**Accurate near full-length HIV-1 proviral genomes from Oxford Nanopore reads.**

Nanopore sequencing is a cheap way to read a 9 kb HIV amplicon end to end, and its error rate collides with HIV's own diversity: map noisy reads to a distant reference and real variation gets lost. NanoHIV solves this by iterating. Each round maps the reads, polishes a consensus from the raw signal, and uses that consensus as the next round's reference, so the reference converges on the genome actually in the tube.

Across nine single-genome amplicons sequenced on both platforms, NanoHIV consensus sequences agreed with Illumina at 99.4% on average.

> Wright IA, Delaney KE, Katusiime MGK, Botha JC, Engelbrecht S, Kearney MF, van Zyl GU.
> **NanoHIV: A Bioinformatics Pipeline for Producing Accurate, Near Full-Length HIV Proviral Genomes Sequenced Using the Oxford Nanopore Technology.**
> *Cells* 10(10), 2577 (2021). [doi:10.3390/cells10102577](https://doi.org/10.3390/cells10102577)

## How it works

| Round | Map reads to | Gap-open penalty | Output |
|---|---|---|---|
| 1 | Your reference | standard (4) | Consensus 1 |
| 2 | Consensus 1 | lowered (3), to recover true indels | Consensus 2 |
| 3 | Consensus 2 | standard (4) | Final consensus |

Mapping uses `minimap2`; consensus calling uses `nanopolish` on the raw FAST5 signal.

## Run

NanoHIV is a single script. Install its dependencies:

```
sudo apt install python3 python3-click samtools minimap2 nanopolish
```

Then:

```
python nanohiv.py \
    --reference HXB2.fasta \
    --reads sample.fastq \
    --fast5 fast5_dir/ \
    --output sample_consensus.fasta
```

| Option | Meaning |
|---|---|
| `--reference` | Starting reference genome (FASTA) |
| `--reads` | Basecalled reads for one amplicon (FASTQ) |
| `--fast5` | Directory of raw signal files for those reads |
| `--output` | Where to write the final consensus (FASTA) |
| `--standard-gap-penalty` | Gap-open penalty for rounds 1 and 3 (default 4) |
| `--lower-gap-penalty` | Gap-open penalty for round 2 (default 3) |

## Pairs with

Feed the consensus genomes to [HIVIntact](https://github.com/ramics/HIVIntact) to classify each provirus as intact or defective.
