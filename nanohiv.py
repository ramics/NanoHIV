import click

import os
import subprocess
import sys
import tempfile

# run minimap2 with nanopore settings
def map(ref, fq, gap_open, f_out):

    with open(f_out, 'w') as out:
        subprocess.call(["minimap2",
                        "-x",
                        "map-ont",
                        "-O", str(gap_open),
                        "-E", str(gap_open/2.0),
                        "-a",
                        fq], 
                        stdout=out
    )

# take a minimap2 samfile, compress, sort and index it
def prepare_bam(sam, bam):
    
    with open(bam, 'w') as out:

        # compress
        temp1, temp_name1 = tempfile.mkstemp()
        subprocess.call(["samtools",
                        "view",
                        "-b",
                        sam],
                        stdout=temp1
        )
        os.close(temp1)

        # sort
        temp2, temp_name2 = tempfile.mkstemp()
        subprocess.call(["samtools",
                        "sort",
                        temp_name1],
                        stdout=temp2
        )
        os.close(temp2)

        # index
        subprocess.call(["samtools",
                        "index",
                        temp_name2],
                        stdout=out
        )

        os.unlink(temp_name1)
        os.unlink(temp_name2)

# generate consensus sequence from bam file,
# reference and fastq
def generate_consensus(ref, fq, bam, fasta):

    temp, temp_name = tempfile.mkstemp()
    subprocess.call(["nanopolish",
                    "variants",
                    "--max-haplotypes",
                    "3000", 
                    "--faster",
                    "--fix-homopolymers",
                    "--snps",
                    "-t", 
                    "32",
                    "-p",
                    "1",
                    "--reads", fq,
                    "--bam", bam,
                    "--genome", ref],
                    stdout=temp
    )
    os.close(temp)

    with open(fasta, "w") as out:
        subprocess.call(["nanopolish",
                        "vcf2fasta",
                        "-g", ref,
                        temp_name],
                        stdout=out
        )
    os.unlink(temp_name)

@click.command()
@click.option(
    '--reference',
     type=click.Path(exists=True, readable=True, resolve_path=True)
)
@click.option(
    '--reads',
     type=click.Path(exists=True, readable=True, resolve_path=True)
)
@click.option(
    '--output',
     type=click.Path(resolve_path=True)
)
@click.option(
    '--standard-gap-penalty', default=4
)
@click.option(
    '--lower-gap-penalty', default=1
)

def nanohiv(reference, reads, output, standard_gap_penalty, lower_gap_penalty):
    """
    Check consensus sequences for intactness.
    """

    if reference is None:
        raise ValueError('Need a --reference.')
    
    if reads is None:
        raise ValueError('Need a fastq file --reads.')
    
    if output is None:
        raise ValueError('Need an --output file.')

    sam, sam_name = tempfile.mkstemp()
    bam, bam_name = tempfile.mkstemp()
    vcf, vcf_name = tempfile.mkstemp()

    # generate a consensus
    print("Generating initial consensus.")
    consensus1, consensus1_name = tempfile.mkstemp()
    map(reference, reads, standard_gap_penalty, sam_name)
    prepare_bam(sam_name, bam_name)
    generate_consensus(reference, reads, bam_name, consensus1_name)

    # remap to consensus with lower gap penalty
    print("Remapping reads to consensus with lower gap penalty.")
    consensus2, consensus2_name = tempfile.mkstemp()
    map(consensus1_name, reads, lower_gap_penalty, sam_name)
    prepare_bam(sam_name, bam_name)
    generate_consensus(consensus1_name, reads, bam_name, consensus2_name)

    # final remapping to consensus
    print("Final remapping to consensus.")
    map(consensus2_name, reads, standard_gap_penalty, sam_name)
    prepare_bam(sam_name, bam_name)
    generate_consensus(consensus2_name, reads, bam_name, output)

    os.unlink(sam_name)
    os.unlink(bam_name)
    os.unlink(vcf_name)
    os.unlink(consensus1_name)
    os.unlink(consensus2_name)

    print("Complete.")

if __name__ == '__main__':
   nanohiv()