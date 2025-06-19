![GitHub CI](https://github.com/Juke34/RAIN/actions/workflows/main.yml/badge.svg)

# RAIN - RNA Alterations Investigation using Nextflow

RAIN is a Nextflow workflow designed for epitranscriptomic analyses, enabling the detection of RNA modifications in comparison to a reference genome.
Its primary goal is to distinguish true RNA editing events from genomic variants such as SNPs, with a particular emphasis on identifying A-to-I (Adenosine-to-Inosine) editing.

<img src="doc/img/IRD.png" width="300" height="100" /> <img src="doc/img/MIVEGEC.png" width="150" height="100" />

<img src="doc/img/baargin_flowchart.jpg" width="900" height="500" />

## Table of Contents

   * [Foreword](#foreword)
   * [Flowchart](#flowchart)
   * [Installation](#installation)
      * [Nextflow](#nextflow)
      * [Container platform](#container-platform)
        * [Docker](#docker)
        * [Singularity](#singularity)  
   * [Usage and test](#usage)
   * [Parameters](#parameters)
   * [Output](#output)
   * [Author](#author-and-contributors)
   * [Contributing](#contributing)


## Foreword

...

## Flowchart

...

## Installation

The prerequisites to run the pipeline are:  

  * [Nextflow](https://www.nextflow.io/)  >= 22.04.0
  * [Docker](https://www.docker.com) or [Singularity](https://sylabs.io/singularity/)  

### Nextflow 

  * Via conda 

    <details>
      <summary>See here</summary>
      
      ```bash
      conda create -n nextflow
      conda activate nextflow
      conda install bioconda::nextflow
      ```  
    </details>

  * Manually
    <details>
      <summary>See here</summary>
      Nextflow runs on most POSIX systems (Linux, macOS, etc) and can typically be installed by running these commands:

      ```bash
      # Make sure 11 or later is installed on your computer by using the command:
      java -version
      
      # Install Nextflow by entering this command in your terminal(it creates a file nextflow in the current dir):
      curl -s https://get.nextflow.io | bash 
      
      # Add Nextflow binary to your user's PATH:
      mv nextflow ~/bin/
      # OR system-wide installation:
      # sudo mv nextflow /usr/local/bin
      ```
    </details>

### Container platform

To run the workflow you will need a container platform: docker or singularity.

### Docker

Please follow the instructions at the [Docker website](https://docs.docker.com/desktop/)

### Singularity

Please follow the instructions at the [Singularity website](https://docs.sylabs.io/guides/latest/admin-guide/installation.html)

## Usage

### Help

You can first check the available options and parameters by running:

```bash
nextflow run Juke34/RAIN -r v1.5.0 --help
```

### Profile

To run the workflow you must select a profile according to the container platform you want to use:   
- `singularity`, a profile using Singularity to run the containers
- `docker`, a profile using Docker to run the containers

The command will look like that: 

```bash
nextflow run Juke34/RAIN -r vX.X.X -profile docker <rest of paramaters>
```

Another profile is available (/!\\ actually not yet implemented):

- `slurm`, to add if your system has a slurm executor (local by default) 

The use of the `slurm` profile  will give a command like this one:

```bash
nextflow run Juke34/RAIN -r vX.X.X -profile singularity,slurm <rest of paramaters>
```

### Test

With nextflow and docker available you can run (where vX.X.X is the release version you wish to use):

```bash
nextflow run -profile docker,test Juke34/RAIN -r vX.X.X
```

Or via a clone of the repository: 

```
git clone https://github.com/Juke34/rain.git
cd rain
nextflow run -profile docker,test rain.nf
```

## Parameters

```
RAIN - RNA Alterations Investigation using Nextflow - v0.1

        Usage example:
    nextflow run rain.nf -profile docker --genome /path/to/genome.fa --annotation /path/to/annotation.gff3 --reads /path/to/reads_folder --output /path/to/output --aligner hisat2

        Parameters:
    --help                      Prints the help section

        Input sequences:
    --annotation                Path to the annotation file (GFF or GTF)
    --reads                     path to the reads file, folder or csv. If a folder is provided, all the files with proper extension in the folder will be used. You can provide remote files (commma separated list).
                                    file extension expected : <.fastq.gz>, <.fq.gz>, <.fastq>, <.fq> or <.bam>. 
                                                              for paired reads extra <_R1_001> or <_R2_001> is expected where <R> and <_001> are optional. e.g. <sample_id_1.fastq.gz>, <sample_id_R1.fastq.gz>, <sample_id_R1_001.fastq.gz>)
                                    csv input expects 6 columns: sample, fastq_1, fastq_2, strandedness and read_type. 
                                    fastq_2 is optional and can be empty. Strandedness, read_type expects same values as corresponding RAIN parameter; If a value is provided via RAIN paramter, it will override the value in the csv file.
                                    Example of csv file:
                                        sample,fastq_1,fastq_2,strandedness,read_type
                                        control1,path/to/data1.fastq.bam,,auto,short_single
                                        control2,path/to/data2_R1.fastq.gz,path/to/data2_R2.fastq.gz,auto,short_paired
    --genome                    Path to the reference genome in FASTA format.
    --read_type                 Type of reads among this list [short_paired, short_single, pacbio, ont] (no default)

        Output:
    --output                    Path to the output directory (default: result)

       Optional input:
    --aligner                   Aligner to use [default: hisat2]
    --edit_site_tool            Tool used for detecting edited sites. Default: reditools3
    --strandedness              Set the strandedness for all your input reads (default: null). In auto mode salmon will guess the library type for each fastq sample. [ 'U', 'IU', 'MU', 'OU', 'ISF', 'ISR', 'MSF', 'MSR', 'OSF', 'OSR', 'auto' ]
    --edit_threshold            Minimal number of edited reads to count a site as edited (default: 1)
    --aggregation_mode          Mode for aggregating edition counts mapped on genomic features. See documentation for details. Options are: "all" (default) or "cds_longest"
    --clipoverlap               Clip overlapping sequences in read pairs to avoid double counting. (default: false)

        Nextflow options:
    -profile                    Change the profile of nextflow both the engine and executor more details on github README [debug, test, itrop, singularity, local, docker]
```

## Output

Here the description of typical ouput you will get from RAIN:  

```
└── rain_results                                         # Output folder set using --outdir. Default: <alignment_results>
    │
    ├── AliNe                                            # Folder containing AliNe alignment pipeline result (see https://github.com/Juke34/AliNe)
    │   ├── alignment                                    # bam alignment used by RAIN
    │   ├── salmon_strandedness                          # strandedness collected by AliNe in case auto mode was in used for fastq files
    │   └── ...      
    │
    ├── bam_indicies                                     # bam and indices bam.bai
    │
    ├── FastQC                                           # bam and indices bam.bai
    │
    ├── gatk_markduplicates                              # metrics and bam after markduplicates
    │
    └── Reditools2/Reditools3/Jacusa/sapin/              # Editing output from corresponding tool
    │
    └── feature_edits                                    # Editing computed at different level (genomic features, chromosome, global)

## Author and contributors

Eduardo Ascarrunz (@eascarrunz)
Jacques Dainat  (@Juke34)

## Contributing

Contributions from the community are welcome ! See the [Contributing guidelines](https://github.com/Juke34/rain/blob/main/CONTRIBUTING.md)