# Contacts comparison between conditions

For two conditions, from the CSV files describing the amino acids contacts between a Region of Interest and the other 
regions of a protein during a Molecular Dynamics simulation, the script will produce Multiple Sequences Alignments 
(MSA) annotated with the number of contacts.
In red, the number of contacts that are different between the two conditions and in blue the common contacts.

The input CSV data are produced by the [plot_contacts](https://github.com/njeanne/plot_contacts/tree/main) script.

## Conda environment

A [conda](https://docs.conda.io/projects/conda/en/latest/index.html) YAML environment file is provided: 
`conda_env/contacts_comparison_env.yml`.
The file contains all the dependencies to run the script.
The conda environment is generated using the command:
```shell script
# create the environment
conda env create -f conda_env/contacts_comparison_env.yml

# activate the environment
conda activate contacts_comparison
```

## Usage

The script can be tested with the test data provided in the `data` directory, which contains a CSV file describing the 
different conditions and the location of the directory containing the CSV output files from the [plot_contacts.py](https://github.com/njeanne/plot_contacts) 
script.

The input CSV file must be a comma separated file with a header as in the following example:

| condition    | path |
|--------------|---|
| insertions   | data/plot_contacts_outputs/insertions |
| duplications | data/plot_contacts_outputs/duplications |
| WT           | data/plot_contacts_outputs/WT |

An alignment in fasta format of all the sequences of the `plot_contacts.py` must be performed before.

The command:
```shell script
conda activate contacts_comparison

./contacts_comparison.py --out results --domains data/AB248520-3e_WT_ORF1_domains.csv \
--aln data/alignment_whole_sequences_AA.fa \
--md-time 1020 --roi HVR data/common_contacts_conditions_HVR_annotations_goulet.csv

conda deactivate
```

The option `--group` can be used
to regroup under the same condition `ìnsertions` and `duplications` in the test dataset with the following synthax, 
`--group insertions duplications`.


## Outputs

The outputs are:
- two CSV files by compared conditions, one for the common positions, one for the different positions.
  An example of the different positions between the equences with insertions and the Wild types ones.
 
|position alignment|number of contacts|domain                    |number of samples with contacts|number of samples|original positions                                |
|------------------|------------------|--------------------------|-------------------------------|-----------------|--------------------------------------------------|
|17                |2                 |Putative Capping Pore nsP1|1                              |9                |17:HEPAC-154_KIF1B_ORF1                           |
|51                |1                 |Putative Capping Pore nsP1|1                              |9                |51:HEPAC-154_KIF1B_ORF1                           |
|94                |1                 |Putative Capping Pore nsP1|1                              |9                |94:HEPAC-154_KIF1B_ORF1                           |
|103               |1                 |Putative Capping Pore nsP1|1                              |9                |103:HEPAC-154_KIF1B_ORF1                          |
|294               |2                 |Putative Capping Pore nsP1|2                              |9                |294:HEPAC-100_GATM_ORF1 &#124; 294:HEPAC-154_KIF1B_ORF1|
|...           |...                 |...               |...                              |...               |...                           |
|1973              |2                 |RdRp nsP5                 |1                              |9                |1738:HEPAC-100_GATM_ORF1                          |
|1975              |1                 |RdRp nsP5                 |1                              |9                |1737:HEPAC-93_EEF1A1_ORF1                         |
|1988              |2                 |after RdRp nsP5           |1                              |9                |1777:HEPAC-64_ZNF787_ORF1                         |
|1989              |1                 |after RdRp nsP5           |1                              |9                |1778:HEPAC-64_ZNF787_ORF1                         |

- the annotated MSAs by domain with the contacts between the region of interest and the other domains.
  The numbers of common contacts' positions are highlighted in blue and the number of different contacts' position in red.
  ![MSA](doc/_static/msa.svg)
