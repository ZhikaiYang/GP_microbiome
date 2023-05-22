[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)


## Data

#### Genotype ID
- Genotype names lookup table: `data/BG_genotype_names.csv` 
 - (not needed, because the data files I provide just use the same genotype name, i.e. `GX_name` in the table)
 - TODO: need a master file (n_MM_name = 232, n_GX_name = 277, n_overlap_with_SNP_file_geno_name = 212)
 
 - Population structure
  - PCs: `data/hmp321_282_agpv4_maf005_miss03_pruned.eigenvec`

#### Genome-sequencing data
- Whole genome sequencing data (n=271) of HapMap V3.2.1 (with imputation, AGPv4) obtained from the Panzea database https://www.panzea.org/genotypes (for more information : [Bukowski et. al., 2018](https://doi.org/10.1093/gigascience/gix134))
- SNP data (20 million): ``
- SNP data (50k SNPs): `largedata/Zhikai/sinfo_micobiome_geno.txt`



#### Omics data
- Microbiome (3626 ASVs): log relative abundance (`largedata/Zhikai/sinfo_micobiome_geno.txt`), count (`largedata/Zhikai/sinfo_micobiome_geno_count.txt`) (for more information : [Meier et. al., 2022](https://doi.org/10.7554/eLife.75790))
  - Sample info: `data/sample_info_3626asvs.txt` 
  - For information about N treatment level of each ceil_id corresponding to two rows with the same genotype


#### Phenomics data
- UAV Phenotype: Canopy Coverage (`data/ppj220030-sup-0002-tables1.csv`), Vegetation Indices (`data/ppj220030-sup-0003-tables2.csv`) 
  - (for more information: [Rodene et. al., 2022](https://acsess.onlinelibrary.wiley.com/doi/full/10.1002/ppj2.20030))


## Data Processing and Sharing:

`profiling/1.A.3_yanglab_data_sharing.Rmd`


# Project Guideline

- To guide group members having a better sense about the project layout, here we briefly introduce the specific purposes of the [dir system](https://jyanglab.github.io/2017-01-07-project/). The layout of dirs is based on the idea borrowed from [ProjectTemplate](http://projecttemplate.net/architecture.html).

- The guideline for the collaborative [workflow](https://jyanglab.github.io/2017-01-10-project-using-github/).

- Check out progress and things [to-do](TODO.md) and throw ideas via the wiki page.


## License
This is an ongoing research project. It was intended for internal lab usage. It has not been extensively tested. Use at your own risk.
It is a free and open source software, licensed under [GPLv3](LICENSE).
