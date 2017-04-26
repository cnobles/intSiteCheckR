# intSiteChecker: summarize specimen clonal abundance
Generate a short summary from any processed specimen for the most abundant and/or specific clones. 

## Installation
intSiteChecker runs completely in R environments and can be installed by cloning from this repository. Further, it depends on a number of other packages including: argparse, yaml, pander, RMySQL, GenomicRanges, dplyr, hiAnnotator, sonicLength, and gintools ("https://github.com/cnobles/gintools"). 
```
git clone https://github.com/cnobles/intSiteChecker.git
```

## Brief Overview
Using a specimen or list of specimens, intSiteChecker downloads all unique integration site data from a specified INSPIIRED database (either MySQL and/or SQLite, modify in configuration file). The data is then processed by breakpoint refinement methods and integration standardization methods before calculationg clonal abundances. The data is then filtered to the most abundant clones and returned in a short table to the user.

## Using intSiteChecker
There are various features to intSiteChecker than allow it to serve multiple rolls. By default, given a single specimen number, intSiteChecker will return the top 10 most abundant clones in a short summary table, as show below:
```
Rscript ~/path/to/intSiteChecker.R GTSP0927
...
%
% pPatient - m??? - PBMC - GTSP0927
Most Abundant Clones within Specimen GTSP0927

---------------------------------------------------
   gene_id         posID       estAbund   relAbund
------------- --------------- ---------- ----------
   ARK5 *~     chr19-5937498     578       0.076

   CCT4 *     chr3-120521420     492       0.065

  SHELUP3 ~    chr3+2926709      444       0.059

   BRC1 ~     chr18-47291645     309       0.041

  ATCK3 *~    chr21+65492731     304       0.040
---------------------------------------------------
*Note: This data has been scrambled from the original output, only used as an example.
```
Further, several specimen numbers can be given to intSiteChecker, which will return a series of tables with the most abundant clones from each specimen. 
```
Rscript ~/path/to/intSiteChecker.R GTSP0397 GTSP0398 GTSP0399
```
Additionally, a position id (in the format of chr#[+-]position, ie. chr12-4270498) or list of position ids and be added to make sure integration sites matching the position ids or within a "fuzzy window" from a position id will be tracked in the data summary. If the integration sites given in position ids are the only desired results to return, then they can be exclusively processed, through this will remove the relative abundance of the sites with respect to the sample. 
```
Rscript ~/path/to/intSiteChecker.R GTSP0397 -p chr19-28833957
Rscript ~/path/to/intSiteChecker.R GTSP0397 -p chr19-28833957 chr17+2825892 
Rscript ~/path/to/intSiteChecker.R GTSP0397 -p chr19-28833957 chr17+2825892 -e
```

## Arguments
**specimen** [positional] specifies which specimen(s) to summarize. Multiple specimens can be separated by a space.

**help** [-h, --help] returns usage and argument infomation.

**position ids** [-p, --position_ids] specifies which genomic location(s) to focus and return in the summary along with most abundant clones. Multiple position ids can be separated by a space after the "-p" flag.

**position ids exclusively** [-e, --position_ids_exclusively] If included, analysis will only focus on genomic locations specified by position ids and will drop the rest of the integration sites of the specimen(s).

**abundance method** [-a, --abundance_method] Method for calculating clonal abundance. Defaults to config option. Available options include: readCount, uniqFragLength, sonicAbundance.

**number of top clones** [-n, --num_of_clones] Number of most abundant clones to return in summary.

**fuzzy window** [-f, --fuzzy_window] A number to allow for +/- differences in position ids. Given 5, the genomic locations specified by position ids will be given a +/- 5 bp window for filtering integration sites.

**reference genome** [-r, --ref_genome] Reference genome. i.e. hg38.

**intsite database** [-i, --intsite_database] Group name for INSPIIRED integration site database.

**specimen database** [-s, --specimen_database] Group name for a specimen management database with specimen information.

**cores** [-c, --cores] Number of cores to parallel processes on. A value of 0 will not utilize the parallel package in R, while values greater than 0 will split appropriate processes to that number of maximum cores (up to the number of cores on the machine).

**output file** [-o, --output] By specifying an output file name, the analysis and summary data will be saved as a *.rds file, which can be quickly read into R using base::readRDS() and assigned to a variable. ie. -o GTSP1322_check -> GTSP1322_check.rds.

## Additional configurations
Please refer to the config.yml file for additional configurations and documentation regarding the options.
