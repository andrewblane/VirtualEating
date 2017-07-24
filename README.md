# VirtualEating
Andrew Lane, University of California, Berkeley

## Changelog:
20150930 (v. 2.0): Major overhaul. The iPython notebook examples/implementation have been
deprecated (moved to the archive folder) and the core logic abstracted into the
`eating` package.
20150530 (v. 1.0): Initial release.

## Description
The `eating` package is a set of functions designed to help in planning, implementing and
analyzing CRISPR-EATING enzymatically generated sgRNA libraries for genome labeling
or screening as decribed in [Lane et al., 2015](http://www.cell.com/developmental-cell/abstract/S1534-5807(15)00392-5)

Given that CRISPR-EATING protocol produces sgRNAs from any source of DNA, it's essential to be 
able to predict, tune and analyze the sgRNA output of the protocol. 

**Please contact andylane@gmail.com if you would like access to a multithreaded Amazon Web Services-deployable version. This version has modfications that separate database and parallelized scoring/BLAST threads, making it much faster and more scaleable when scoring large numbers of guides but somewhat more complex to deploy. SQLite databases with scores for all possible guides in hg19/GRCh37 are also available.*
        
For each function, parameters are described in the function's docstring.

The `eating` package addresses this as follows:
- Prediction is covered by `eating.base`:
    - Given an input DNA sequence, `base` provides functions to predict the output of CRISPR-EATING
when the molecular biology protocol is applied to that sequence. 
    - `base.digest_target` simulates PAM selection, cutting and MmeI-20mer trimming of input DNA to produce a list of sgRNAs
    - `base.score_guides` uses the BioPython BLAST interface on a local BLAST installation and database to generate target-genome-wide off-target analysis of sgRNAs and processes the resulting hit data into a per-guide score by implementing the sgRNA off-target scoring algorithm described in [Hsu et al., 2013](https://www.nature.com/nbt/journal/v31/n9/full/nbt.2647.html). 
        - The implementation memoizes and stores the calculated scores in an sqlite database that can be re-used in other workflows
    - `base.count_non_overlapping_guides` is a helper function to aggregate a linear cluster of guides and count the true non-overlapping count of guides within that region. This is useful for applications such as CRISPR-imaging, where guide tiling density is the desired characteristic.
    
- Tuning is covered by `eating.designpcr`
    - `designpcr.find_amplicons` searches for regions of an input genome that, when PCR amplified, will be digested into 
       tens or hundreds of highly-specific sgRNAs. The resulting `Amplicon` class stores information about the location of the guides
       within the genome, the identities of the individual guides, their scores and, importantly, how close these clusters are to 
       low-scoring guides that must be avoided in downstream PCR.
    - `designpcr.primer_search` proposes and specificity-screens primers to amplify chosen amplicons
        - Primers are proposed using the fast [primer3-py](https://github.com/libnano/primer3-py) interface
        - `designpcr.screen_primer_in_silico_pcr` implements an in-silico PCR strategy to determine if primers proposed by primer3-py
        are likely to be successful on the target genome. Each primer is  BLAST-ed against the template genome and the distance and orientation between matching sites in the genome is analyzed to ensure that off-target binding of primers cannot produce an amplicon. In practice, around *85%* of the predicted-good primers output by this tool produce a single band in PCR using the protocol described the CRISPR-EATING paper.
    - `designpcr.collect_good_primers` outputs primers and amplicon attributes, including size of amplicon and number of guides generated from an individual amplicon into CSV file for ordering of oligonucleotides from IDT.
    - Sample Output:
    ```
    Gene	Sequence_id	forward_seq	forward_start	forward_length	forward_tm	forward_gc	reverse_seq	reverse_start	reverse_length	reverse_tm	reverse_gc	input_seq_length	PCR_product_length	Guides_Contained	Non-overlapping Guide Count
    GNB2L1	6282	GGGATAGGGACGGGGAGAAC	432	20	61.12018408	65	ACCCTCCGGAAGCACAGTT	1600	19	61.1473871	57.89473684	1596	1168	41	36
    EF1A	3947	ACAGAAGCAACCAAAAATCAAACTT	157	25	58.82737049	32	TCCCTTCCAGGCGGCCTC	1948	18	63.80873271	72.22222222	1942	1791	43	36
    EEF1G	13188	AATGCCACTCTCCAGGATGA	202	20	58.41176058	50	AGGAGGTGGGAGGGACAG	1984	18	59.54788591	66.66666667	1989	1782	59	52
    ```
    
- Analysis is covered by a limited number of functions (`eating.visualize`). These are designed to operate on FASTQ files resulting from HiSeq 2000 sequencing of CRISPR-EATING libraries.
    - `visualize.plot_complexity` attempts to describe the diversity of the resulting sgRNA library by plotting unique sequences against
    sequence count. This is intended to be a quick screen for gross PCR "jackpotting" (aka the overamplification of one or a few sequences in PCR).
    - `visualize.string2feat` annotates Bio.Seq objects with the locations of sgRNAs observed in sequencing. This is useful for troubleshooting, e.g., failures of the Mung-Bean nuclease step where observed guides are preferentially found at the termini of amplicons. 


### Installation:
Copy the `eating` folder from this git repo into to your python
site-packages directory or your project working directory.

## Usage
```python
import eating
```
Consult individual function docstrings for instructions or contact andylane@gmail.com for code orientation.
