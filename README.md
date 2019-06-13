# ASVCLR
Accurate Structural Variation Caller for Long Reads

-------------------
ASVCLR is an accurate Structural Variation Caller for Long Reads, such as PacBio sequencing and Oxford Nanopore sequencing. ASVCLR can detect both short indels (e.g. <50 bp) and long structural varitions (e.g. >=50 bp), such as duplications, inversions and translocations, and producing fewer false positives with high precise variant margins.  

## Prerequisites
ASVCLR depend on the following libraries and tools:
* HTSlib (http://www.htslib.org/download/)
* Canu (v1.7.1) (https://github.com/marbl/canu/releases)
* BLAT (http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/)
* g++.
The above library and tools should be installed before compiling ASVCLR. Canu and BLAT should be globally accessible in the computer system, these executable files `canu` and 'blat' should be placed or linked into the `$PATH` directory.
Note that: Canu v1.7.1 is about 5 to 10 folds faster than Canu v1.8, but it may cannot construct the assembly results in some genomic regions, therefore if you care more about the accuracy of the results than the running time, please use Canu v1.8 instead.
According to our human chromosome 1 simulated results, the Recall was increased from xxx to xxx, Precision was increased from xxx to xxx, and F1 score was increased from xxx to xxx when switching Canu v1.7.1 to v.18. It is a very small improvement.



## Compiling ASVCLR

The binary file can be generated by typing:
```sh
$ ./autogen.sh
```
and the binary file 'asvclr' will be output into the folder 'bin' in the package directory.


## Quick usage

ASVCLR can be run by typing the `all` command simply:
```sh
$ asvclr all -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam
```
Then, the following commands `detect`, `assemble` and `call` will be performed in turn. The help information can be shown:
```sh
$ asvclr
Program: asvclr (Accurate Structural Variation Caller for Long Reads)
Version: 0.1.0 (using htslib 1.9)

Usage:  asvclr  <command> [options]

Commands:
     detect       detect indel signatures in aligned reads
     assemble     assemble candidate regions and align assemblies
                  back to reference
     call         call indels by alignments of local genome assemblies
     all          run the above commands in turn
```


## Usage

There are three steps to run ASVCLR: detect, assemble and call, which are described as following:  
* __`detect`__: Detect structural variation regions according to abnormal features, including insertions, deletions, duplications, inversions and translocations. Note that these detected regions usually are not accurate enough as well as other SV callers, and may need further processing.  
* __`assemble`__: Perform local assembly for the detected variant regions using Canu, and extract the corresponding local reference.  
* __`call`__: Align the assembly result (contigs) to its local reference using BLAT to generate the sim4 formated alignments, and call each type variations using the BLAT alignments.  

The reference and an sorted BAM file will be the input of ASVCLR, and the variants stored in the VCF file format and translocations in BEDPE file format will be generated as the output.

Usage:
* `detect` step: 
```sh
$ asvclr detect -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam
```
* `assemble` step:
```sh
$ asvclr assemble -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam
```
* `call` step:
```sh
$ asvclr call -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam
```

### Detect Step


### Assemble Step


### Call Step


## Output Result Description



------------------
If you have problems or some suggestions, please contact: xzhu@hrbnu.edu.cn without hesitation. 

---- Enjoying!!! -----

