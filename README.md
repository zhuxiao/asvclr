# ASVCLR
Accurate Structural Variation Caller for Long Reads

------------------------------------------------------------------------------
ASVCLR is an accurate Structural Variation Caller for Long Reads, such as PacBio sequencing and Oxford Nanopore sequencing. ASVCLR can detect indels (insertions and deletions) both for short indels (e.g. <50 bp) and long structural varitions (e.g. >=50 bp), and producing precise variant margins with high presision and fewer false positives.  

There are three steps to run ASVCLR: detect, assemble and call, which are described as following:  
>    (1) detect -- Detect structural variation regions according to abnormal features, including insertions, deletions, duplications, inversions and translocations. Note that these detected regions usually are not accurate enough as well as other SV callers, and may need further processing.  
>    (2) assemble -- Perform local assembly for the detected variant regions using Canu, and extract the corresponding local reference.  
>    (3) call -- Align the assembly result (contigs) to its local reference using BLAT to generate the sim4 formated alignments, and call each type variations using the BLAT alignments.  


The reference and an sorted BAM file can be the input of ASVCLR, and the variants stored in the VCF file format is generated as the output.  

Usage: 
>    (1) detect step:  
>        asvlr detect -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam  
>    (2) assemble step:  
>        asvlr assemble -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam  
>    (3) call step:  
>        asvlr call -t 14 -c 20000 -f hg38.fa hg38_ngmlr_sorted.bam  

------------------------------------------------------------------------------
If you have problems or some suggestions, please connect: xzhu@hrbnu.edu.cn  

---- Enjoying!!! -----

