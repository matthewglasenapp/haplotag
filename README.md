Tutorial in haplotag.md

The data is located on Hummingbird in ```/hb/groups/cornejo_lab/matt/sample_bams_vcfs/```

Each sample has its own directory containing multiple files. Here is an example of the files for HG002.

```
HG002.haplotag.bam
HG002.haplotag.bam.bai
HG002.haplotag.hap1.bam
HG002.haplotag.hap1.bam.bai
HG002.haplotag.hap2.bam
HG002.haplotag.hap2.bam.bai
HG002.haplotag.no_hp.bam
HG002.haplotag.no_hp.bam.bai
HG002.vcf.gz
HG002.vcf.gz.tbi
```

```HG002.haplotag.bam``` is the alignment file containing both haplotypes. Begin your investigation with this file. I have also split the BAM by haplotype into ```HG002.haplotag.hap1.bam``` and ```HG002.haplotag.hap2.bam```. Reads that could not be assigned a haplotype are in ```HG002.haplotag.no_hp.bam```. ```HG002.vcf.gz``` is the VCF file containing the genotypes called by DeepVariant. 

It is about 2.6 GB of data.

Download it to your local device with the following command. Substitute <username> with your username and specify custom path instead of ~/Desktop/

```
rsync -avz --progress --chmod=ugo+rwx <username>@hb.ucsc.edu:/hb/groups/cornejo_lab/matt/sample_bams_vcfs ~/Desktop/
```
