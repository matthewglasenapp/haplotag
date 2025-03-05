# Phase/Haplotag Tutorial

Look at VCF records in the HLA-A gene (chr6  29941259  29949572)
```
bcftools view -r chr6:29941259-29949572 HG002.vcf.gz | less
```

Example of phased genotype, indicated by "|" and PS tag (29921307)
```
chr6    29941293        .       T       C       59.2    PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:56:57:30,27:0.473684:59,0,59:29921307
```

Example of unphased genotype, indicated by "/" and no PS tag
Homozygous genotypes (i.e, non-ref/non-ref) are not phased by WhatsHap nor HiPhase
```
chr6    29941404        .       A       C       59.8    PASS    .       GT:GQ:DP:AD:VAF:PL      1/1:55:62:0,62:1:59,56,0
```

FORMAT Fields definitions:
  GT: Genotype
  GQ: Conditional genotype quality
  DP: Read depth
  AD: Read depth for each allele
  PL: Phred-scaled genotype likelihoods rounded to the closest integer
  PS: Phase set identifier

Look at BAM records in the HLA-A gene (chr6  29941259  29949572)
```
samtools view HG002.haplotag.bam chr6:29941259-29949572 | less
```

Example FORMAT fields from a haplotagged read
```
RG:Z:m84039_240622_113450_s1	qs:i:0	qe:i:7293	mg:f:99.3829	NM:i:193	HP:i:2	PC:i:3720	PS:i:29921307
```
This read belongs to haplotype 2 in phase set 29921307

Split BAM file by HP tag (I've done this already for each sample)
```
samtools view -h HG002.haplotag.bam | grep -E 'HP:i:1|^@' | samtools view -b -o HG002.haplotag.hap1.bam
samtools view -h HG002.haplotag.bam | grep -E 'HP:i:2|^@' | samtools view -b -o HG002.haplotag.hap2.bam
samtools view -h HG002.bam | grep -E -v 'HP:i:[12]' | samtools view -b -o HG002.haplotag.no_hp.bam
```

Index each new BAM file (I've done this already for each sample)
```
samtools index HG002.haplotag.hap1.bam
samtools index HG002.haplotag.hap2.bam
samtools index G002.haplotag.no_hp.bam
```

Load each file separately in IGV!

In IGV, right click on reads, click "Color alignments by" -> "tag" and type "HP"

Here is an example of great phasing in HLA-A. 

![hla-a](https://github.com/user-attachments/assets/03cfb635-60cc-401e-9a4b-2f890e1c9ff1)

You can also load the BAM from before splitting by HP tag. You can color alignments by PS or HP tag. "PS" stands for "Phase Set" and represents the haploblocks defined by HiPhase. 




