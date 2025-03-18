# Phase/Haplotag Tutorial

Research Questions:
For the gene x sample combinations that could not be fully phased, 

1. Were there any overlapping haploblocks? If so, how many?
2. Were there unphased heterozygous genotypes? If so, where in the gene are they located?
3. How does the IGV look at these sites? Are there structural variants or coverage depth anomalies? Do the IGVs for unphased samples look different from phased samples?

Here is the latest heap map of the phasing results.
![hiphase_heat_map](https://github.com/user-attachments/assets/719edbc6-cc93-4b05-9a0f-a7674159a9e0)

The heat map shows that we could not phase across the entire gene for the following gene x sample combinations:
```
DRB1: HG002, HG003, HG005, HG01106, HG01928, HG02055, HG02630, HG03579, IHW09049, IHWW09071, IHW09118, IHW09175, IHW09245, IHW09251, IHW09359, IHW09364, NA20129, NA24694, NA24695
DRB5: HG01258, HG02055
DQB1: IHW09071, IHW09118
DQB2: IHW09118
DPB1: IHW09224, NA19240
```
Let's look at HG002 HLA-DRB1. 

Here is the IGV with reads colored by tag "PS" (phase set).
![igv_snapshot](https://github.com/user-attachments/assets/72e271c9-7d7c-4c99-9893-7c8077b572b9)

It looks like there are two haploblocks (pink, blue). Let's zoom in on the breakpoint between the two haploblocks. 
![igv_snapshot2](https://github.com/user-attachments/assets/5847b11b-104e-4f16-84c6-d62e64d92eb8)

All of the reads are soft-clipped at this breakpoint. And there are 6 bases without any coverage. Why?

Let's record the problematic coordinates (chr6:32,584,977). This region is in HLA-DRB1 intron 1. 

Look at VCF records in the HLA-A gene (chr6  29941259  29949572)
```
bcftools view -r chr6:29941259-29949572 HG002.vcf.gz | less
```

Example of phased genotype, indicated by "1|0" and PS tag (29921307)
```
chr6    29941293        .       T       C       59.2    PASS    .       GT:GQ:DP:AD:VAF:PL:PS   1|0:56:57:30,27:0.473684:59,0,59:29921307
```

Example of unphased genotype, indicated by "/" and no PS tag.
Homozygous genotypes (i.e, non-ref/non-ref) are not phased by WhatsHap or HiPhase.
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

![igv_snapshot](https://github.com/user-attachments/assets/13cd7d05-3efa-4f87-aac5-5014f10969c1)

![igv_snapshot](https://github.com/user-attachments/assets/9bda6f23-aa94-4411-a85f-6319956f72c8)


