This procedure describes what I do to analyze the results of screens
performed with lentiCRISPR v2 (https://www.addgene.org/52961/).

The genomic DNA is amplified with 4 staggered forward primers. This is to
avoid having the same sequence in every read (the constant sequence right
in front of the gRNA). Also, the primers contain a random stretch of BVHD
nucleotides, the purpose of which is to create some diversity at the
beginning of the read (they carry no information).

The constant part in front of the gRNA is `ACTCTTGTGGAAAGGACGAAACACCG`.
The sequence after the gRNA in the plasmid is `TTTTAGAGCTAGAAAT`.
The P5 sequence is `AATGATACGGCGACCACCGAGATCT`.
The P7 sequence is `CAAGCAGAAGACGGCATACGAGAT`.
The forward Illumina sequencing primer is `ACACTCTTTCCCTACACGACGCTCTTCCGATCT`.
The reverse Illmunina sequencing primer is `GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT`.

The sequences of the 4 forward primers are
1. `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTBVHDBVHDBVACTCTTGTGGAAAGGACGAAACACCG`
2. `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTVHDBVHDBVHCACTTGTGGAAAGGACGAAACACCG`
3. `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTDBVHDBVHDBTTGTGGAAAGGACGAAACACCG`
4. `AATGATACGGCGACCACCGAGATCTACACTCTTTCCCTACACGACGCTCTTCCGATCTHDBVHDTTGTGGAAAGGACGAAACACCG`

The reverse primer is
`CAAGCAGAAGACGGCATACGAGATNNNNNNNNGTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTCCAATTCCCACTCCTTTCAAGACCT`

#### Extract gRNAs

`python extract_gRNA.py file.fatsq | starcode -d2 > file.stc`

#### Map reads in the genome
`bwa mem -T15 -L5,0 index <(python bcd_to_fasta.py file.stc) | samtools view -Sq6 - > file.sam`

#### Get coordinates of exons
`curl  -s "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/knownGene.txt.gz" | gunzip -c | awk '{n=int($8); split($9,S,/,/);split($10,E,/,/); for(i=1;i<=n;++i) {printf("%s,%s,%s,%s,%s\n",$1,$2,$3,S[i],E[i]);} }' > exons.txt`

#### Map gRNAs to exons
`python map_gRNA_to_exons.py exons.txt file.sam > file.scores`

#### Merge results into a single file
`python merge_scores.py *.scores > scores.txt`

#### Example output

Below is an example of output.

    sgRNA                   gene     file1   file2
    CAGATTCCTCATGTTTGTGC    SNAI2    237        12
    GCAGACTTTGCGGCTATGCG    KIF20A   89        231
    ACCTGTCAAGCAGTATCCCA    KHDRBS1  86          0
    CATGTGCTGGTCCTGTAGGC    OR52B4   69          0


The first column is the sequence of the sgRNA. The second column is the
name of the gene where it maps. The following columns are the number of
reads in each file listed in input.
