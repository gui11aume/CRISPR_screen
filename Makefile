# genecode_v24.txt.gz downloaded from https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=685033455_QRKe8pfLxa6ILnybGSr1VGIW1Equ&clade=mammal&org=Human&db=hg38&hgta_group=genes&hgta_track=knownGene&hgta_table=0&hgta_regionType=genome&position=chr1%3A11102837-11267747&hgta_outputType=primaryTable&hgta_outFileName=geneid.txt.gz
# gencode.v24.annotation.gtf.gz downloaded from https://www.gencodegenes.org/releases/24.html
SHELL=bash

TARGETS= \
	K015.scores \
	K016.scores \
	K017.scores \
	K018.scores \
	K019.scores \
	K020.scores \
	K021.scores \
	K022.scores \
	K023.scores \
	K024.scores \
	K025.scores \
	K026.scores \
	K027.scores \
	K028.scores \
	K030.scores \
	K031.scores

PARAMS= -T15 -L5,0
INDEX= /mnt/shared/seq/bwa/GRCh38/hg38.fasta

all: $(TARGETS)

%.stc: %_*.fastq.gz
	python extract_gRNA.py $< | starcode -d2 > $@

%.sam: %.stc
	python bcd_to_fasta.py $< > $*.fasta
	bwa mem $(PARAMS) $(INDEX) $*.fasta | samtools view -Sq6 - > $@
	rm $*.fasta

%.scores: %.sam exons.txt
	python map_gRNA_to_exons.py gencode.v24.annotation.gtf.gz exons.txt $< > $@

exons.txt:
	zcat genecode_v24.txt.gz | awk '{n=int($$8); split($$9,S,/,/);split($$10,E,/,/); for(i=1;i<=n;++i) {printf("%s,%s,%s,%s,%s\n",$$NF,$$2,$$3,S[i],E[i]);} }' > exons.txt
