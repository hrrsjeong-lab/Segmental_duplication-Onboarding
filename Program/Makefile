# Variables
SHELL = /bin/bash -O globstar
CPUS = 10
MEMS = 100
MEM_PER_CPU = $(shell echo '$(MEMS)/$(CPUS)' | bc)

# Programs
SLURM = srun --chdir=$(abspath .) --cpus-per-task=$(CPUS) --export=ALL --job-name="SDs_$(@F)" --mem=$(MEMS)G
PYTHON = /home/Live/jaewoong/JeongLab_Onboarding/3_Segmental_duplications/.venv/bin/python3 -B
BISER = /home/Live/jaewoong/JeongLab_Onboarding/3_Segmental_duplications/.venv/bin/biser
SEDEF = /home/Share/Tools/sedef-1.1/sedef.sh
BGZIP = /usr/bin/bgzip
BEDTOOLS = /usr/bin/bedtools
SAMTOOLS = /usr/bin/samtools
rwildcard = $(foreach d,$(wildcard $(1:=/*)),$(call rwildcard,$d,$2) $(filter $(subst *,%,$2),$d))
product = $(foreach x,$1,$(addprefix $x.,$2))

# References
HUMAN_FASTA = /home/Share/References/chm13/chm13v2.0.maskedY.fasta
HUMAN_GENE_GTF = /home/Share/References/chm13/hs1.ncbiRefSeq.gtf.gz
HUMAN_LEN_TXT = /home/Share/References/chm13/hs1.chrom.sizes.txt
HUMAN_SD_BED = /home/Share/References/chm13/chm13v2.0_SD.full.bed

# Variables
SPECIES = Human# Chimp Bonob
TYPES = dip mat pat

# Options
LOGGING = 1>$(abspath $@).log 2>&1
LOGGING_DIR = 1>$(abspath $(@D)).log 2>&1

all: step01 step02 step03 step04 step05
.PHONY: all

clean:
.PHONY: clean

test:
.PHONY: test

# Step 01 (Softmask)
01_Softmask:
	mkdir -p $@

01_Softmask/Human.dip.fasta: /home/Store/HG002v1.1/hg002v1.1.fasta.gz /home/Store/HG002v1.1/hg002v1.1.dip.rmsk.3col.bed | 01_Softmask
	$(SLURM) $(BGZIP) --decompress --stdout --threads $(CPUS) $(filter %.fasta.gz,$^) > $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))
	$(SLURM) $(BEDTOOLS) maskfasta -fi $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^)) -bed $(filter %.bed,$^) -fo $(abspath $@) -soft $(LOGGING)
	$(SLURM) $(SAMTOOLS) faidx $@
	rm -fv $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))

01_Softmask/Human.mat.fasta: /home/Store/HG002v1.1/hg002v1.1.mat.fasta.gz /home/Store/HG002v1.1/hg002v1.1.mat.rmsk.3col.bed | 01_Softmask
	$(SLURM) $(BGZIP) --decompress --stdout --threads $(CPUS) $(filter %.fasta.gz,$^) > $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))
	$(SLURM) $(BEDTOOLS) maskfasta -fi $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^)) -bed $(filter %.bed,$^) -fo $(abspath $@) -soft $(LOGGING)
	$(SLURM) $(SAMTOOLS) faidx $@
	rm -fv $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))

01_Softmask/Human.pat.fasta: /home/Store/HG002v1.1/hg002v1.1.pat.fasta.gz /home/Store/HG002v1.1/hg002v1.1.pat.rmsk.3col.bed | 01_Softmask
	$(SLURM) $(BGZIP) --decompress --stdout --threads $(CPUS) $(filter %.fasta.gz,$^) > $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))
	$(SLURM) $(BEDTOOLS) maskfasta -fi $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^)) -bed $(filter %.bed,$^) -fo $(abspath $@) -soft $(LOGGING)
	$(SLURM) $(SAMTOOLS) faidx $@
	rm -fv $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))

01_Softmask/Chimp.%.fasta: /home/Store/T2T-Primates/Chimpanzee/mPanTro3.%.cur.20231122.fasta.gz /home/Store/T2T-Primates/Chimpanzee/mPanTro3_v2.0.ideogram.bed | 01_Softmask
	$(SLURM) $(BGZIP) --decompress --stdout --threads $(CPUS) $(filter %.fasta.gz,$^) > $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))
	$(SLURM) $(BEDTOOLS) maskfasta -fi $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^)) -bed $(filter %.bed,$^) -fo $(abspath $@) -soft $(LOGGING)
	$(SLURM) $(SAMTOOLS) faidx $@
	rm -fv $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))

01_Softmask/Chimp.mat.fasta: 01_Softmask/Chimp.hap1.fasta
	ln -v $(abspath $^ $@) $(LOGGING)
	$(SLURM) $(SAMTOOLS) faidx $@ $(LOGGING)

01_Softmask/Chimp.pat.fasta: 01_Softmask/Chimp.hap2.fasta
	ln -v $(abspath $^ $@)
	$(SLURM) $(SAMTOOLS) faidx $@

01_Softmask/Bonob.%.fasta: /home/Store/T2T-Primates/Bonobo/mPanPan1.%.cur.20231122.fasta.gz /home/Store/T2T-Primates/Bonobo/mPanPan1_v2.0.ideogram.bed | 01_Softmask
	$(SLURM) $(BGZIP) --decompress --stdout --threads $(CPUS) $(filter %.fasta.gz,$^) > $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))
	$(SLURM) $(BEDTOOLS) maskfasta -fi $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^)) -bed $(filter %.bed,$^) -fo $(abspath $@) -soft $(LOGGING)
	$(SLURM) $(SAMTOOLS) faidx $@
	rm -fv $(subst .fasta.gz,.fasta,$(filter %.fasta.gz,$^))

01_Softmask/check.txt: Python/count_soft.py $(HUMAN_FASTA) $(addprefix 01_Softmask/,$(addsuffix .fasta,$(call product,$(SPECIES),$(TYPES))))
	$(PYTHON) $^ $@

step01: $(addprefix 01_Softmask/,$(addsuffix .fasta,$(call product,$(SPECIES),$(TYPES))))
.PHONY: step01

# Step 02 (BISER)
02_BISER:
	mkdir -p $@

02_BISER/chm13.bedpe: $(HUMAN_FASTA) | 02_BISER
	$(SLURM) $(BISER) --gc-heap $(MEM_PER_CPU)G --threads $(CPUS) --keep-contigs --max-chromosome-size 10000000 --temp $(abspath .) --output $(abspath $@) $^ $(LOGGING)

02_BISER/%.bedpe: 01_Softmask/%.fasta | 02_BISER
	$(SLURM) $(BISER) --gc-heap $(MEM_PER_CPU)G --threads $(CPUS) --keep-contigs --max-chromosome-size 10000000 --temp $(abspath .) --output $(abspath $@) $^ $(LOGGING)

step02: 02_BISER/chm13.bedpe $(addprefix 02_BISER/,$(addsuffix .bedpe,$(call product,$(SPECIES),$(TYPES))))
.PHONY: step02

# Step 03 (SEDEF)
03_SEDEF:
	mkdir -p $@

03_SEDEF/chm13/final.bed: $(HUMAN_FASTA) | 03_SEDEF
	$(SLURM) /usr/bin/script --quiet --return --command '$(SEDEF) -o $(abspath $(@D)) -j $(CPUS) -f $^' $(LOGGING_DIR)

03_SEDEF/%/final.bed: 01_Softmask/%.fasta | 03_SEDEF
	$(SLURM) $(SEDEF) -o $(abspath $(@D)) -j $(CPUS) -f $^ $(LOGGING_DIR)

step03: 03_SEDEF/chm13/final.bed $(addprefix 03_SEDEF/,$(addsuffix /final.bed,$(call product,$(SPECIES),$(TYPES))))
.PHONY: step03

# Step 04 (Figure of SEDEF)
04_Figure_sedef:
	mkdir -p $@

04_Figure_sedef/%.zip: CPUS = 1
04_Figure_sedef/%.zip: MEMS = 50
04_Figure_sedef/%.zip: Python/read_sedef.py 03_SEDEF/%/final.bed $(HUMAN_SD_BED) $(HUMAN_LEN_TXT) | 04_Figure_sedef
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

step04: 04_Figure_sedef/chm13.zip
.PHONY: step04

# Step 05 (Figure 0f BISER)
05_Figure_biser:
	mkdir -p $@

05_Figure_biser/%.zip: CPUS = 1
05_Figure_biser/%.zip: MEMS = 50
05_Figure_biser/%.zip: Python/read_biser.py 02_BISER/%.bedpe $(HUMAN_SD_BED) $(HUMAN_LEN_TXT) | 05_Figure_biser
	$(SLURM) $(PYTHON) $(abspath $^ $@) $(LOGGING)

step05: 05_Figure_biser/chm13.zip
.PHONY: step05
