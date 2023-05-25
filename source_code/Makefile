
SCRIPTS = ./
OUTPUT_DIR = ../outputs
ifeq ($(PREFIX),)
override PREFIX = output
endif
GENOME_LEN = $(shell cat $(GENOME_DIR)/chrLength.txt)

single_end: $(OUTPUT_DIR)/$(PREFIX)-report.html $(OUTPUT_DIR)/$(PREFIX)-insertion_visualization.txt $(OUTPUT_DIR)/$(PREFIX)-collapsed.sam
	@touch single_end
	@echo "\nDONE\n"

# run R script for insertion visualization
$(OUTPUT_DIR)/$(PREFIX)-insertion_visualization.txt: $(OUTPUT_DIR)/$(PREFIX)-insertions.txt
	Rscript $(SCRIPTS)/visualize_insertions.R $(OUTPUT_DIR)/$(PREFIX) $(GENOME_LEN)

# run filtering and collapsing reads
$(OUTPUT_DIR)/$(PREFIX)-collapsed.sam: $(OUTPUT_DIR)/$(PREFIX).Aligned.out.sam
	$(SCRIPTS)/filter_reads_match+mismatch.py --input $(OUTPUT_DIR)/$(PREFIX).Aligned.out.sam --output $(OUTPUT_DIR)/$(PREFIX)-filtered.sam
	$(SCRIPTS)/collapse_reads.py --input $(OUTPUT_DIR)/$(PREFIX)-filtered.sam --output $(OUTPUT_DIR)/$(PREFIX)-collapsed.sam
	igvtools sort $(OUTPUT_DIR)/$(PREFIX)-collapsed.sam $(OUTPUT_DIR)/$(PREFIX)-collapsed.sorted.sam
	igvtools index $(OUTPUT_DIR)/$(PREFIX)-collapsed.sorted.sam

# run R script for generating plots
$(OUTPUT_DIR)/$(PREFIX)-report.html: $(OUTPUT_DIR)/$(PREFIX)-coverage.txt $(OUTPUT_DIR)/$(PREFIX)-insertions.txt $(OUTPUT_DIR)/$(PREFIX)-deletions.txt
	Rscript $(SCRIPTS)/plot.R $(OUTPUT_DIR)/$(PREFIX) $(GENOME_LEN)

# run python script creating text files with insertions, deletions and coverage
$(OUTPUT_DIR)/$(PREFIX)-coverage.txt $(OUTPUT_DIR)/$(PREFIX)-insertions.txt $(OUTPUT_DIR)/$(PREFIX)-deletions.txt: $(OUTPUT_DIR)/$(PREFIX).Aligned.out.sam
	$(SCRIPTS)/count_reference_coverage.py --input $(OUTPUT_DIR)/$(PREFIX).Aligned.out.sam --outputCov $(OUTPUT_DIR)/$(PREFIX)-coverage.txt --outputIns $(OUTPUT_DIR)/$(PREFIX)-insertions.txt --outputDel $(OUTPUT_DIR)/$(PREFIX)-deletions.txt

# run STAR aligner
$(OUTPUT_DIR)/$(PREFIX).Aligned.out.sam: $(OUTPUT_DIR)/$(PREFIX)-_trimmed.fq	
ifeq ("$(wildcard $(GENOME_DIR)/Genome)","")
	@echo "\n\n ! Generated genome indexes in '$(GENOME_DIR)' doesn't exist \n\n"
	@exit 1
else
	STAR --genomeDir $(GENOME_DIR) --readFilesIn $(OUTPUT_DIR)/$(PREFIX)-_trimmed.fq --outFileNamePrefix $(OUTPUT_DIR)/$(PREFIX). --outSAMattributes NH HI AS nM MD --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --outFilterScoreMinOverLread 0.5
endif

# run Trim Galore
$(OUTPUT_DIR)/$(PREFIX)-_trimmed.fq:
ifeq ("$(wildcard $(SAMPLE))","")
	@echo "\n\n ! Input sample file '$(SAMPLE)' doesn't exist \n\n"
	@exit 1
else
	trim_galore -o $(OUTPUT_DIR) --no_report_file --basename $(PREFIX)- $(SAMPLE)
endif



paired_end: $(OUTPUT_DIR)/$(PREFIX)_paired-report.html $(OUTPUT_DIR)/$(PREFIX)_paired-insertion_visualization.txt $(OUTPUT_DIR)/$(PREFIX)_paired-collapsed.sam
	@touch paired_end
	@echo "\nDONE\n"

# run R script for insertion visualization
$(OUTPUT_DIR)/$(PREFIX)_paired-insertion_visualization.txt: $(OUTPUT_DIR)/$(PREFIX)_paired-insertions.txt
	Rscript $(SCRIPTS)/visualize_insertions.R $(OUTPUT_DIR)/$(PREFIX)_paired $(GENOME_LEN)

# run filtering and collapsing reads
$(OUTPUT_DIR)/$(PREFIX)_paired-collapsed.sam: $(OUTPUT_DIR)/$(PREFIX)_paired.Aligned.out.sam
	$(SCRIPTS)/filter_reads_match+mismatch.py --input $(OUTPUT_DIR)/$(PREFIX)_paired.Aligned.out.sam --output $(OUTPUT_DIR)/$(PREFIX)_paired-filtered.sam
	$(SCRIPTS)/collapse_reads.py --input $(OUTPUT_DIR)/$(PREFIX)_paired-filtered.sam --output $(OUTPUT_DIR)/$(PREFIX)_paired-collapsed.sam
	igvtools sort $(OUTPUT_DIR)/$(PREFIX)_paired-collapsed.sam $(OUTPUT_DIR)/$(PREFIX)_paired-collapsed.sorted.sam
	igvtools index $(OUTPUT_DIR)/$(PREFIX)_paired-collapsed.sorted.sam

# run R script for generating plots
$(OUTPUT_DIR)/$(PREFIX)_paired-report.html: $(OUTPUT_DIR)/$(PREFIX)_paired-coverage.txt $(OUTPUT_DIR)/$(PREFIX)_paired-insertions.txt $(OUTPUT_DIR)/$(PREFIX)_paired-deletions.txt
	Rscript $(SCRIPTS)/plot.R $(OUTPUT_DIR)/$(PREFIX)_paired $(GENOME_LEN)

# run python script creating text files with insertions, deletions and coverage
$(OUTPUT_DIR)/$(PREFIX)_paired-coverage.txt $(OUTPUT_DIR)/$(PREFIX)_paired-insertions.txt $(OUTPUT_DIR)/$(PREFIX)_paired-deletions.txt: $(OUTPUT_DIR)/$(PREFIX)_paired.Aligned.out.sam
	$(SCRIPTS)/count_reference_coverage.py --input $(OUTPUT_DIR)/$(PREFIX)_paired.Aligned.out.sam --outputCov $(OUTPUT_DIR)/$(PREFIX)_paired-coverage.txt --outputIns $(OUTPUT_DIR)/$(PREFIX)_paired-insertions.txt --outputDel $(OUTPUT_DIR)/$(PREFIX)_paired-deletions.txt

# run STAR aligner
$(OUTPUT_DIR)/$(PREFIX)_paired.Aligned.out.sam: $(OUTPUT_DIR)/$(PREFIX)-merged.assembled.fastq
ifeq ("$(wildcard $(GENOME_DIR)/Genome)","")
	@echo "\n\n ! Generated genome indexes in '$(GENOME_DIR)' doesn't exist \n\n"
	@exit 1
else
	STAR --genomeDir $(GENOME_DIR) --readFilesIn $(OUTPUT_DIR)/$(PREFIX)-merged.assembled.fastq --outFileNamePrefix $(OUTPUT_DIR)/$(PREFIX)_paired. --outSAMattributes NH HI AS nM MD --scoreGapNoncan 0 --scoreGapGCAG 0 --scoreGapATAC 0 --outFilterScoreMinOverLread 0.5
endif

# run Pear
$(OUTPUT_DIR)/$(PREFIX)-merged.assembled.fastq: $(OUTPUT_DIR)/$(PREFIX)-_val_1.fq $(OUTPUT_DIR)/$(PREFIX)-_val_2.fq
	pear -f $(OUTPUT_DIR)/$(PREFIX)-_val_1.fq -r $(OUTPUT_DIR)/$(PREFIX)-_val_2.fq -o $(OUTPUT_DIR)/$(PREFIX)-merged

# run Trim Galore
$(OUTPUT_DIR)/$(PREFIX)-_val_1.fq $(OUTPUT_DIR)/$(PREFIX)-_val_2.fq:
ifeq ("$(wildcard $(FORWARD))","")
	@echo "\n\n ! Input sample file '$(FORWARD)' doesn't exist \n\n"
	@exit 1
endif
ifeq ("$(wildcard $(REVERSE))","")
	@echo "\n\n ! Input sample file '$(REVERSE)' doesn't exist \n\n"
	@exit 1
endif
	trim_galore -o $(OUTPUT_DIR) --no_report_file --basename $(PREFIX)- --paired $(FORWARD) $(REVERSE)

