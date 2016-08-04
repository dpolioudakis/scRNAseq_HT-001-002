fingerprint_cmds.txt
	while read sampleID; do echo python fingerprint.py ../data/fastq/Merged_For_Variant_Calling/${sampleID}_R1.fastq ../data/fastq/Merged_For_Variant_Calling/${sampleID}_R2.fastq ../data/fastq/Merged_For_Variant_Calling/${sampleID} >> fingerprint_cmds.txt; done < ../analysis/tables/Human_Only_Capture_Sites_10^5Hs_10^5Mm.txt
	Made for inputing samples into fingerprint.py for variant calling
	Run fingerprint.py on hoffman2
