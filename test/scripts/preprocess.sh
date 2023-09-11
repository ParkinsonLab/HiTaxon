TESTDIR=$1 

#For each simulated instance, process the aggregated reads
for i in {1..10}; do
    #Convert from FASTA to FASTQ; needed for fastq-join
    seqtk seq -F '#' "${TESTDIR}/simulated_reads/$i/merged1.fa" > "${TESTDIR}/simulated_reads/$i/merged1.fq"
    seqtk seq -F '#' "${TESTDIR}/simulated_reads/$i/merged2.fa" > "${TESTDIR}/simulated_reads/$i/merged2.fq"
    #Merge paired-end reads
    fastq-join "${TESTDIR}/simulated_reads/$i/merged1.fq" "${TESTDIR}/simulated_reads/$i/merged2.fq" -o "${TESTDIR}/simulated_reads/$i/merged.fq"
    #Convert merged FASTQ to FASTA
    seqtk seq -a "${TESTDIR}/simulated_reads/$i/merged.fqjoin" > "${TESTDIR}/simulated_reads/$i/merged.fa"
done