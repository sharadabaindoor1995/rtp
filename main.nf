#!/usr/bin/env nextflow

if(has_extension(params.input, ".csv")){
   
   csv_file = file(params.input, checkIfExists: true)
   ch_input = extract_data(csv_file)

}else{

   ch_input = Channel.fromFilePairs(params.input, checkIfExists: true)

}

( ch_qc_reads, ch_raw_reads) = ch_input.into(2)

ch_fasta = Channel.value(file(params.fasta))
ch_gtf = Channel.value(file(params.gtf))

process FASTQC{
    publishDir "${params.outdir}/quality_control/fastqc", mode: 'copy'

    input:
    tuple val(base), file(reads) from ch_qc_reads

    output:
    file("*.{html,zip}") into ch_multiqc

    script:
    """
    fastqc -q $reads
    """
}

process TX{
    publishDir "${params.outdir}/reference", mode: 'copy'

    input:
    file(fasta) from ch_fasta
    file(gtf) from ch_gtf

    output:
    file("${fasta.baseName}.tx.fa") into transcriptome_created

    script:
    """
    gffread -F -w "${fasta.baseName}.tx.fa" -g $fasta $gtf
    """
}


process INDEX{
    publishDir "${params.outdir}/results/reference", mode: 'copy'

    input:
    file(tx) from transcriptome_created

    output:
    file("*.idx") into index_created

    script:
    """
    kallisto index -i ${tx.simpleName}.idx $tx
    """
}

process KALLISTO_QUANT{
    publishDir "${params.outdir}/results/kallisto",mode:'copy'

    input:
    tuple val(base),file(reads) from ch_raw_reads
    file(index) from index_created

    output:
    tuple val(base), file("${base}") into kallisto_out
    file("${base}.kallisto_log") into kallisto_logs

    script:
    """
    kallisto quant \
    -i $index \
    -t 2 \
    -o ${base}/ \
    --bias \
    $reads &> "${base}.kallisto_log"

    """
}

process MULTIQC{
    publishDir "${params.outdir}/quality_control/multiqc", mode: 'copy'

    input:
    file(htmls) from ch_multiqc.collect()
    file(kallisto_logs) from kallisto_logs.collect()

    output:
    file("*.html") into ch_out

    script:
    """
    multiqc .
    """
}
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "error:  Invalid CSV input - malformed row (e.g. missing column) in ${row}, consult documentation."
    return true
}

def return_file(it) {
    if (!file(it).exists()) exit 1, "error: Cannot find supplied FASTQ input file. Check file: ${it}"
    return file(it)
}

def has_extension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}


def extract_data(csvFile){
    Channel
        .fromPath(csvFile)
        .splitCsv(header: true, sep: ',')
        .map{ row ->

        def expected_keys = ["Sample_ID", "Read1", "Read2"]
        if(!row.keySet().containsAll(expected_keys)) exit 1, "error: Invalid CSV input - malformed column names. Please use the column names 'Sample_ID', 'Read1', 'Read2'."

        checkNumberOfItem(row, 3)

        def samples = row.Sample_ID
        def read1 = row.Read1.matches('NA') ? 'NA' : return_file(row.Read1)
        def read2 = row.Read2.matches('NA') ? 'NA' : return_file(row.Read2)

        if( samples == '' || read1 == '' || read2 == '' ) exit 1, "error: a field does not contain any information. Please check your CSV file"
        if( !has_extension(read1, "fastq.gz") && !has_extension(read1, "fq.gz") && !has_extension(read1, "fastq") && !has_extension(read1, "fq")) exit 1, "error: A R1 file has a non-recognizable FASTQ extension. Check: ${r1}"
        if( !has_extension(read2, "fastq.gz") && !has_extension(read2, "fq.gz") && !has_extension(read2, "fastq") && !has_extension(read2, "fq")) exit 1, "error: A R2 file has a non-recognizable FASTQ extension. Check: ${r2}"

        // output tuple mimicking fromFilePairs
        [ samples, [read1, read2] ]

        }
}
