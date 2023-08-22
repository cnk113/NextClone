#!/bin/bash nextflow

include { 
    ngs_trim_reads;
    ngs_filter_reads;
    ngs_prepare_clone_barcodes;
    ngs_map_barcodes;
    ngs_collapse_barcodes
} from "./modules/extract_ngs_barcodes"

include { 
    sc_get_unmapped_reads;
    sc_remove_low_qual_reads;
    sc_split_unmapped_reads;
    sc_map_unmapped_reads;
    sc_merge_barcodes 
} from "./modules/extract_sc_clone_barcodes"

workflow {

    if (params.mode == 'NGS') {
        ch_barcode_chunks = Channel.fromPath(params.ngs_fastq_files) | 
            ngs_trim_reads |
            ngs_filter_reads |
            ngs_prepare_clone_barcodes
        
        ch_barcode_mappings = ngs_map_barcodes(ch_barcode_chunks.flatten())
        ngs_collapse_barcodes(ch_barcode_mappings.collect())

    } 
    
    if (params.mode == 'scRNAseq') {
        ch_unmapped_fastas = Channel.fromPath(params.possorted_aligned_reads_bam_file) | 
            sc_get_unmapped_reads |
            sc_remove_low_qual_reads |
            sc_split_unmapped_reads
        
        ch_mapped_fastas = sc_map_unmapped_reads(ch_unmapped_fastas[0].flatten())
        sc_merge_barcodes(ch_mapped_fastas.collect())
    }
}