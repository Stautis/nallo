include { SNIFFLES as RUN_SNIFFLES } from '../../../modules/nf-core/sniffles/main.nf'
include { SNIFFLES as RUN_SNIFFLES_JOINT } from '../../../modules/nf-core/sniffles/main.nf'

workflow SNIFFLES {
    take:
    ch_bam_bai        // channel: [mandatory] [ val(meta), path(bam), path(bai) ]
    ch_fasta          // channel: [mandatory] [ val(meta), path(fasta) ]
    ch_tandem_repeats // channel:  [optional] [ val(meta), path(bed) ]

    main:
    RUN_SNIFFLES(
        ch_bam_bai.map { meta, bam, bai -> [meta, bam, bai, [], []] },
        ch_fasta,
        ch_tandem_repeats,
        false,
        true,
    )

    ch_family_snfs = RUN_SNIFFLES.out.snf
        .map { meta, snf -> [[id: meta.family_id, family_id: meta.family_id], [meta.id, snf]] }
        .groupTuple()
        .map { family_meta, id_snf_pairs ->
            def sorted_pairs = id_snf_pairs.sort { left, right ->
                left[0] <=> right[0]
            }
            [
                family_meta + [sample_ids: sorted_pairs.collect { it[0] }],
                sorted_pairs.collect { _id, snf -> snf },
            ]
        }

    RUN_SNIFFLES_JOINT(
        ch_family_snfs.map { meta, snfs -> [meta, snfs, [], [], []] },
        ch_fasta,
        ch_tandem_repeats,
        true,
        false,
    )

    emit:
    vcf = RUN_SNIFFLES_JOINT.out.vcf.map { meta, vcf -> [meta + [sv_caller: 'sniffles', needs_reheader: false, needs_header_repair: true, skip_vep_prep: false], vcf] } // channel: [ val(meta), path(vcf) ]
}
