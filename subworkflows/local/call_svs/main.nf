include { DEBREAK   } from '../debreak/main'
include { HIFICNV   } from '../hificnv/main'
include { SAWFISH   } from '../sawfish/main'
include { SEVERUS   } from '../severus/main'
include { SNIFFLES1 } from '../sniffles1/main'
include { SNIFFLES  } from '../sniffles/main'
include { SNIFFLES as SNIFFLES_JOINT } from '../snifflesJoint/main'

workflow CALL_SVS {
    take:
    ch_bam_bai                              // channel: [ val(meta), path(bam), path(bai) ]
    ch_tandem_repeats                       // channel: [ val(meta), path(bed) ]
    ch_snvs                                 // channel: [ val(meta), path(vcf) ]
    ch_fasta                                // channel: [ val(meta), path(fasta) ]
    ch_expected_xy_bed                      // channel: [ val(meta), path(bed) ]
    ch_expected_xx_bed                      // channel: [ val(meta), path(bed) ]
    ch_exclude_bed                          // channel: [ val(meta), path(bed) ]
    sv_callers_to_run                       //    List: [ 'caller1', 'caller2', ... ]
    force_sawfish_joint_call_single_samples //    bool: Force joint-calling with Sawfish even for single samples
    create_hificnv_maf_track                //    bool: Should we create a MAF track for HiFiCNV calls?
    create_sawfish_maf_track                //    bool: Should we create a MAF track for Sawfish calls?

    main:
    ch_sv_calls = channel.empty()

    if (sv_callers_to_run.contains('sniffles')) {
        ch_family_counts = ch_bam_bai
            .map { meta, _bam, _bai -> [meta.family_id, meta.id] }
            .groupTuple()
            .map { family_id, sample_ids -> [family_id, sample_ids.unique().size()] }

        ch_bams = ch_bam_bai
            .map { meta, bam, bai -> [meta.family_id, meta, bam, bai] }
            .combine(ch_family_counts, by: 0)

        ch_bams_branched = ch_bams.branch { _family_id, _meta, _bam, _bai, family_size ->
            single: family_size == 1
            joint: family_size > 1
        }

        SNIFFLES(
            ch_bams_branched.single.map { _family_id, meta, bam, bai, _family_size ->
                [meta, bam, bai]
            },
            ch_fasta,
            ch_tandem_repeats,
        )

        SNIFFLES_JOINT(
            ch_bams_branched.joint.map { _family_id, meta, bam, bai, _family_size ->
                [meta, bam, bai]
            },
            ch_fasta,
            ch_tandem_repeats,
        )

        ch_sv_calls = ch_sv_calls
            .mix(SNIFFLES.out.vcf.map { meta, vcf -> [meta, vcf, []] })
            .mix(SNIFFLES_JOINT.out.vcf.map { meta, vcf -> [meta, vcf, []]})
    }

    if (sv_callers_to_run.contains('sniffles1')) {
        SNIFFLES1(
            ch_bam_bai
        )
        ch_sv_calls = ch_sv_calls.mix(SNIFFLES1.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('severus')) {
        SEVERUS(
            ch_bam_bai,
            ch_tandem_repeats,
        )
        ch_sv_calls = ch_sv_calls.mix(SEVERUS.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('debreak')) {
        DEBREAK(
            ch_bam_bai,
            ch_fasta,
        )
        ch_sv_calls = ch_sv_calls.mix(DEBREAK.out.vcf.map { meta, vcf -> [meta, vcf, []] })
    }

    if (sv_callers_to_run.contains('hificnv')) {
        HIFICNV(
            ch_bam_bai,
            ch_snvs,
            ch_fasta,
            ch_expected_xy_bed,
            ch_expected_xx_bed,
            ch_exclude_bed,
            create_hificnv_maf_track,
        )
        ch_sv_calls = ch_sv_calls.mix(
            HIFICNV.out.vcf.join(HIFICNV.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        )
    }

    if (sv_callers_to_run.contains('sawfish')) {
        SAWFISH(
            ch_bam_bai,
            ch_snvs,
            ch_fasta,
            ch_expected_xy_bed,
            ch_expected_xx_bed,
            ch_exclude_bed,
            create_sawfish_maf_track,
            force_sawfish_joint_call_single_samples,
        )
        ch_sv_calls = ch_sv_calls.mix(
            SAWFISH.out.vcf.join(SAWFISH.out.tbi, failOnMismatch: true, failOnDuplicate: true)
        )
    }

    emit:
    sv_calls                           = ch_sv_calls                                                                                      // channel: [ val(meta), path(vcf), path(tbi) ] — tbi is [] for callers with skip_vep_prep: false
    hificnv_depth                      = sv_callers_to_run.contains('hificnv') ? HIFICNV.out.depth : channel.empty()                      // channel: [ val(meta), path(bw) ]
    hificnv_copynum                    = sv_callers_to_run.contains('hificnv') ? HIFICNV.out.copynum : channel.empty()                    // channel: [ val(meta), path(bedgraph) ]
    hificnv_maf                        = sv_callers_to_run.contains('hificnv') ? HIFICNV.out.maf : channel.empty()                        // channel: [ val(meta), path(bw) ]
    sawfish_depth_bw                   = sv_callers_to_run.contains('sawfish') ? SAWFISH.out.depth_bw : channel.empty()                   // channel: [ val(meta), path(bw) ]
    sawfish_copynum_bedgraph           = sv_callers_to_run.contains('sawfish') ? SAWFISH.out.copynum_bedgraph : channel.empty()           // channel: [ val(meta), path(bedgraph) ]
    sawfish_gc_bias_corrected_depth_bw = sv_callers_to_run.contains('sawfish') ? SAWFISH.out.gc_bias_corrected_depth_bw : channel.empty() // channel: [ val(meta), path(bw) ]
    sawfish_maf_bw                     = sv_callers_to_run.contains('sawfish') ? SAWFISH.out.maf_bw : channel.empty()                     // channel: [ val(meta), path(bw) ]
}
