version 1.0

######################################################################################
## A pipeline for running the clipper script
######################################################################################

workflow RunClipper {
    input {
        File? aligned_bam
        File? intermediate_file
        String prefix

        File? reference_fasta
        File? reference_fai
        Boolean stream_bam = true
        String? billing_project
        Int min_clip = 20

        Int? min_dist
        Int? min_cluster_count
        Int? max_unique_breakends
        Int? max_cluster_dist

        File? mei_vcf
        String? mei_sample
        Boolean mei_all_sites = false
        String mei_info_key = "MEI_type"
        String? mei_classes
        String mei_exclude_svtypes = "DEL"
        Int mei_window = 100
        Float mei_partner_frac = 0.5

        String docker = "quay.io/ymostovoy/clipper:2.1"
    }

    parameter_meta {
        aligned_bam:            "aligned BAM or CRAM (not needed if intermediate_file is given)"
        intermediate_file:      "existing intermediate file from a previous run (pre-2.0 files are fine); skips BAM processing"
        prefix:                 "e.g. sample name"
        reference_fasta:        "reference FASTA (only needed for CRAM input)"
        reference_fai:          "index for reference_fasta"
        stream_bam:             "stream the BAM/CRAM from GCS instead of copying it to disk first; the stream is checked against the object's size and CRC32C (default=true)"
        billing_project:        "Google project to bill for requester-pays buckets"
        min_clip:               "minimum length (bp) of a clipped read end with no split-read partner to write to unaligned_clips (default=20)"
        min_dist:               "Minimum distance between where read splits align on the reference (bp) (default=1)"
        min_cluster_count:      "Minimum number of reads in a cluster (default=5)"
        max_unique_breakends:   "Maximum number of unique breakends for a cluster (default=10)"
        max_cluster_dist:       "Max distance between supplementary breakends to cluster them together (default=50)"
        mei_vcf:                "SV VCF/BCF with MEI-annotated calls for this sample; clusters explained by MEIs get FILTER=MEI_SITE/MEI_PARTNER"
        mei_sample:             "sample name in mei_vcf (needed if mei_vcf has more than one sample)"
        mei_all_sites:          "use all MEI sites in mei_vcf regardless of genotype (default=false)"
        mei_info_key:           "INFO key marking MEI calls (default=MEI_type)"
        mei_classes:            "comma-separated values of mei_info_key to treat as MEIs, e.g. ALU,LINE,SVA (default: any non-missing value)"
        mei_exclude_svtypes:    "comma-separated SVTYPEs to ignore in mei_vcf, or NONE (default=DEL)"
        mei_window:             "max distance (bp) between a cluster or partner breakend and an MEI call (default=100)"
        mei_partner_frac:       "min fraction of a cluster's reads with partners at MEI calls to flag it MEI_PARTNER (default=0.5)"
    }

    if (!defined(intermediate_file)) {
        call ProcessBam {
            input:
                aligned_bam = select_first([aligned_bam]),
                prefix = prefix,
                reference_fasta = reference_fasta,
                reference_fai = reference_fai,
                stream_bam = stream_bam,
                billing_project = billing_project,
                min_clip = min_clip,
                docker = docker
        }
    }

    call GetClusters {
        input:
            intermediate_file = select_first([intermediate_file, ProcessBam.intermediate_file]),
            prefix = prefix,
            min_dist = select_first([min_dist, 1]),
            min_cluster_count = select_first([min_cluster_count, 5]),
            max_unique_breakends = select_first([max_unique_breakends, 10]),
            max_cluster_dist = select_first([max_cluster_dist, 50]),
            mei_vcf = mei_vcf,
            mei_sample = mei_sample,
            mei_all_sites = mei_all_sites,
            mei_info_key = mei_info_key,
            mei_classes = mei_classes,
            mei_exclude_svtypes = mei_exclude_svtypes,
            mei_window = mei_window,
            mei_partner_frac = mei_partner_frac,
            docker = docker
    }

    output {
        File clipped_vcf = GetClusters.clustervcf
        File? clipped_intermediate = ProcessBam.intermediate_file
        File? unaligned_clips = ProcessBam.unaligned_clips
    }
}

task ProcessBam {
    input {
        File aligned_bam
        String prefix
        File? reference_fasta
        File? reference_fai
        Boolean stream_bam = true
        String? billing_project
        Int min_clip = 20
        Int threads = 4
        String docker
        RuntimeAttr? runtime_attr_override
    }

    parameter_meta {
        aligned_bam: { description: "aligned BAM or CRAM", localization_optional: true }
    }

    # streaming only needs room for the intermediate file
    Int disk_size = if stream_bam then 20 else ceil(size(aligned_bam, "GB")) + 20

    command <<<
        set -euo pipefail

        BAM="~{aligned_bam}"
        if [[ "$BAM" == gs://* ]]; then
            if ~{stream_bam}; then
                # gcloud storage cat doesn't checksum, and a truncated stream looks complete to processbam:
                # check the bytes against the object's size and CRC32C
                gcloud storage objects describe ~{"--billing-project=" + billing_project} "$BAM" \
                    --format='value(size,crc32c_hash)' > object_checksum.txt
                read -r SIZE CRC < object_checksum.txt
                echo "Expecting $SIZE bytes with CRC32C $CRC" >&2
                gcloud storage cat ~{"--billing-project=" + billing_project} "$BAM" \
                    | python3 /opt/bigclipper/bigclipper_verifystream.py --size "$SIZE" --crc32c "$CRC" \
                    | python3 /opt/bigclipper/bigclipper_processbam.py - -o "~{prefix}" -d . -t ~{threads} -m ~{min_clip} ~{"-r " + reference_fasta}
            else
                # gcloud storage cp validates the download's checksum itself
                gcloud storage cp ~{"--billing-project=" + billing_project} "$BAM" local_alignments
                python3 /opt/bigclipper/bigclipper_processbam.py local_alignments -o "~{prefix}" -d . -t ~{threads} -m ~{min_clip} ~{"-r " + reference_fasta}
                rm local_alignments
            fi
        else
            python3 /opt/bigclipper/bigclipper_processbam.py "$BAM" -o "~{prefix}" -d . -t ~{threads} -m ~{min_clip} ~{"-r " + reference_fasta}
        fi
    >>>

    output {
        File intermediate_file = "~{prefix}_intermediate.bed"
        File unaligned_clips = "~{prefix}_unaligned_clips.bed"
    }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          threads,
        mem_gb:             8,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        1,
        docker:             docker
    }

    RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
    runtime {
        cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
        memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
        disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " SSD"
        bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
        preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
        maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
        docker:                 select_first([runtime_attr.docker,            default_attr.docker])
    }
}

task GetClusters {
    input {
        File intermediate_file
        String prefix
        Int min_dist
        Int min_cluster_count
        Int max_unique_breakends
        Int max_cluster_dist
        File? mei_vcf
        String? mei_sample
        Boolean mei_all_sites = false
        String mei_info_key = "MEI_type"
        String? mei_classes
        String mei_exclude_svtypes = "DEL"
        Int mei_window = 100
        Float mei_partner_frac = 0.5
        String docker
        RuntimeAttr? runtime_attr_override
    }

    Int disk_size = 3*ceil(size(intermediate_file, "GB") + size(mei_vcf, "GB")) + 10
    String outname = "~{prefix}_clipped_reads_d~{min_dist}_c~{min_cluster_count}_s~{max_cluster_dist}_u~{max_unique_breakends}.vcf"

    command <<<
        set -euo pipefail

        python3 /opt/bigclipper/bigclipper_getclusters.py ~{intermediate_file} \
            -d ~{min_dist} -c ~{min_cluster_count} -s ~{max_cluster_dist} -u ~{max_unique_breakends} \
            -o "~{outname}" \
            ~{"--mei_vcf " + mei_vcf} ~{"--mei_sample " + mei_sample} ~{if mei_all_sites then "--mei_all_sites" else ""} \
            ~{"--mei_classes " + mei_classes} --mei_exclude_svtypes ~{mei_exclude_svtypes} \
            --mei_info_key ~{mei_info_key} --mei_window ~{mei_window} --mei_partner_frac ~{mei_partner_frac}
    >>>

    output { File clustervcf = outname }

    #########################
    RuntimeAttr default_attr = object {
        cpu_cores:          1,
        mem_gb:             4,
        disk_gb:            disk_size,
        boot_disk_gb:       10,
        preemptible_tries:  3,
        max_retries:        0,
        docker:             docker
    }

     RuntimeAttr runtime_attr = select_first([runtime_attr_override, default_attr])
     runtime {
         cpu:                    select_first([runtime_attr.cpu_cores,         default_attr.cpu_cores])
         memory:                 select_first([runtime_attr.mem_gb,            default_attr.mem_gb]) + " GiB"
         disks: "local-disk " +  select_first([runtime_attr.disk_gb,           default_attr.disk_gb]) + " HDD"
         bootDiskSizeGb:         select_first([runtime_attr.boot_disk_gb,      default_attr.boot_disk_gb])
         preemptible:            select_first([runtime_attr.preemptible_tries, default_attr.preemptible_tries])
         maxRetries:             select_first([runtime_attr.max_retries,       default_attr.max_retries])
         docker:                 select_first([runtime_attr.docker,            default_attr.docker])
     }
}

struct RuntimeAttr {
    Float? mem_gb
    Int? cpu_cores
    Int? disk_gb
    Int? boot_disk_gb
    Int? preemptible_tries
    Int? max_retries
    String? docker
}
