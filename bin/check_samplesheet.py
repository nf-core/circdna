#!/usr/bin/env python
# This script is based on the example at: https://raw.githubusercontent.com/nf-core/test-datasets/viralrecon/samplesheet/samplesheet_test_illumina_amplicon.csv
#
# circdna.nf samplesheet validator (v4.7 / circdnalr routing model)
#
# Routing fields (platform x assay x datatype x concatemer x read_type):
#   platform    illumina | pacbio | ont                     -> pre-processing chain
#   assay       wgs | circleseq | rca | ciderseq | enriched -> engine set
#   datatype    gdna | eccdna                                -> background vs detection
#   pair        string (nullable)                            -> study grouping key
#   concatemer  true | false                                 -> long-read RCA TideHunter switch
#   read_type   hifi | clr | ont | pe | se                   -> alignment preset
#   enrichment  none | circleseq | mobilome | ciderseq | t7 | exov | other (annotation only)
#
# Legacy:
#   * --protocol short_read|pacbio|ont is only used as a PLATFORM fallback when a
#     row does not carry an explicit `platform` column (short_read -> illumina).
#   * A `protocol` column with short_read|long_read values is mapped the same way.
#
# Mixed tables are allowed: one file may hold Illumina Circle-seq + PacBio WGS +
# ONT RCA rows; each row is routed on its own fields.

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat nf-core/circdna samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT> <INPUT_FORMAT> <PROTOCOL> [MODE]"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")
    parser.add_argument("INPUT_FORMAT", help="'FASTQ' or 'BAM' File Format.")
    parser.add_argument("PROTOCOL", help="'short_read', 'pacbio', or 'ont' legacy platform fallback.")
    parser.add_argument("MODE", nargs="?", default="eccdna", help="'reference' or 'eccdna' (datatype backfill when row omits datatype).")
    return parser.parse_args(args)


def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception


def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)


# ---------------------------------------------------------------------------
# Valid values / defaults
# ---------------------------------------------------------------------------
VALID_DATATYPES   = ["gdna", "eccdna"]
VALID_PLATFORMS   = ["illumina", "pacbio", "ont"]
VALID_ASSAYS      = ["wgs", "circleseq", "rca", "ciderseq", "enriched"]
VALID_READ_TYPES  = ["hifi", "clr", "ont", "pe", "se"]
VALID_CONCATEMER  = ["true", "false"]
VALID_ENRICHMENT  = ["none", "circleseq", "mobilome", "ciderseq", "t7", "exov", "other"]
VALID_ENTRYPOINTS = ["cleaned_fastq", "raw_fastq", "subreads", "hifi_bam"]
LEGACY_PROTOCOLS  = ["short_read", "long_read"]   # only mapped, never routed on
PLATFORM_FROM_PROTOCOL = {"short_read": "illumina", "long_read": None, "pacbio": "pacbio", "ont": "ont"}

OUTPUT_HEADER = [
    "sample", "single_end", "fastq_1", "fastq_2", "input_bam", "entrypoint",
    "platform", "assay", "datatype", "pair", "concatemer", "read_type", "enrichment",
]


def assert_valid_combination(platform, assay, datatype, enrichment):
    """Validate the §1.3 combination table."""
    if platform == "illumina":
        if assay == "wgs" and datatype == "gdna":
            return
        if assay == "circleseq" and datatype == "eccdna":
            return
        if assay == "wgs" and datatype == "eccdna" and enrichment != "none":
            return  # short-read WGS library that is actually enriched (e.g. circSeq mislabels)
        print_error(
            "Invalid combination for illumina: assay='{}' datatype='{}'. Allowed: "
            "wgs+gdna, circleseq+eccdna (Illumina RCA must be written as circleseq).".format(assay, datatype)
        )
    else:  # pacbio / ont (long-read)
        if assay == "ciderseq":
            if datatype != "eccdna":
                print_error("Invalid combination: ciderseq requires datatype='eccdna' (got '{}').".format(datatype))
            return
        if assay == "rca":
            if datatype != "eccdna":
                print_error("Invalid combination: rca requires datatype='eccdna' (got '{}').".format(datatype))
            return
        if assay == "enriched":
            if datatype != "eccdna":
                print_error("Invalid combination: enriched requires datatype='eccdna' (got '{}').".format(datatype))
            return
        if assay == "wgs":
            if datatype != "gdna":
                print_error(
                    "Invalid combination: long-read wgs requires datatype='gdna' (got '{}'). "
                    "Un-enriched WGS cannot be used as eccDNA detection.".format(datatype)
                )
            return
        print_error(
            "Invalid assay '{}' for platform '{}'. Allowed long-read assays: wgs, rca, ciderseq, enriched.".format(
                assay, platform
            )
        )


def check_samplesheet(file_in, file_out, input_format, protocol_arg, mode="eccdna"):
    """
    Unified samplesheet checker.

    FASTQ input (short + long + mixed):
        sample,fastq_1,fastq_2,[input_bam,entrypoint],platform,assay,datatype,pair,concatemer,read_type,enrichment
    BAM input (short-read only):
        sample,bam
    """

    sample_mapping_dict = {}
    OPTIONAL_FIELDS = [
        "lane", "datatype", "data_type", "platform", "protocol", "group", "pair",
        "assay", "concatemer", "read_type", "enrichment", "input_bam", "entrypoint",
        "library_source",
    ]
    legacy_platform_fallback = PLATFORM_FROM_PROTOCOL.get(protocol_arg, None)
    default_datatype = "gdna" if mode == "reference" else "eccdna"

    with open(file_in, "r") as fin:
        header = [x.strip('"') for x in fin.readline().strip().split(",")]

        if header[:1] != ["sample"]:
            print_error("header must start with 'sample'", "Header", ",".join(header))

        has_lane          = "lane" in header
        has_datatype      = ("datatype" in header or "data_type" in header)
        has_platform      = "platform" in header
        has_protocol_col  = "protocol" in header
        has_group         = "group" in header
        has_pair          = "pair" in header
        has_assay         = "assay" in header
        has_concatemer    = "concatemer" in header
        has_read_type     = "read_type" in header
        has_enrichment    = "enrichment" in header
        has_input_bam     = "input_bam" in header
        has_entrypoint    = "entrypoint" in header
        has_library_source = "library_source" in header
        is_bam_format     = input_format == "BAM"

        if is_bam_format:
            # Short-read BAM input: sample,bam[,optional routing fields]
            if "bam" not in header:
                print_error("BAM input requires a 'bam' column", "Header", ",".join(header))

        for line in fin:
            lspl = [x.strip().strip('"') for x in line.strip().split(",")]
            if not any(lspl):
                continue

            if is_bam_format:
                if len(lspl) < 2:
                    print_error("Invalid number of columns (minimum = 2)!", "Line", line)
                row = dict(zip(header, lspl + [''] * (len(header) - len(lspl))))
                sample = row["sample"].replace(" ", "_")
                if not sample:
                    print_error("Sample entry has not been specified!", "Line", line)
                bam = row["bam"]
                if bam and bam.find(" ") != -1:
                    print_error("BAM file contains spaces!", "Line", line)
                if bam and not bam.endswith(".bam"):
                    print_error("Bam file does not have extension '.bam'!", "Line", line)
                if not bam:
                    print_error("Bam file has not been specified!", "Line", line)

                meta = build_meta(row, protocol_arg, input_format, is_short_read=True, single_end=False, default_datatype=default_datatype)
                sample_info = ["0", "", "", "", "", "", meta["platform"], meta["assay"], meta["datatype"],
                                meta["pair"], meta["concatemer"], meta["read_type"], meta["enrichment"]]
                sample_mapping_dict.setdefault(sample, []).append(sample_info)
                continue

            # ------------------------------------------------------------------
            # FASTQ input (short + long + mixed)
            # ------------------------------------------------------------------
            if len(lspl) < 2:
                print_error("Invalid number of columns (minimum = 2)!", "Line", line)
            row = dict(zip(header, lspl + [''] * (len(header) - len(lspl))))

            sample = row["sample"].replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            fastq_1    = row.get("fastq_1", "")
            fastq_2    = row.get("fastq_2", "")
            input_bam  = row.get("input_bam", "")
            entrypoint = row.get("entrypoint", "")

            # TRANSCRIPTOMIC source (circRNA etc.) is out of scope for circdna.nf
            if has_library_source and row.get("library_source", "").upper() == "TRANSCRIPTOMIC":
                print_error("library_source=TRANSCRIPTOMIC (circRNA) is not accepted by circdna.nf!", "Line", line)

            # platform: row value > legacy `protocol` column > CLI --protocol fallback
            platform = ""
            if has_platform and row.get("platform", ""):
                platform = row["platform"].strip().lower()
            elif has_protocol_col and row.get("protocol", ""):
                mapped = PLATFORM_FROM_PROTOCOL.get(row["protocol"].strip().lower(), None)
                if mapped:
                    platform = mapped
            if not platform:
                platform = legacy_platform_fallback
            if not platform:
                print_error(
                    "platform missing and no legacy --protocol fallback available! "
                    "Set the 'platform' column (illumina|pacbio|ont).",
                    "Line", line,
                )
            if platform not in VALID_PLATFORMS:
                print_error(
                    "Invalid platform '{}'! Must be one of: {}".format(platform, ", ".join(VALID_PLATFORMS)),
                    "Line", line,
                )

            is_short_read = platform == "illumina"

            # entrypoint validation (long-read only)
            if not is_short_read and entrypoint:
                if entrypoint not in VALID_ENTRYPOINTS:
                    print_error(
                        "Invalid entrypoint '{}'! Must be one of: {}".format(entrypoint, ", ".join(VALID_ENTRYPOINTS)),
                        "Line", line,
                    )
            if not is_short_read and entrypoint == "":
                entrypoint = "cleaned_fastq"

            # file validation
            for fastq in [fastq_1, fastq_2]:
                if fastq:
                    if fastq.find(" ") != -1:
                        print_error("FastQ file contains spaces!", "Line", line)
                    if not fastq.endswith(".fastq.gz") and not fastq.endswith(".fq.gz"):
                        print_error(
                            "FastQ file does not have extension '.fastq.gz' or '.fq.gz'!",
                            "Line", line,
                        )
            if input_bam:
                if input_bam.find(" ") != -1:
                    print_error("BAM file contains spaces!", "Line", line)
                if not input_bam.endswith(".bam"):
                    print_error("BAM file does not have extension '.bam'!", "Line", line)

            if is_short_read:
                if not fastq_1:
                    print_error("Illumina rows require fastq_1!", "Line", line)
            else:
                if not fastq_1 and not input_bam:
                    print_error("Either fastq_1 or input_bam must be provided for long-read rows!", "Line", line)

            single_end = (not fastq_2) or (is_short_read is False and not input_bam and fastq_2 == "")

            meta = build_meta(row, protocol_arg, input_format, is_short_read=is_short_read, single_end=single_end, default_datatype=default_datatype)
            # 长读无 assay 时留空，由 input_check 用 params.assay 回填后再断言；此处不提前报错
            if meta["assay"]:
                assert_valid_combination(meta["platform"], meta["assay"], meta["datatype"], meta["enrichment"])

            sample_info = [
                "1" if single_end else "0",
                fastq_1, fastq_2, input_bam, entrypoint,
                meta["platform"], meta["assay"], meta["datatype"],
                meta["pair"], meta["concatemer"], meta["read_type"], meta["enrichment"],
            ]
            if has_lane:
                lane = row.get("lane", "").strip()
                sample_info.append(lane)

            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)
        make_dir(out_dir)
        with open(file_out, "w") as fout:
            out_header = list(OUTPUT_HEADER)
            if has_lane:
                out_header.append("lane")
            fout.write(",".join(out_header) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                for idx, val in enumerate(sample_mapping_dict[sample]):
                    # val = [single_end, fastq_1, fastq_2, input_bam, entrypoint,
                    #        platform, assay, datatype, pair, concatemer, read_type, enrichment(, lane)]
                    row_vals = [sample] + list(val)
                    if not has_lane:
                        # avoid silent column drift; numbered suffix when no lane to keep unique
                        pass
                    fout.write(",".join(row_vals) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))


def build_meta(row, protocol_arg, input_format, is_short_read, single_end, default_datatype="eccdna"):
    """Build + validate routing fields for one row; apply §5.1 defaults."""
    meta = {}

    # platform ------------------------------------------------------------
    platform = ""
    if "platform" in row and row.get("platform", ""):
        platform = row["platform"].strip().lower()
    elif "protocol" in row and row.get("protocol", ""):
        platform = PLATFORM_FROM_PROTOCOL.get(row["protocol"].strip().lower(), None) or ""
    if not platform:
        platform = PLATFORM_FROM_PROTOCOL.get(protocol_arg, None) or ""
    meta["platform"] = platform

    # datatype ------------------------------------------------------------
    datatype = ""
    if "datatype" in row and row.get("datatype", ""):
        datatype = row["datatype"].strip().lower()
    elif "data_type" in row and row.get("data_type", ""):
        datatype = row["data_type"].strip().lower()
    if datatype and datatype not in VALID_DATATYPES:
        print_error(
            "Invalid datatype '{}'! Must be one of: {}".format(datatype, ", ".join(VALID_DATATYPES)),
            "Sample: {}".format(row.get("sample")),
        )
    if not datatype:
        datatype = default_datatype
    meta["datatype"] = datatype

    # assay ---------------------------------------------------------------
    assay = row.get("assay", "").strip().lower() if "assay" in row else ""
    if assay and assay not in VALID_ASSAYS:
        print_error(
            "Invalid assay '{}'! Must be one of: {}".format(assay, ", ".join(VALID_ASSAYS)),
            "Sample: {}".format(row.get("sample")),
        )
    if not assay:
        if is_short_read:
            # 短读无 assay：datatype=gdna→wgs，eccdna→circleseq
            assay = "wgs" if datatype == "gdna" else "circleseq"
        else:
            # 长读无 assay：不猜（尤其禁止猜 ciderseq），留给 input_check 用 params.assay 回填或报错
            assay = ""
    meta["assay"] = assay

    # pair / group ----------------------------------------------------------
    pair = ""
    if "pair" in row and row.get("pair", ""):
        pair = row["pair"].strip()
    elif "group" in row and row.get("group", ""):
        pair = row["group"].strip()
    meta["pair"] = pair

    # concatemer ------------------------------------------------------------
    concatemer = ""
    if "concatemer" in row and row.get("concatemer", ""):
        concatemer = row["concatemer"].strip().lower()
        if concatemer not in VALID_CONCATEMER:
            print_error(
                "Invalid concatemer '{}'! Must be one of: {}".format(concatemer, ", ".join(VALID_CONCATEMER)),
                "Sample: {}".format(row.get("sample")),
            )
    if not concatemer:
        # 缺省：assay=rca → true，其余 false（仅当 assay 已知）
        concatemer = "true" if assay == "rca" else "false"
    meta["concatemer"] = concatemer

    # read_type --------------------------------------------------------------
    read_type = row.get("read_type", "").strip().lower() if "read_type" in row else ""
    if read_type and read_type not in VALID_READ_TYPES:
        print_error(
            "Invalid read_type '{}'! Must be one of: {}".format(read_type, ", ".join(VALID_READ_TYPES)),
            "Sample: {}".format(row.get("sample")),
        )
    if not read_type:
        if is_short_read:
            read_type = "se" if single_end else "pe"
        else:
            # CLR 必须显式写 read_type=clr；缺省 ont→ont，pacbio→hifi
            read_type = "ont" if meta["platform"] == "ont" else "hifi"
    meta["read_type"] = read_type

    # enrichment (annotation only) -------------------------------------------
    enrichment = row.get("enrichment", "").strip().lower() if "enrichment" in row else ""
    if enrichment and enrichment not in VALID_ENRICHMENT:
        print_error(
            "Invalid enrichment '{}'! Must be one of: {}".format(enrichment, ", ".join(VALID_ENRICHMENT)),
            "Sample: {}".format(row.get("sample")),
        )
    if not enrichment:
        enrichment = "none"
    meta["enrichment"] = enrichment

    return meta


def main(args=None):
    args = parse_args(args)
    check_samplesheet(args.FILE_IN, args.FILE_OUT, args.INPUT_FORMAT, args.PROTOCOL, args.MODE)


if __name__ == "__main__":
    sys.exit(main())
