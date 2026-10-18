#!/usr/bin/env python
"""
metadata_to_samplesheet.py — map the representative metadata.csv into a runnable
circdna.nf samplesheet using the §1.2 routing rules.

Columns (from plan §7):
    sample,fastq_1,fastq_2,input_bam,entrypoint,platform,assay,datatype,pair,concatemer,read_type,enrichment

Rules (plan §1.2):
    * workflow=circrna 或 library_source=TRANSCRIPTOMIC → rejected
    * platform: ILLUMINA→illumina, PACBIO_SMRT→pacbio, OXFORD_NANOPORE→ont
    * analysis_role=reference_genome → datatype=gdna ; primary_eccDNA → eccdna
    * experimental_type=CIDER-seq → assay=ciderseq, enrichment=ciderseq
    * long-read + enriched + RCA + non-CIDER → assay=rca
    * long-read + enriched + no RCA → assay=enriched
    * short-read + enriched → assay=circleseq
    * un-enriched WGS / genome-skimming / size fractionation → assay=wgs
    * amplification_method=RCA on illumina → still circleseq
    * notes contain linearized / T7 debranching → concatemer=false
    * long-read RCA without linearization → concatemer=true
    * assay in (wgs, enriched) → concatemer=false
    * PacBio HiFi → read_type=hifi ; PacBio RS/RS II → read_type=clr ; ONT → ont ; Illumina → pe/se
    * pair ← bioproject (or sra_study)

Usage:
    python bin/metadata_to_samplesheet.py samplesheets/metadata.csv -o samplesheets/circdna_representative.csv \
        --root /data1/users/siyangming/PlanteccDNADB/eccDNA
"""
import argparse
import csv
import sys

OUTPUT_COLS = ["sample", "fastq_1", "fastq_2", "input_bam", "entrypoint",
                "platform", "assay", "datatype", "pair", "concatemer", "read_type", "enrichment"]

EXPERIMENTAL_TO_ASSAY = {
    "CIDER-seq": ("ciderseq", "ciderseq"),
}

SHORT_ENRICHMENT_TO_ENRICHMENT = {
    "mobilome-seq": "mobilome",
    "Mobilome-seq": "mobilome",
    "circSeq": "circleseq",
    "eccDNA-seq": "circleseq",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("metadata", help="metadata.csv (representative set)")
    p.add_argument("-o", "--out", default="samplesheets/circdna_representative.csv")
    p.add_argument("--root", default="/data1/users/siyangming/PlanteccDNADB/eccDNA",
                    help="Base directory holding per-species fastq folders")
    p.add_argument("--reject-file", default=None, help="Optional: write rejected rows here")
    p.add_argument("--full-master", default=None,
                    help="Legacy master long-read samplesheet to migrate (adds platform/assay/datatype/... columns per row)")
    p.add_argument("--sra", default="samplesheets/SraRunInfo_eccDNA_all2.csv",
                    help="SRA Run Info for platform/strategy inference in --full-master mode")
    p.add_argument("--species-filter", default=None,
                    help="Only emit rows whose species dir matches (e.g. Arabidopsis_thaliana); used for per-species files")
    return p.parse_args()


def species_dir(species):
    return species.replace(" ", "_")


def infer(meta):
    """Return routing dict per §1.2, or None to reject."""
    if meta.get("workflow") == "circrna":
        return None
    if meta.get("library_source", "").upper() == "TRANSCRIPTOMIC":
        return None

    platform_raw = meta.get("platform", "").upper()
    platform = {"ILLUMINA": "illumina", "PACBIO_SMRT": "pacbio", "OXFORD_NANOPORE": "ont"}.get(platform_raw)
    if not platform:
        return None

    role = meta.get("analysis_role", "")
    datatype = "gdna" if role == "reference_genome" else "eccdna"
    exp_type = meta.get("experimental_type", "")

    notes = (meta.get("notes") or "")
    layout = meta.get("layout", "")
    read_length = (meta.get("read_length") or "").lower()

    # experimental type override (CIDER)
    if exp_type in EXPERIMENTAL_TO_ASSAY:
        assay, enrichment = EXPERIMENTAL_TO_ASSAY[exp_type]
        concatemer = "false"
        # PacBio RS / RS II → clr（§1.2：CLR 必须显式写）
        if platform == "ont":
            read_type = "ont"
        elif "RS" in (notes or "") or "RS" in (meta.get("sequencing_generation") or ""):
            read_type = "clr"
        else:
            read_type = "hifi"
        return dict(platform=platform, assay=assay, datatype=datatype, concatemer=concatemer,
                    read_type=read_type, enrichment=enrichment)

    enriched = (meta.get("eccDNA_enrichment") or "").strip().lower() in ("yes", "true", "1")
    is_long = platform in ("pacbio", "ont")

    if enriched:
        if is_long:
            amplification = (meta.get("circular_dna_amplification") or "").upper()
            has_rca = "RCA" in amplification or "RCA" in (meta.get("amplification_method") or "").upper()
            linearized = ("linearized" in notes.lower() or "debranch" in notes.lower())
            if has_rca:
                assay = "rca"
                concatemer = "false" if linearized else "true"
            else:
                assay = "enriched"   # 无 RCA 富集（如 T7）
                concatemer = "false"
            enrichment = "ciderseq" if assay == "ciderseq" else "other"
            if "mobilome" in (meta.get("eccDNA_enrichment_method") or "").lower():
                enrichment = "mobilome"
            read_type = "ont" if platform == "ont" else ("hifi" if "hifi" in read_length else "clr")
        else:
            assay = "circleseq"
            concatemer = "false"
            method = meta.get("eccDNA_enrichment_method") or ""
            enrichment = "circleseq"
            for k, v in SHORT_ENRICHMENT_TO_ENRICHMENT.items():
                if k.lower() in method.lower():
                    enrichment = v
                    break
            if "T7" in method:
                enrichment = "t7"
            elif "Exonuclease" in method or "ExoV" in method:
                enrichment = "exov"
            read_type = "se" if layout == "SINGLE" else "pe"
    else:
        assay = "wgs"
        concatemer = "false"
        enrichment = "none"
        read_type = "se" if layout == "SINGLE" else "pe"
        if is_long:
            read_type = "ont" if platform == "ont" else ("hifi" if "hifi" in read_length else "clr")

    return dict(platform=platform, assay=assay, datatype=datatype, concatemer=concatemer,
                read_type=read_type, enrichment=enrichment)


def fastq_paths(root, meta, read_type):
    """Build fastq_1/fastq_2 paths following the existing convention."""
    run = meta.get("run_id", "")
    sp = species_dir(meta.get("species", ""))
    base = "{}/{}/{}".format(root.rstrip("/"), sp, run)
    if read_type in ("pe",):
        return base + "_1.fastq.gz", base + "_2.fastq.gz"
    return base + ".fastq.gz", ""


def load_sra_runinfo(path):
    """run_id -> dict(Platform, Model, LibraryStrategy, BioProject)."""
    info = {}
    try:
        with open(path) as fh:
            for row in csv.DictReader(fh):
                run = row.get("Run") or row.get("run_id") or ""
                if run:
                    info[run] = row
    except FileNotFoundError:
        sys.stderr.write("WARN: SRA Run Info {} not found; --full-master inference limited\n".format(path))
    return info


def infer_from_sra(sra_row):
    """SRA-based fallback inference (only for runs NOT in metadata.csv)."""
    platform_raw = (sra_row.get("Platform") or "").upper()
    platform = {"ILLUMINA": "illumina", "PACBIO_SMRT": "pacbio", "OXFORD_NANOPORE": "ont"}.get(platform_raw)
    if not platform:
        return None
    strategy = (sra_row.get("LibraryStrategy") or "").upper()
    model = sra_row.get("Model") or ""
    if platform != "illumina":
        if strategy == "WGS":
            assay, datatype = "wgs", "gdna"
        else:
            # AMPLICON / OTHER on long-read: enriched circular-DNA detection
            assay, datatype = "rca", "eccdna"
        concatemer = "true" if assay == "rca" else "false"
        read_type = "ont" if platform == "ont" else ("clr" if "RS" in model.upper() else "hifi")
        enrichment = "none" if assay == "wgs" else "other"
    else:
        assay = "wgs" if strategy == "WGS" else "circleseq"
        datatype = "gdna" if assay == "wgs" else "eccdna"
        concatemer = "false"
        read_type = "se" if (sra_row.get("LibraryLayout") or "").upper() == "SINGLE" else "pe"
        enrichment = "none" if assay == "wgs" else "circleseq"
    return dict(platform=platform, assay=assay, datatype=datatype, concatemer=concatemer,
                read_type=read_type, enrichment=enrichment)


def main():
    args = parse_args()

    # authoritative metadata mapping (representative rows)
    auth = {}
    with open(args.metadata) as fh:
        for meta in csv.DictReader(fh):
            run = meta.get("run_id") or ""
            routed = infer(meta)
            if routed and run:
                auth[run] = routed

    if args.full_master:
        migrate_master(args, auth)
        return

    rejected = []
    rows = []
    with open(args.metadata) as fh:
        for meta in csv.DictReader(fh):
            routed = infer(meta)
            if routed is None:
                rejected.append(meta.get("run_id") or meta.get("sample_id") or "?")
                continue
            fq1, fq2 = fastq_paths(args.root, meta, routed["read_type"])
            rows.append({
                "sample": meta["run_id"],
                "fastq_1": fq1,
                "fastq_2": fq2,
                "input_bam": "",
                "entrypoint": "cleaned_fastq" if routed["platform"] != "illumina" else "",
                "platform": routed["platform"],
                "assay": routed["assay"],
                "datatype": routed["datatype"],
                "pair": (meta.get("bioproject") or meta.get("sra_study") or "").strip(),
                "concatemer": routed["concatemer"],
                "read_type": routed["read_type"],
                "enrichment": routed["enrichment"],
            })

    write_rows(rows, args.out)
    sys.stderr.write("metadata_to_samplesheet: wrote {} rows to {}; rejected {}\n".format(
        len(rows), args.out, len(rejected)))
    if rejected:
        sys.stderr.write("  rejected: {}\n".format(", ".join(sorted(rejected))))
    if args.reject_file:
        with open(args.reject_file, "w") as fh:
            fh.write("\n".join(sorted(rejected)) + "\n")


def migrate_master(args, auth):
    """Migrate a legacy long-read samplesheet (sample,fastq_1,fastq_2) to the routing schema."""
    sra = load_sra_runinfo(args.sra)
    rows = []
    skipped = []
    with open(args.full_master) as fh:
        for row in csv.DictReader(fh):
            sample = row.get("sample", "").strip()
            fq1 = row.get("fastq_1", "").strip()
            fq2 = row.get("fastq_2", "").strip()
            if not sample:
                continue
            # species filter (by path segment, e.g. Arabidopsis_thaliana)
            if args.species_filter and args.species_filter not in fq1:
                continue

            if sample in auth:
                routed = auth[sample]
                pair = ""
            elif sample in sra:
                routed = infer_from_sra(sra[sample])
                if not routed:
                    skipped.append(sample)
                    continue
                pair = (sra[sample].get("BioProject") or "").strip()
            else:
                skipped.append(sample)
                sys.stderr.write("WARN: no metadata/SRA row for {} — skipped\n".format(sample))
                continue

            rows.append({
                "sample": sample,
                "fastq_1": fq1,
                "fastq_2": fq2,
                "input_bam": "",
                "entrypoint": "cleaned_fastq",
                "platform": routed["platform"],
                "assay": routed["assay"],
                "datatype": routed["datatype"],
                "pair": pair,
                "concatemer": routed["concatemer"],
                "read_type": routed["read_type"],
                "enrichment": routed["enrichment"],
            })

    write_rows(rows, args.out)
    sys.stderr.write("metadata_to_samplesheet: migrated {} rows to {}; skipped {}\n".format(
        len(rows), args.out, len(skipped)))
    if skipped:
        sys.stderr.write("  skipped: {}\n".format(", ".join(sorted(skipped))))


def write_rows(rows, out_path):
    with open(out_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUTPUT_COLS)
        w.writeheader()
        for r in sorted(rows, key=lambda x: x["sample"]):
            w.writerow(r)


if __name__ == "__main__":
    sys.exit(main())
