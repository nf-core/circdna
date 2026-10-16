#!/usr/bin/env python3
"""Regenerate CircleSeeker CI fixtures from the upstream simulator.

Provenance
----------
Upstream simulator (CircleSeeker v1.6.1)::

    https://github.com/YaoxinBio/CircleSeeker/blob/v1.6.1/tests/simulation/eccdna_simulator.py

CI adaptation
-------------
The official UeccDNA simulator emits single-copy circular sequences
(``copy_number=1.0``, no tandem concatemer). TideHunter / CtcReads-Core needs
RCA-like tandem repeats to produce Confirmed calls on tiny fixtures, so this
script:

1. Runs the upstream ``SimulationConfig`` / ``run_simulation`` with CI-scale
   parameters (``num_mecc=0`` so ``auto_scale_genome`` does not inflate the
   reference to 1 Mb).
2. Expands each simulated Uecc read into ``TANDEM_COPIES`` tandem repeats and
   rewrites the FASTA header copy-number field.

Canonical published copies (after review) live in nf-core/test-datasets::

    data/genomics/homo_sapiens/pacbio/fasta/circleseeker_test_reads.fa
    data/genomics/homo_sapiens/pacbio/fasta/circleseeker_test_ref.fa

Until merged, SiYangming/test-datasets@circleseeker-testdata hosts the same
paths. Local copies under ``tests/fixtures/`` are kept for offline regeneration
checks and checksums in ``tests/fixtures/MANIFEST.json``.

Usage
-----
::

    # Uses CIRCLESEEKER_SRC or ../CircleSeeker next to this monorepo checkout
    python tests/generate_testdata.py

    python tests/generate_testdata.py \\
        --circle-seeker-src /path/to/CircleSeeker \\
        --out tests/fixtures
"""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

# ---------------------------------------------------------------------------
# Pinned generation parameters (do not change without bumping MANIFEST)
# ---------------------------------------------------------------------------
CIRCLESEEKER_TAG = "v1.6.1"
CIRCLESEEKER_SIMULATOR_REL = "tests/simulation/eccdna_simulator.py"
SEED = 42
NUM_UECC = 3
NUM_MECC = 0
NUM_CECC = 0
NUM_CHR = 1
CHR_LENGTH = 8000
UECC_LENGTH_RANGE = (300, 800)
TANDEM_COPIES = 3
AUTO_SCALE_GENOME = False

READS_NAME = "circleseeker_test_reads.fa"
REF_NAME = "circleseeker_test_ref.fa"


def _default_circle_seeker_src() -> Path | None:
    env = Path(__file__).resolve()
    # Prefer an explicit sibling checkout used in this workspace layout
    candidates = [
        Path("/Users/siyangming/nextflow_nf_core/CircleSeeker"),
        env.parents[5] / "CircleSeeker" if len(env.parents) > 5 else None,
        Path.cwd() / "CircleSeeker",
    ]
    for c in candidates:
        if c and (c / CIRCLESEEKER_SIMULATOR_REL).is_file():
            return c
    return None


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def _write_tandem_reads(src_reads: Path, dst_reads: Path, tandem: int) -> int:
    records: list[tuple[str, str]] = []
    name: str | None = None
    seq: list[str] = []
    for line in src_reads.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(seq)))
            name = line[1:]
            seq = []
        else:
            seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))

    with dst_reads.open("w") as fh:
        for hdr, s in records:
            parts = hdr.split("|")
            if len(parts) >= 4:
                parts[2] = str(len(s))
                parts[-1] = f"{tandem}.0"
                new_hdr = "|".join(parts)
            else:
                new_hdr = f"{hdr}|tandem|{len(s)}|{tandem}.0"
            tandem_seq = s * tandem
            fh.write(f">{new_hdr}\n")
            for i in range(0, len(tandem_seq), 80):
                fh.write(tandem_seq[i : i + 80] + "\n")
    return len(records)


def generate(circle_seeker_src: Path, out_dir: Path) -> dict:
    sim_path = circle_seeker_src / CIRCLESEEKER_SIMULATOR_REL
    if not sim_path.is_file():
        raise FileNotFoundError(
            f"Upstream simulator not found: {sim_path}\n"
            f"Clone https://github.com/YaoxinBio/CircleSeeker @ {CIRCLESEEKER_TAG} "
            "and pass --circle-seeker-src"
        )

    sys.path.insert(0, str(sim_path.parent))
    from eccdna_simulator import SimulationConfig, run_simulation  # type: ignore

    work = out_dir / "_simulator_workdir"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    cfg = SimulationConfig(
        output_dir=str(work),
        num_uecc=NUM_UECC,
        num_mecc=NUM_MECC,
        num_cecc=NUM_CECC,
        num_chromosomes=NUM_CHR,
        chr_length=CHR_LENGTH,
        seed=SEED,
        auto_scale_genome=AUTO_SCALE_GENOME,
        uecc_length_range=UECC_LENGTH_RANGE,
    )
    run_simulation(cfg)

    out_dir.mkdir(parents=True, exist_ok=True)
    ref_dst = out_dir / REF_NAME
    reads_dst = out_dir / READS_NAME
    shutil.copyfile(work / "reference.fa", ref_dst)
    n_reads = _write_tandem_reads(work / "simulated_reads.fa", reads_dst, TANDEM_COPIES)

    # Keep upstream artefacts for audit
    for name in ("ground_truth.csv", "simulation_summary.json", "simulated_reads.fa"):
        src = work / name
        if src.is_file():
            shutil.copyfile(src, out_dir / f"upstream_{name}")

    manifest = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "circle_seeker_tag": CIRCLESEEKER_TAG,
        "simulator_relpath": CIRCLESEEKER_SIMULATOR_REL,
        "circle_seeker_src_basename": circle_seeker_src.resolve().name,
        "parameters": {
            "seed": SEED,
            "num_uecc": NUM_UECC,
            "num_mecc": NUM_MECC,
            "num_cecc": NUM_CECC,
            "num_chromosomes": NUM_CHR,
            "chr_length": CHR_LENGTH,
            "uecc_length_range": list(UECC_LENGTH_RANGE),
            "auto_scale_genome": AUTO_SCALE_GENOME,
            "tandem_copies": TANDEM_COPIES,
        },
        "ci_adaptation": (
            "Each upstream Uecc read (copy_number=1.0) is expanded to "
            f"{TANDEM_COPIES} tandem copies so CtcReads-Core/TideHunter has RCA-like signal."
        ),
        "outputs": {
            REF_NAME: {"bytes": ref_dst.stat().st_size, "md5": _md5(ref_dst)},
            READS_NAME: {
                "bytes": reads_dst.stat().st_size,
                "md5": _md5(reads_dst),
                "n_reads": n_reads,
            },
        },
        "test_datasets_paths": [
            f"data/genomics/homo_sapiens/pacbio/fasta/{READS_NAME}",
            f"data/genomics/homo_sapiens/pacbio/fasta/{REF_NAME}",
        ],
    }
    (out_dir / "MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    shutil.rmtree(work)
    return manifest


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        "--circle-seeker-src",
        type=Path,
        default=None,
        help=f"Path to YaoxinBio/CircleSeeker checkout at {CIRCLESEEKER_TAG}",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path(__file__).resolve().parent / "fixtures",
        help="Output directory for fixtures + MANIFEST.json",
    )
    args = parser.parse_args()

    src = args.circle_seeker_src or _default_circle_seeker_src()
    if src is None:
        print(
            "ERROR: CircleSeeker source not found. Pass --circle-seeker-src "
            f"pointing at a checkout of {CIRCLESEEKER_TAG}.",
            file=sys.stderr,
        )
        return 1

    manifest = generate(src, args.out)
    print(json.dumps(manifest, indent=2))
    print(f"\nWrote fixtures to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
