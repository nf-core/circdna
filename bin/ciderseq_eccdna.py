#!/usr/bin/env python
"""CIDER-Seq2 eccDNA detection wrapper (cs-eccDNA.py).

cs-eccDNA.py hard-codes relative 'input/' and 'output/' directories and a
tab-separated list file.  This wrapper stages the (deconcat fasta, stat) pairs
found in --input_dir (paired by file stem) into an `input/` directory, writes
the list file, runs cs-eccDNA.py against the host genome and copies the
resulting candidate fasta files into --outdir.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile


def pair_files(input_dir):
    """Pair deconcat fasta/stat files in a directory by file stem."""
    pairs = []
    for fa in sorted(glob.glob(os.path.join(input_dir, '*.fa')) + glob.glob(os.path.join(input_dir, '*.fasta'))):
        stem = os.path.splitext(os.path.basename(fa))[0]
        stat = os.path.join(input_dir, stem + '.stat')
        pairs.append((fa, stat if os.path.isfile(stat) else None))
    return pairs


def main():
    ap = argparse.ArgumentParser(description='CIDER-Seq2 eccDNA detection (cs-eccDNA.py)')
    ap.add_argument('--genome', required=True, help='host genome FASTA file')
    ap.add_argument('--input_dir', required=True, help='directory with deconcat fasta/stat file pairs')
    ap.add_argument('--outdir', required=True, help='output directory for candidate fasta files')
    ap.add_argument('--blast_threads', type=int, default=4, help='threads for blastn (default: 4)')
    ap.add_argument('--gap_window', type=int, default=150, help='gap window for reordering (default: 150)')
    args = ap.parse_args()

    exe = shutil.which('cs-eccDNA.py')
    if exe is None:
        sys.exit('[ciderseq_eccdna] cs-eccDNA.py not found in PATH')

    pairs = pair_files(args.input_dir)
    if not pairs:
        sys.exit('[ciderseq_eccdna] no fasta files found in %s' % args.input_dir)

    os.makedirs(args.outdir, exist_ok=True)

    staging = tempfile.mkdtemp(prefix='cider_eccdna_')
    os.makedirs(os.path.join(staging, 'input'))
    os.makedirs(os.path.join(staging, 'output'))

    list_lines = []
    for i, (fa, st) in enumerate(pairs):
        base = 'sample%d' % (i + 1)
        fa_name = base + '.fasta'
        st_name = base + '.stat'
        shutil.copy(fa, os.path.join(staging, 'input', fa_name))
        if st:
            shutil.copy(st, os.path.join(staging, 'input', st_name))
            list_lines.append('input/%s\tinput/%s' % (fa_name, st_name))
        else:
            list_lines.append('input/%s' % fa_name)

    list_file = os.path.join(staging, 'input.list')
    with open(list_file, 'wt') as fh:
        fh.write('\n'.join(list_lines) + '\n')

    genome = os.path.abspath(args.genome)
    cmd = [sys.executable, exe, '--blast_threads', str(args.blast_threads),
            '--gap_window', str(args.gap_window), list_file, genome]
    sys.stderr.write('[ciderseq_eccdna] %s\n' % ' '.join(cmd))
    subprocess.run(cmd, cwd=staging, check=True)

    candidates = sorted(glob.glob(os.path.join(staging, 'output', '*.fasta')))
    if not candidates:
        sys.exit('[ciderseq_eccdna] no candidate fasta files produced in output/')

    for cand in candidates:
        dst = os.path.join(args.outdir, os.path.basename(cand))
        shutil.copy(cand, dst)
        sys.stderr.write('[ciderseq_eccdna] %s -> %s\n' % (cand, dst))

    shutil.rmtree(staging, ignore_errors=True)


if __name__ == '__main__':
    main()
