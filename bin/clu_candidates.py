#!/usr/bin/env python3
"""
clu_candidates.py - Extract eccDNA candidate clusters from RepeatExplorer2 (seqclust) output.

Extracted from ECCsplorer/lib/eccClusterer.py:eccClusterer.

Input:
    --clu_dir    : seqclust output root (directory that contains seqclust/clustering/)
    --pre_a      : readID prefix dataset A (default: TR)
    --pre_b      : readID prefix dataset B (default: CO)
Output:
    --out_candidates : cluster_candidates.csv (clusters enriched in dataset A)
    --out_summary    : comparative_cluster_table.csv (all clusters with scores)
Options:
    --ecc-prop   : minimum TR proportion to call a cluster eccDNA candidate (default: 0.8)

Notes (fixed vs the image copy):
    * seqclust writes COMPARATIVE_ANALYSIS_COUNTS.csv via R `write.table(col.names=TRUE)`.
    Its header row starts with "supercluster" and the sample-count columns are the
    SQLite column names wrapped in brackets, e.g. `[TR]` / `[CO]`. The previous parser
    only looked for a row whose first cell == 'cluster' and exact-matched 'TR'/'CO',
    which never matched the real file. The parser here normalizes brackets/quotes and
    locates the cluster-id column by name, then joins counts to CLUSTER_TABLE rows by
    cluster id (not by row order).
    * seqclust may die while generating its final HTML report (Rserve bug:
    "Error in index_html(imagemap) : object 'imagemap' not found"). In that case
    CLUSTER_TABLE.csv is never written, but COMPARATIVE_ANALYSIS_COUNTS.csv already
    is. We therefore build the cluster list from the comparative table when
    CLUSTER_TABLE.csv is missing, instead of failing.
"""

import argparse
import csv
import os
import sys

import numpy as np


def _norm(x):
    """Normalize a TSV cell: strip whitespace, surrounding quotes and SQLite []."""
    return x.strip().strip('"').strip('[]')


def read_cluster_rows(path):
    """Parse CLUSTER_TABLE.csv with python csv module (handles quotes/newlines),
    skipping the leading statistics rows and the header."""
    rows = []
    with open(path, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or not row[0].strip():
                continue
            first = row[0].strip().strip('"')
            if first.startswith('Number_of'):
                continue
            if first == 'Cluster':
                continue
            rows.append(row)
    return rows


def main():
    parser = argparse.ArgumentParser(description='Extract eccDNA candidate clusters from seqclust output')
    parser.add_argument('--clu_dir', required=True, help='seqclust output root directory')
    parser.add_argument('--pre_a', default='TR', help='ReadID prefix dataset A (default: TR)')
    parser.add_argument('--pre_b', default='CO', help='ReadID prefix dataset B (default: CO)')
    parser.add_argument('--out_candidates', required=True, help='Output cluster candidates CSV')
    parser.add_argument('--out_summary', required=True, help='Output comparative cluster table CSV')
    parser.add_argument('--ecc-prop', type=float, default=0.8, help='Min TR proportion for candidate (default: 0.8)')
    args = parser.parse_args()

    clustering_dir = args.clu_dir
    cluster_tab = os.path.join(clustering_dir, 'CLUSTER_TABLE.csv')
    comp_tab = os.path.join(clustering_dir, 'COMPARATIVE_ANALYSIS_COUNTS.csv')

    # --- parse comparative table first (needed in both branches) ---
    if not os.path.isfile(comp_tab):
        print(f'ERROR: seqclust outputs not found under {clustering_dir}', file=sys.stderr)
        print(f'  CLUSTER_TABLE.csv exists: {os.path.isfile(cluster_tab)}', file=sys.stderr)
        print(f'  COMPARATIVE_ANALYSIS_COUNTS.csv exists: {os.path.isfile(comp_tab)}', file=sys.stderr)
        sys.exit(1)

    comp_header = None
    comp_rows = []
    with open(comp_tab, 'r', newline='') as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if not row or not row[0].strip():
                continue
            first = _norm(row[0])
            if first.startswith('#'):
                continue
            # Header row: contains a cluster-id column (cluster/clusterindex) and
            # at least one of the sample-prefix columns.
            if comp_header is None and any(_norm(c) in ('cluster', 'clusterindex') for c in row):
                comp_header = [_norm(c) for c in row]
                continue
            comp_rows.append(row)

    if comp_header is None:
        print(f'ERROR: no header row found in {comp_tab}', file=sys.stderr)
        sys.exit(1)

    # locate cluster-id column and TR/CO columns
    try:
        cl_id_col = comp_header.index('cluster')
    except ValueError:
        try:
            cl_id_col = comp_header.index('clusterindex')
        except ValueError:
            cl_id_col = 0
    tr_col = co_col = None
    for i, h in enumerate(comp_header):
        if h == args.pre_a and tr_col is None:
            tr_col = i
        elif h == args.pre_b and co_col is None:
            co_col = i
    if tr_col is None or co_col is None:
        print(f'ERROR: comparative columns for {args.pre_a}/{args.pre_b} not found '
                f'(header={comp_header})', file=sys.stderr)
        sys.exit(1)

    comp_cl = np.array(comp_rows) if comp_rows else np.zeros((0, len(comp_header)), dtype=str)

    # --- cluster list: prefer CLUSTER_TABLE.csv, fall back to comparative table ---
    if os.path.isfile(cluster_tab):
        rows = read_cluster_rows(cluster_tab)
        if not rows:
            print('No clusters found.', file=sys.stderr)
            open(args.out_candidates, 'w').close()
            open(args.out_summary, 'w').close()
            sys.exit(0)
        cl_tab = np.array(rows)
        cl_ids = cl_tab[:, 0].astype(int)
    else:
        # seqclust crashed in its final HTML-report step; cluster list is taken
        # from the comparative counts table (cluster id column only).
        print(f'WARN: {cluster_tab} not found; building cluster list from '
                f'{os.path.basename(comp_tab)}', file=sys.stderr)
        ids = []
        seen = set()
        for r in comp_cl:
            if len(r) <= cl_id_col:
                continue
            tok = r[cl_id_col].strip().strip('"').strip('[]')
            if not tok.lstrip('-').isdigit():
                continue
            cid = int(tok)
            if cid not in seen:
                seen.add(cid)
                ids.append(cid)
        if not ids:
            print('No clusters found.', file=sys.stderr)
            open(args.out_candidates, 'w').close()
            open(args.out_summary, 'w').close()
            sys.exit(0)
        cl_ids = np.array(sorted(ids))
        cl_tab = cl_ids.reshape(-1, 1)

    def add_cnt_col(tab, i):
        # join by cluster id, not by row order
        m = {int(r[cl_id_col]): r[i] for r in comp_cl if len(r) > max(cl_id_col, i)}
        col = np.array([m.get(int(cid), '0') for cid in cl_ids]).reshape(-1, 1)
        return np.hstack((tab, col))

    cl_tab = add_cnt_col(cl_tab, tr_col)
    cl_tab = add_cnt_col(cl_tab, co_col)

    # score + candidate filtering
    tr_idx = cl_tab.shape[1] - 2
    co_idx = cl_tab.shape[1] - 1
    tr_vals = cl_tab[:, tr_idx].astype(float)
    co_vals = cl_tab[:, co_idx].astype(float)
    score = np.round(tr_vals / (tr_vals + co_vals + 1e-9), 3)
    cl_tab = np.hstack((cl_tab, score.reshape(-1, 1)))

    np.savetxt(args.out_summary, cl_tab, delimiter='\t', fmt='%s')
    cand = cl_tab[score >= args.ecc_prop]
    np.savetxt(args.out_candidates, cand, delimiter='\t', fmt='%s')

    print(f'Clusters total: {cl_tab.shape[0]}, candidates (score>={args.ecc_prop}): {cand.shape[0]}', file=sys.stderr)


if __name__ == '__main__':
    main()
