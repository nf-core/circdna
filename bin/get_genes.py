#!/usr/bin/env python3

from ac_util import *

# This file is imported by amplicon_classifier.py to get genes

def merge_intervals(feature_dict, tol=1):
    for item, usort_intd in feature_dict.items():
        for chrom, usort_ints in usort_intd.items():
            # sort ints
            sort_ints = sorted(usort_ints)
            # merge sorted ints
            mi = [sort_ints[0]]
            for ival in sort_ints[1:]:
                if ival[0] <= mi[-1][1] + tol:
                    ui = (mi[-1][0], ival[1])
                    mi[-1] = ui

                else:
                    mi.append(ival)

            feature_dict[item][chrom] = mi


def get_gene_ends(gs, ge, strand, a, b):
    has5p, has3p = False, False
    if strand == "-":
        gs, ge = ge, gs

    if a <= gs <= b:
        has5p = True

    if a <= ge < b:
        has3p = True

    return has5p, has3p


def get_genes_from_intervals(gene_lookup, feature_dict, gseg_cn_d):
    feat_to_gene_trunc = defaultdict(lambda: defaultdict(set))
    feat_to_gene_cn = defaultdict(lambda: defaultdict(float))
    for feat_name, curr_fd in feature_dict.items():
        for chrom, intlist in curr_fd.items():
            for a, b in intlist:
                ogenes_ints = gene_lookup[chrom][a:b]
                for gp in ogenes_ints:
                    gname, strand = gp.data
                    has5p, has3p = get_gene_ends(gp.begin, gp.end, strand, a, b)
                    gsegs_hit = gseg_cn_d[chrom][gp.begin:gp.end]
                    gene_valid_cns = [x.data for x in gsegs_hit if x.end - x.begin > 1000]
                    if gene_valid_cns:
                        gene_cn = max(gene_valid_cns)
                    else:
                        gene_cn = "unknown"
                    feat_to_gene_cn[feat_name][gname] = gene_cn
                    if has5p:
                        feat_to_gene_trunc[feat_name][gname].add("5p")
                    if has3p:
                        feat_to_gene_trunc[feat_name][gname].add("3p")

    return feat_to_gene_trunc, feat_to_gene_cn


def get_gseg_cns(graphf, add_chr_tag):
    gseg_cn_d = defaultdict(IntervalTree)
    with open(graphf) as infile:
        for line in infile:
            fields = line.rsplit()
            if line.startswith("sequence"):
                c, s, e = fields[1].rsplit(":")[0], int(fields[1].rsplit(":")[1][:-1]), int(fields[2].rsplit(":")[1][:-1])+1
                cn = float(fields[3])
                if add_chr_tag and not c.startswith('chr'):
                    c = "chr" + c

                gseg_cn_d[c].addi(s, e, cn)

    return gseg_cn_d


def extract_gene_list(sname, ampN, gene_lookup, cycleList, segSeqD, bfb_cycle_inds, ecIndexClusters,
                      invalidInds, bfbStat, ecStat, ampClass, graphf, add_chr_tag, prefix):
    feature_dict = {}
    gseg_cn_d = get_gseg_cns(graphf, add_chr_tag)
    invalidSet = set(invalidInds)
    all_used = invalidSet.union(bfb_cycle_inds)
    used_segs = defaultdict(IntervalTree)
    graph_cns = get_graph_cns(graphf, add_chr_tag)
    if bfbStat:
        # collect unmerged genomic intervals comprising the feature
        bfb_interval_dict = defaultdict(list)
        for b_ind in bfb_cycle_inds:
            if b_ind not in invalidSet:
                for c_id in cycleList[b_ind]:
                    chrom, l, r = segSeqD[abs(c_id)]
                    if chrom:
                        bfb_interval_dict[chrom].append((l, r))
                        used_segs[chrom].addi(l, r+1)

        feature_dict["BFB_1"] = bfb_interval_dict

    if ecStat:
        # collect unmerged genomic intervals comprising the feature
        for amp_ind, ec_cycle_inds in enumerate(ecIndexClusters):
            ec_interval_dict = defaultdict(list)
            for e_ind in ec_cycle_inds:
                all_used.add(e_ind)
                if e_ind not in invalidSet:
                    for c_id in cycleList[e_ind]:
                        chrom, l, r = segSeqD[abs(c_id)]
                        if chrom:
                            used_segs[chrom].addi(l, r+1)
                            # chop out low cn regions
                            seg_t = IntervalTree([Interval(l, r+1)])
                            olapping_low_cns = [x for x in graph_cns[chrom][l:r+1] if x.data < 4]
                            for x in olapping_low_cns:
                                seg_t.chop(x.begin, x.end+1)

                            for x in seg_t:
                                ec_interval_dict[chrom].append((x.begin, x.end))

            feature_dict["ecDNA_" + str(amp_ind + 1)] = ec_interval_dict

    if ampClass != "No amp/Invalid":
        other_interval_dict = defaultdict(list)
        for o_ind in range(len(cycleList)):
            if o_ind not in all_used:
                for c_id in cycleList[o_ind]:
                    if abs(c_id) not in used_segs:
                        chrom, l, r = segSeqD[abs(c_id)]
                        if not used_segs[chrom][l:r] and not chrom is None:
                            other_interval_dict[chrom].append((l, r))

        if not ecStat and not bfbStat:
            feature_dict[ampClass + "_1"] = other_interval_dict
        else:
            feature_dict["unknown_1"] = other_interval_dict

    # merge all the intervals in each list of intervals
    # tot_init_intervals = sum([len(ilist) for fd in feature_dict.values() for ilist in fd.values()])
    merge_intervals(feature_dict)
    # tot_final_intervals = sum([len(ilist) for fd in feature_dict.values() for ilist in fd.values()])
    # print("Feature extraction: started with " + str(tot_init_intervals) + " unmerged intervals, finished with " + str(
    #     tot_final_intervals) + " intervals")

    write_interval_beds(prefix, sname, ampN, feature_dict)
    return get_genes_from_intervals(gene_lookup, feature_dict, gseg_cn_d)
