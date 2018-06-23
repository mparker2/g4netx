import os
import regex as re
from operator import methodcaller
import networkx as nx
import click
from .fileutils import FastaReader, BedWriter, sort_bed_file
from . import g4filter as filt

def get_cluster(g_runs, max_loop_length, allow_bulges=False):
    try:
        match = next(g_runs)
        g = match.span(0)
        is_bulged = test_if_bulged(match, allow_bulges)
    except StopIteration:
        # no G runs in sequence
        return
    cluster_end = g[1]
    cluster = [(*g, is_bulged),]
    while True:
        try:
            match = next(g_runs)
            g = match.span(0)
            is_bulged = test_if_bulged(match, allow_bulges)
        except StopIteration:
            if cluster:
                yield cluster
            break
        if g[0] - cluster_end <= max_loop_length:
            cluster_end = g[1]
            cluster.append((*g, is_bulged))
        else:
            yield cluster
            cluster_end = g[1]
            cluster = [(*g, is_bulged),]


def test_if_bulged(match, allow_bulges):
    if not allow_bulges:
        return False
    else:
        return any(match.groups())


def run_combinations(g_runs, max_loop_length):
    for i, g1 in enumerate(g_runs):
        for g2 in g_runs[i + 1:]:
            dist = g2[0] - g1[1]
            if dist <= 0:
                continue
            elif dist > max_loop_length:
                break
            else:
                yield g1, g2, dist


def create_bulged_regex(base, g_run_length, max_bulge_length):
    if base == 'G':
        bulge = '[ACT]{1,' + str(max_bulge_length) + '}'
    else:
        bulge = '[AGT]{1,' + str(max_bulge_length) + '}'
    regex = base
    bulge_backrefs = []
    for i in range(1, g_run_length):
        regex += '('
        for j in bulge_backrefs:
            regex += '(?(' + str(j) + ')|'
        regex += bulge
        regex += ')' * len(bulge_backrefs)
        regex += ')?'
        bulge_backrefs.append(i)
        regex += base
    return regex


def build_g_run_graphs(seq, run_length=3, max_loop_length=7,
                       allow_bulges=False, max_bulge_length=5,
                       soft_mask=False):
    for strand, base in [('+', 'G'), ('-', 'C')]:
        if allow_bulges:
            if run_length == 2:
                # no bulges for 2 tetrad
                regex = '{base}{{{run_length}}}'.format(
                    base=base, run_length=run_length)
            else:
                regex = create_bulged_regex(
                    base, run_length, max_bulge_length)
        else:
            regex = '{base}{{{run_length}}}'.format(
                base=base, run_length=run_length)
        if soft_mask:
            re_args = (regex.IGNORECASE)
        else:
            re_args = ()
        g_runs = re.finditer(regex, seq, *re_args, overlapped=True)
        clusterer = get_cluster(g_runs, run_length, max_loop_length)
        while True:
            graph = nx.Graph()
            try:
                g_run_cluster = next(clusterer)
            except StopIteration:
                break
            for g1, g2, dist in run_combinations(
                    g_run_cluster, max_loop_length):
                graph.add_edge(g1, g2, dist=dist)
            yield graph, strand


def enumerate_subgraphs(g, k):
    if len(g) >= k:
        for v in g:
            v_ext = [u for u in g.neighbors(v) if v < u]
            yield from extend_subgraph(g, [v], v_ext, k, 0)


def extend_subgraph(g, v_sg, v_ext, k, i):
    if len(v_sg) == k:
        yield v_sg
        return
    j = 0
    while len(v_ext) != 0:
        w = v_ext.pop(0)
        v_ext_ = [u for u in g.neighbors(w)
                 if u not in v_sg and u > w]
        yield from extend_subgraph(g, v_sg + [w, ], v_ext_, k, i + 1)


def score_g4(tetrads, n_bulges, length):
    not_bulged = 4 - n_bulges
    tetrad_score = tetrads ** 2
    score = tetrad_score * not_bulged
    score += (tetrad_score * n_bulges) / 2
    return score / (length - 4 * tetrads)


def g4netx(seq, run_length, max_loop_length=7,
           max_bulges=False, max_bulge_length=5,
           soft_mask=False):
    allow_bulges = bool(max_bulges)
    for g, strand in build_g_run_graphs(seq, run_length,
                                        max_loop_length,
                                        allow_bulges,
                                        max_bulge_length,
                                        soft_mask):
        for sg in enumerate_subgraphs(g, 4):
            n_bulged = sum(tet[2] for tet in sg)
            if n_bulged <= max_bulges:
                score = score_g4(run_length, n_bulged, sg[-1][1] - sg[0][0])
                yield sg, score, strand


def as_bed12(ranges_iter, ref_name,
             bed_name='g4netx',
             colour='85,118,209'):
    bed12 = '\t'.join(('{}',) * 12)
    for g4_ranges, score, strand in ranges_iter:
        start = g4_ranges[0][0]
        end = g4_ranges[-1][1]
        count = len(g4_ranges)
        lengths = [str(e - s) for s, e, _ in g4_ranges]
        starts = [str(s - start) for s, *_ in g4_ranges]
        yield bed12.format(
            ref_name, start, end, bed_name,
            '{:.2f}'.format(score),
            strand, start, start,
            colour, count,
            ','.join(lengths),
            ','.join(starts)
        )

        
@click.command()
@click.option('--fasta', '-f', required=False, default='-')
@click.option('--bed', '-b', required=False, default='-')
@click.option('--min-g-run-length', '-i', required=False,
              default=3, type=click.IntRange(2, 6))
@click.option('--max-g-run-length', '-j', required=False,
              default=6, type=click.IntRange(2, 6))
@click.option('--max-loop-length', '-l', required=False,
              default=7, type=click.IntRange(1, 100))
@click.option('--max-bulges', '-u', required=False, default=0, type=click.IntRange(0, 4))
@click.option('--max-bulge-length', '-v', required=False, default=5,
              type=click.IntRange(1, 7))
@click.option('--filter-overlapping', '-F',
              is_flag=True, required=False, default=False)
@click.option('--merge-overlapping', '-M',
              is_flag=True, required=False, default=False)
@click.option('--soft-mask', '-s',
              is_flag=True, required=False, default=False)
def g4netx_cli(fasta, bed,
               min_g_run_length,
               max_g_run_length,
               max_loop_length,
               max_bulges,
               max_bulge_length,
               filter_overlapping,
               merge_overlapping,
               soft_mask):
    g_run_lengths = list(range(min_g_run_length, max_g_run_length + 1))
    with FastaReader(fasta) as fa, BedWriter() as bed_unsorted:
        for ref_name, seq in fa.parse_fasta():
            for run_len in g_run_lengths:
                ranges = g4netx(seq, run_len,
                                max_loop_length,
                                max_bulges,
                                max_bulge_length,
                                soft_mask)
                for record in as_bed12(ranges, ref_name):
                    bed_unsorted.write(record)
    sorted_bed_iter = sort_bed_file(bed_unsorted.fn)
    if filter_overlapping:
        f = filt.filter_overlapping
    elif merge_overlapping:
        f = filt.merge_overlapping
    else:
        f = None
    with BedWriter(bed) as bed_out:
        for record in filt.apply_filter_method(sorted_bed_iter, f):
            bed_out.write(record)
    os.remove(bed_unsorted.fn)

if __name__ == '__main__':
    g4netx_cli()