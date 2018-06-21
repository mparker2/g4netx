import os
import re
from operator import methodcaller
import networkx as nx
import click
from fileutils import FastaReader, BedWriter, sort_bed_file

def run_combinations(g_runs, run_length, max_loop_length):
    g_runs = tuple(g_runs)
    for i, g1 in enumerate(g_runs):
        for g2 in g_runs[i + 1:]:
            dist = (g2 - (g1 + run_length))
            if dist <= 0:
                continue
            elif dist > max_loop_length:
                break
            else:
                yield g1, g2, dist


def build_g_run_graph(seq, run_length=3, max_loop_length=7, rev_comp=False):
    graphs = {'+': nx.Graph(), '-': nx.Graph()}
    for strand, base in [('+', 'G'), ('-', 'C')]:
        regex = '(?={base}{{{run_length}}})(?i)'.format(
            base=base, run_length=run_length)
        g_runs = map(methodcaller('start'), re.finditer(regex, seq))
        for g1, g2, dist in run_combinations(
                g_runs, run_length, max_loop_length):
            graphs[strand].add_edge(g1, g2, dist=dist)
    return graphs


def enumerate_subgraphs(g, k):
    for sg in nx.connected_components(g):
        if len(sg) >= k:
            for v in sg:
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


def subgraph_to_range(subgraph_iter, run_length, strand):
    for sg in subgraph_iter:
        yield [(start, start + run_length) for start in sg], strand


def g4netx(seq, run_length, max_loop_length=7):
    graphs = build_g_run_graph(seq, run_length, max_loop_length)
    for strand in ('+', '-'):
        subgraph_iter = enumerate_subgraphs(graphs[strand], 4)
        yield from subgraph_to_range(subgraph_iter, run_length, strand)


def as_bed12(ranges_iter, ref_name,
             bed_name='g4netx',
             colour='85,118,209'):
    bed12 = '\t'.join(('{}',) * 12)
    for g4_ranges, strand in ranges_iter:
        start = g4_ranges[0][0]
        end = g4_ranges[-1][1]
        count = len(g4_ranges)
        lengths = [str(e - s) for s, e in g4_ranges]
        starts = [str(s - start) for s, _ in g4_ranges]
        yield bed12.format(
            ref_name, start, end, bed_name,
            0, strand, start, start,
            colour, count,
            ','.join(lengths), ','.join(starts)
        )

        
@click.command()
@click.option('--fasta', '-f', required=False, default='-')
@click.option('--bed', '-b', required=False, default='-')
@click.option('--g-run-length', '-g', required=False, default=3)
@click.option('--max-loop-length', '-l', required=False, default=7)
def g4netx_ci(fasta, bed, g_run_length, max_loop_length):
    with FastaReader(fasta) as fa, BedWriter() as bed_unsorted:
        for ref_name, seq in fa.parse_fasta():
            ranges = g4netx(seq, g_run_length, max_loop_length)
            for record in as_bed12(ranges, ref_name):
                bed_unsorted.write(record)
    sorted_bed_iter = sort_bed_file(bed_unsorted.fn)
    with BedWriter(bed) as bed_out:
        for record in sorted_bed_iter:
            bed_out.write(record)
    os.remove(bed_unsorted.fn)

if __name__ == '__main__':
    g4netx_ci()