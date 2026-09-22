"""Экспорт сохранённых графов (.npz) в читаемые форматы.

Примеры:
    # один граф -> edge list
    python export_graph.py results_v2/graphs/pl_dmin2_g2.5_n1000_d2_r0.npz \
        --format edgelist --out graph.txt

    # первые 20 графов из папки -> CSV (source,target) в export/
    python export_graph.py results_v2/graphs --format csv --limit 20 --out export/

    # весь planted-ансамбль -> GraphML (с метками сообществ)
    python export_graph.py results_v2/graphs --format graphml \
        --glob "planted_*" --out export_graphml/

Форматы: edgelist, csv, graphml, gml, json.
Для planted-графов метки сообществ сохраняются как атрибут узла `community`.
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import os
import sys

import networkx as nx

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from src.generators import load_graph

FORMATS = ("edgelist", "csv", "graphml", "gml", "json")
EXT = {"edgelist": "txt", "csv": "csv", "graphml": "graphml", "gml": "gml", "json": "json"}


def write_graph(g: nx.Graph, labels, path: str, fmt: str) -> None:
    if labels is not None:
        for v, lab in zip(g.nodes(), labels):
            g.nodes[v]["community"] = int(lab)

    if fmt == "edgelist":
        with open(path, "w", encoding="utf-8") as fh:
            for u, v in g.edges():
                fh.write(f"{u} {v}\n")
    elif fmt == "csv":
        with open(path, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh)
            w.writerow(["source", "target"])
            for u, v in g.edges():
                w.writerow([u, v])
    elif fmt == "graphml":
        nx.write_graphml(g, path)
    elif fmt == "gml":
        nx.write_gml(g, path)
    elif fmt == "json":
        data = {
            "n": g.number_of_nodes(),
            "m": g.number_of_edges(),
            "edges": [[int(u), int(v)] for u, v in g.edges()],
        }
        if labels is not None:
            data["labels"] = [int(x) for x in labels]
        with open(path, "w", encoding="utf-8") as fh:
            json.dump(data, fh)
    else:
        raise ValueError(f"unknown format: {fmt}")


def collect_files(inp: str, pattern: str | None, limit: int | None):
    if os.path.isdir(inp):
        pat = os.path.join(inp, (pattern or "*") + ".npz")
        files = sorted(glob.glob(pat))
    else:
        files = [inp]
    if limit is not None:
        files = files[:limit]
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input", help=".npz файл или папка с .npz")
    ap.add_argument("--format", choices=FORMATS, default="edgelist")
    ap.add_argument("--out", default=None,
                    help="файл (для одного графа) или папка (для набора)")
    ap.add_argument("--glob", default=None, help="маска файлов, напр. 'planted_*'")
    ap.add_argument("--limit", type=int, default=None, help="максимум графов")
    args = ap.parse_args(argv)

    files = collect_files(args.input, args.glob, args.limit)
    if not files:
        print("[export] не найдено .npz файлов")
        return 1

    single = len(files) == 1 and not os.path.isdir(args.input)
    if single and args.out:
        out_path = args.out
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        g, labels = load_graph(files[0])
        write_graph(g, labels, out_path, args.format)
        print(f"[export] {files[0]} -> {out_path} ({g.number_of_nodes()} узлов, "
              f"{g.number_of_edges()} рёбер)")
        return 0

    out_dir = args.out or os.path.join(os.path.dirname(files[0]) or ".", "export")
    os.makedirs(out_dir, exist_ok=True)
    for i, f in enumerate(files, 1):
        g, labels = load_graph(f)
        stem = os.path.splitext(os.path.basename(f))[0]
        dst = os.path.join(out_dir, f"{stem}.{EXT[args.format]}")
        write_graph(g, labels, dst, args.format)
        if i % 200 == 0 or i == len(files):
            print(f"[export] {i}/{len(files)} -> {out_dir}")
    print(f"[export] готово: {len(files)} графов в {out_dir} ({args.format})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
