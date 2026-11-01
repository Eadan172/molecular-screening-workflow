"""命中列表之后的 MoleCode 命令行。不注册到筛选工作流。"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from src.leadopt.structure import StructureError, export_hit_list, mermaid_to_validated, smiles_to_mermaid


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="leadopt.py",
        description="命中之后把 SMILES 转成 MoleCode 图，或把改过的图校验回分子。不参与筛选、对接或模型预测。",
    )
    sub = parser.add_subparsers(dest="command", required=True)

    to_graph = sub.add_parser("smiles-to-mermaid", help="SMILES 转成 MoleCode Mermaid")
    to_graph.add_argument("smiles", help="一个 SMILES")
    to_graph.add_argument("--name", default="Molecule", help="图里的分子名，默认 Molecule")
    to_graph.add_argument("--output", default=None, help="写入这个 .mmd 文件；省略则打印到标准输出")

    to_smiles = sub.add_parser("mermaid-to-smiles", help="MoleCode 图校验并转回标准 SMILES")
    to_smiles.add_argument("--input", required=True, help="Mermaid 文件")
    to_smiles.add_argument(
        "--require-smarts",
        default=None,
        help="可选。编辑后的结构必须含有这个子结构",
    )

    export = sub.add_parser("export-hits", help="把命中表逐条写成 .mmd")
    export.add_argument("hits", help="含 SMILES 列的 CSV，或每行一个 SMILES 的 .smi")
    export.add_argument("--out", required=True, help="输出目录，不会改写命中表")

    args = parser.parse_args(argv)
    try:
        if args.command == "smiles-to-mermaid":
            return _write_graph(args.smiles, args.name, args.output)
        if args.command == "mermaid-to-smiles":
            return _write_smiles(args.input, args.require_smarts)
        return _export(args.hits, args.out)
    except StructureError as exc:
        print(str(exc), file=sys.stderr)
        return 1


def _write_graph(smiles: str, name: str, output: str | None) -> int:
    graph = smiles_to_mermaid(smiles, name=name)
    if output:
        destination = Path(output)
        if destination.parent != Path("."):
            destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(graph.rstrip() + "\n", encoding="utf-8")
        print(output)
        return 0
    print(graph.rstrip())
    return 0


def _write_smiles(source: str, require_smarts: str | None) -> int:
    path = Path(source)
    if not path.is_file():
        raise StructureError(f"找不到文件: {source}")
    validated = mermaid_to_validated(path.read_text(encoding="utf-8"), require_smarts=require_smarts)
    print(validated.smiles)
    return 0


def _export(hits: str, out: str) -> int:
    summary = export_hit_list(hits, out)
    print(
        f"写出 {summary['ok']} 个结构，拒绝 {summary['rejected']} 个，"
        f"重复 {summary['duplicate']} 个。索引: {summary['index']}"
    )
    if summary["ok"] == 0:
        return 1
    return 0
