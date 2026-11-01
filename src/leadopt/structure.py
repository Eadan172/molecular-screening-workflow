"""SMILES 与 MoleCode Mermaid 的往返，并用 RDKit 拒绝不合法结构。

只在筛选命中之后做查看和局部修改。不要从枚举、QSAR、ADMET 或对接里调用。
"""

from __future__ import annotations

import csv
import io
from contextlib import contextmanager, redirect_stderr, redirect_stdout
from dataclasses import dataclass
from pathlib import Path


class StructureError(ValueError):
    """SMILES 或 MoleCode 图没有通过 RDKit 校验。"""


@dataclass
class ValidatedStructure:
    """一次通过校验的结构。mol 是已经 Sanitize 的 RDKit 分子。"""

    smiles: str
    mol: object
    mermaid: str


_SMILES_PRIORITY = (
    "canonical_smiles",
    "标准smiles",
    "smiles",
    "smi",
    "input_smiles",
    "输入smiles",
)


def smiles_to_mermaid(smiles: str, name: str = "Molecule") -> str:
    """把一个 SMILES 写成 MoleCode Mermaid。无法解析时抛出 StructureError。"""
    chem, mol_to_mermaid, _parse = _libraries()
    text = (smiles or "").strip()
    if not text:
        raise StructureError("SMILES 是空的")
    with _quiet_rdkit():
        mol = chem.MolFromSmiles(text)
    if mol is None:
        raise StructureError(f"SMILES 无法解析: {_short(text)}")
    return mol_to_mermaid(mol, name=name or "Molecule")


def mermaid_to_validated(
    mermaid_text: str,
    require_smarts: str | None = None,
) -> ValidatedStructure:
    """把 MoleCode 图转成通过 RDKit 校验的分子和标准 SMILES。

    strict 模式下，无法识别的原子会被拒绝。价态不合法、解析失败，
    或结构里没有 require_smarts 描述的子结构时，同样拒绝。
    """
    chem, _to_mermaid, mermaid_to_mol = _libraries()
    text = _strip_fence(mermaid_text)
    if not text:
        raise StructureError("MoleCode 图是空的")

    captured = io.StringIO()
    try:
        with _quiet_rdkit(), redirect_stdout(captured), redirect_stderr(captured):
            mol = mermaid_to_mol(text, strict=True)
    except Exception as exc:
        raise StructureError(f"MoleCode 图无法解析: {exc}") from exc
    if mol is None:
        raise StructureError("MoleCode 图无法解析，或含有 RDKit 不接受的原子")

    try:
        with _quiet_rdkit():
            chem.SanitizeMol(mol)
            smiles = chem.MolToSmiles(mol)
    except Exception as exc:
        raise StructureError(f"RDKit 校验失败: {exc}") from exc
    if not smiles or mol.GetNumAtoms() == 0:
        raise StructureError("RDKit 校验失败: 没有得到分子")

    pattern = (require_smarts or "").strip()
    if pattern:
        with _quiet_rdkit():
            query = chem.MolFromSmarts(pattern)
        if query is None:
            raise StructureError(f"子结构 SMARTS 无法解析: {_short(pattern)}")
        if not mol.HasSubstructMatch(query):
            raise StructureError(f"编辑后的结构不含要求的子结构: {_short(pattern)}")

    return ValidatedStructure(smiles=smiles, mol=mol, mermaid=text)


def export_hit_list(hits_path: Path | str, out_dir: Path | str) -> dict:
    """把命中表里的 SMILES 写成各自的 .mmd，并留下 index.csv。

    不修改命中表本身。无法解析的行记为 rejected，重复 SMILES 记为 duplicate。
    再次导出到同一目录时，只替换上一次写下的 hit_*.mmd 和 index.csv。
    """
    source = Path(hits_path)
    if not source.is_file():
        raise StructureError(f"找不到命中表: {hits_path}")
    destination = Path(out_dir)
    if destination.exists() and not destination.is_dir():
        raise StructureError(f"输出位置不是目录: {out_dir}")
    rows = _read_hit_smiles(source)
    if not rows:
        raise StructureError("命中表里没有 SMILES")

    destination.mkdir(parents=True, exist_ok=True)
    for stale in destination.glob("hit_*.mmd"):
        stale.unlink()
    index_path = destination / "index.csv"
    if index_path.exists():
        index_path.unlink()

    written: list[dict] = []
    seen: set[str] = set()
    ok = 0
    for row_number, raw in rows:
        smiles = raw.strip()
        record = {"row": row_number, "file": "", "smiles": smiles, "status": "", "reason": ""}
        if smiles in seen:
            record["status"] = "duplicate"
            record["reason"] = "与更早的一行相同"
            written.append(record)
            continue
        seen.add(smiles)
        try:
            graph = smiles_to_mermaid(smiles, name=f"Hit{ok + 1:04d}")
        except StructureError as exc:
            record["status"] = "rejected"
            record["reason"] = str(exc)
            written.append(record)
            continue
        filename = f"hit_{ok + 1:04d}.mmd"
        (destination / filename).write_text(graph.rstrip() + "\n", encoding="utf-8")
        record["file"] = filename
        record["status"] = "ok"
        written.append(record)
        ok += 1

    with index_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["row", "file", "smiles", "status", "reason"])
        writer.writeheader()
        writer.writerows(written)

    rejected = sum(1 for item in written if item["status"] == "rejected")
    duplicate = sum(1 for item in written if item["status"] == "duplicate")
    return {
        "ok": ok,
        "rejected": rejected,
        "duplicate": duplicate,
        "index": index_path.name,
        "rows": written,
    }


def _libraries():
    try:
        from rdkit import Chem
        from molecode import mermaid_to_mol, mol_to_mermaid
    except ImportError as exc:
        raise StructureError(
            "缺少可选依赖 molecode 或 RDKit。"
            "请执行: python -m pip install -r requirements-leadopt.txt"
        ) from exc
    return Chem, mol_to_mermaid, mermaid_to_mol


def _read_hit_smiles(path: Path) -> list[tuple[int, str]]:
    if path.suffix.lower() in {".smi", ".smiles"}:
        rows = []
        for number, line in enumerate(path.read_text(encoding="utf-8-sig").splitlines(), start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            rows.append((number, stripped.split()[0]))
        return rows

    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames:
            raise StructureError("命中表没有表头")
        column = _smiles_column(reader.fieldnames)
        if column is None:
            raise StructureError("命中表里没有 SMILES 列")
        rows = []
        for number, row in enumerate(reader, start=2):
            value = (row.get(column) or "").strip()
            if value:
                rows.append((number, value))
        return rows


def _smiles_column(fieldnames) -> str | None:
    lookup: dict[str, str] = {}
    for name in fieldnames:
        if name is None:
            continue
        lookup.setdefault(str(name).strip().lower(), name)
    for candidate in _SMILES_PRIORITY:
        if candidate in lookup:
            return lookup[candidate]
    return None


def _strip_fence(text: str) -> str:
    body = (text or "").strip()
    if not body.startswith("```"):
        return body
    lines = body.splitlines()
    if lines and lines[0].startswith("```"):
        lines = lines[1:]
    if lines and lines[-1].startswith("```"):
        lines = lines[:-1]
    return "\n".join(lines).strip()


def _short(value: str, limit: int = 120) -> str:
    compact = " ".join(value.split())
    if len(compact) <= limit:
        return compact
    return compact[:limit] + "..."


@contextmanager
def _quiet_rdkit():
    try:
        from rdkit import RDLogger
    except ImportError:
        yield
        return
    RDLogger.DisableLog("rdApp.*")
    try:
        yield
    finally:
        RDLogger.EnableLog("rdApp.*")
