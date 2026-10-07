"""用 RDKit 计算理化性质，并按用户写下的数值要求标注结果。"""

from __future__ import annotations

from src.agent.models import Constraint, RequestSpec


COMPUTED = {
    "MW",
    "LogP",
    "TPSA",
    "HBD",
    "HBA",
    "RotBonds",
    "QED",
    "Fsp3",
    "AromaticRings",
    "Rings",
    "HeavyAtoms",
    "SA",
    "IC50",
    "Ki",
    "Kd",
    "Affinity",
}

_EXTERNAL_COLUMNS = {
    "ic50": "IC50",
    "pred_ic50": "IC50",
    "predicted_ic50": "IC50",
    "ki": "Ki",
    "kd": "Kd",
    "binding_affinity": "Affinity",
    "affinity": "Affinity",
    "docking_score": "Affinity",
    "vina_score": "Affinity",
}

_ALERTS = [
    ("硝基", "[N+](=O)[O-]"),
    ("酰卤", "[CX3](=[OX1])[F,Cl,Br,I]"),
    ("迈克尔受体", "C=CC=O"),
    ("过氧化物", "[OX2][OX2]"),
]


def calculate_molecules(spec: RequestSpec) -> dict:
    rdkit = _load_rdkit()
    summary = {
        "n_input": len(spec.smiles),
        "n_valid": 0,
        "n_passed": 0,
        "n_failed": 0,
        "n_incomplete": 0,
        "n_invalid": 0,
        "rdkit": rdkit is not None,
        "real_sa": bool(rdkit and rdkit["sascorer"] is not None),
        "not_applied": [],
        "text_requirements": [],
        "fail_counts": {},
        "notes": [],
    }
    if not spec.smiles:
        summary["notes"].append("输入里没有 SMILES，已跳过分子性质计算。")
        return {"molecules": [], "summary": summary}

    if rdkit is None:
        summary["notes"].append("当前环境没有 RDKit，无法计算理化性质。请先运行 setup.sh 或 setup.bat。")
        molecules = [
            _empty_row(smiles, valid=False, reason="未安装 RDKit")
            for smiles in spec.smiles
        ]
        return {"molecules": molecules, "summary": summary}

    molecules = [_properties(smiles, rdkit) for smiles in spec.smiles]
    summary["n_valid"] = sum(1 for row in molecules if row["valid"])
    summary["n_invalid"] = sum(1 for row in molecules if not row["valid"])
    if not summary["real_sa"]:
        summary["notes"].append(
            "当前 RDKit 没有 SA_Score 片段库，表中的 SA代理 只表示结构复杂度，不用于 SA 阈值过滤。"
        )
    return {"molecules": molecules, "summary": summary}


def merge_external(molecules: list[dict], tables: list[dict]) -> int:
    """把外部软件或接口返回的 IC50 / 结合能按 SMILES 对齐。返回对齐上的行数。"""
    matched = 0
    for table_row in tables:
        smiles = _row_smiles(table_row)
        if not smiles:
            continue
        target = _find_molecule(molecules, smiles)
        if target is None:
            continue
        updated = False
        for key, value in table_row.items():
            canon = _EXTERNAL_COLUMNS.get(str(key).strip().lower())
            number = _number(value)
            if canon and number is not None:
                previous = target.get(canon)
                # IC50 和结合能都是越低越有利，同一分子多条记录时保留更低值。
                if previous is None or number < previous:
                    target[canon] = number
                updated = True
        if updated:
            matched += 1
    return matched


def annotate(spec: RequestSpec, molecules: list[dict], summary: dict) -> None:
    real_sa = bool(summary.get("real_sa"))
    fail_counts: dict[str, int] = {}
    not_applied = []
    text_requirements = []
    _collect_global_requirements(spec.constraints, {"not_applied": not_applied, "text_requirements": text_requirements}, real_sa)

    for row in molecules:
        failed = []
        pending = []
        if not row.get("valid"):
            row["numeric_status"] = "invalid"
            row["failed"] = []
            row["pending"] = []
            continue
        for constraint in spec.constraints:
            if constraint.kind != "numeric":
                continue
            if constraint.name == "SA" and not real_sa:
                continue
            if constraint.name not in COMPUTED:
                continue
            label = constraint.display()
            if constraint.op == "between":
                if constraint.value is None or constraint.high is None:
                    pending.append(label)
                    continue
            elif constraint.value is None:
                pending.append(label)
                continue
            value = row.get(constraint.name)
            if value is None:
                pending.append(label)
                continue
            if not _passes(value, constraint):
                failed.append(label)
                fail_counts[label] = fail_counts.get(label, 0) + 1
        row["failed"] = failed
        row["pending"] = pending
        if failed:
            row["numeric_status"] = "fail"
        elif pending:
            row["numeric_status"] = "incomplete"
        else:
            row["numeric_status"] = "pass"

    summary["not_applied"] = not_applied
    summary["text_requirements"] = text_requirements
    summary["fail_counts"] = fail_counts
    summary["n_passed"] = sum(1 for row in molecules if row.get("numeric_status") == "pass")
    summary["n_failed"] = sum(1 for row in molecules if row.get("numeric_status") == "fail")
    summary["n_incomplete"] = sum(1 for row in molecules if row.get("numeric_status") == "incomplete")
    summary["n_invalid"] = sum(1 for row in molecules if row.get("numeric_status") == "invalid")
    summary["n_valid"] = sum(1 for row in molecules if row.get("valid"))
    molecules.sort(
        key=lambda row: (
            {"pass": 0, "incomplete": 1, "fail": 2, "invalid": 3}.get(row.get("numeric_status"), 9),
            -(row.get("QED") or -1),
        )
    )


def _collect_global_requirements(constraints: list[Constraint], summary: dict, real_sa: bool) -> None:
    not_applied = summary.setdefault("not_applied", [])
    text_requirements = summary.setdefault("text_requirements", [])
    for constraint in constraints:
        if constraint.kind != "numeric":
            text_requirements.append(constraint.display())
            continue
        if constraint.name == "SA" and not real_sa:
            not_applied.append(constraint.display() + "（没有真实 SA Score，未过滤）")
        elif constraint.name not in COMPUTED:
            not_applied.append(constraint.display() + "（程序不能直接计算这项）")


def _properties(smiles: str, rdkit: dict) -> dict:
    chem = rdkit["Chem"]
    mol = chem.MolFromSmiles(smiles)
    row = _empty_row(smiles, valid=mol is not None, reason="" if mol is not None else "RDKit 无法解析")
    if mol is None:
        return row
    try:
        chem.SanitizeMol(mol)
    except Exception:
        row["valid"] = False
        row["reason"] = "结构无法标准化"
        return row

    descriptors = rdkit["Descriptors"]
    lipinski = rdkit["Lipinski"]
    qed = rdkit["QED"]
    row.update(
        {
            "canonical_smiles": chem.MolToSmiles(mol),
            "valid": True,
            "MW": _number(descriptors.MolWt(mol)),
            "LogP": _number(descriptors.MolLogP(mol)),
            "TPSA": _number(descriptors.TPSA(mol)),
            "HBD": _number(lipinski.NumHDonors(mol)),
            "HBA": _number(lipinski.NumHAcceptors(mol)),
            "RotBonds": _number(lipinski.NumRotatableBonds(mol)),
            "Fsp3": _number(descriptors.FractionCSP3(mol)),
            "QED": _number(qed.qed(mol)),
            "AromaticRings": _number(lipinski.NumAromaticRings(mol)),
            "Rings": _number(lipinski.RingCount(mol)),
            "HeavyAtoms": _number(lipinski.HeavyAtomCount(mol)),
            "alerts": _alerts(mol, chem),
        }
    )
    row["lipinski"] = bool(
        row["MW"] <= 500 and row["LogP"] <= 5 and row["HBD"] <= 5 and row["HBA"] <= 10
    )
    row["veber"] = bool(row["RotBonds"] <= 10 and row["TPSA"] <= 140)
    scorer = rdkit["sascorer"]
    if scorer is not None:
        try:
            row["SA"] = _number(scorer.calculateScore(mol))
        except Exception:
            row["SA"] = None
    else:
        bertz = rdkit["GraphDescriptors"].BertzCT(mol)
        row["SAProxy"] = _number(max(1.0, min(10.0, 1.0 + float(bertz) / 200.0)))
    return row


def _empty_row(smiles: str, valid: bool, reason: str) -> dict:
    return {
        "input_smiles": smiles,
        "canonical_smiles": "",
        "valid": valid,
        "reason": reason,
        "MW": None,
        "LogP": None,
        "TPSA": None,
        "HBD": None,
        "HBA": None,
        "RotBonds": None,
        "Fsp3": None,
        "QED": None,
        "AromaticRings": None,
        "Rings": None,
        "HeavyAtoms": None,
        "SA": None,
        "SAProxy": None,
        "IC50": None,
        "Ki": None,
        "Kd": None,
        "Affinity": None,
        "lipinski": None,
        "veber": None,
        "alerts": [],
        "numeric_status": "invalid" if not valid else "incomplete",
        "failed": [],
        "pending": [],
    }


def _alerts(mol, chem) -> list[str]:
    hits = []
    for name, smarts in _ALERTS:
        pattern = chem.MolFromSmarts(smarts)
        if pattern is not None and mol.HasSubstructMatch(pattern):
            hits.append(name)
    return hits


def _passes(value: float, constraint: Constraint) -> bool:
    if constraint.op == "<=":
        return value <= constraint.value
    if constraint.op == "<":
        return value < constraint.value
    if constraint.op == ">=":
        return value >= constraint.value
    if constraint.op == ">":
        return value > constraint.value
    if constraint.op == "==":
        return abs(value - constraint.value) <= 1e-6
    if constraint.op == "between":
        return constraint.value <= value <= constraint.high
    return False


def _find_molecule(molecules: list[dict], smiles: str) -> dict | None:
    for row in molecules:
        if smiles == row.get("input_smiles") or smiles == row.get("canonical_smiles"):
            return row
    return None


def _row_smiles(row: dict) -> str:
    for key, value in row.items():
        if str(key).strip().lower() in {"smiles", "smi", "canonical_smiles", "input_smiles"} and value:
            return str(value).strip()
    return ""


def _number(value):
    if value is None or value == "":
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number != number:
        return None
    return number


def _load_rdkit():
    try:
        from rdkit import Chem
        from rdkit.Chem import Descriptors, GraphDescriptors, Lipinski, QED
    except ImportError:
        return None
    return {
        "Chem": Chem,
        "Descriptors": Descriptors,
        "GraphDescriptors": GraphDescriptors,
        "Lipinski": Lipinski,
        "QED": QED,
        "sascorer": _load_sascorer(),
    }


def _load_sascorer():
    try:
        from rdkit.Contrib.SA_Score import sascorer

        return sascorer
    except Exception:
        pass
    try:
        import os
        import sys

        from rdkit.Chem import RDConfig

        folder = os.path.join(RDConfig.RDContribDir, "SA_Score")
        if folder not in sys.path:
            sys.path.append(folder)
        import sascorer

        return sascorer
    except Exception:
        return None
