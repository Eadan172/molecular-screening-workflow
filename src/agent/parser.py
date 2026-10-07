"""把用户的需求文本解析成靶点、理化阈值、分子和外部软件。"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from src.agent.models import Constraint, RequestSpec, ToolSpec


_SECTIONS = {
    "基本信息": "basic",
    "基本要求": "basic",
    "项目": "basic",
    "需求": "basic",
    "basic": "basic",
    "理化性质": "properties",
    "理化": "properties",
    "性质": "properties",
    "properties": "properties",
    "physicochemical": "properties",
    "活性与成药性": "activity",
    "活性": "activity",
    "成药性": "activity",
    "admet": "activity",
    "activity": "activity",
    "分子": "molecules",
    "化合物": "molecules",
    "smiles": "molecules",
    "molecules": "molecules",
    "计算软件": "software",
    "软件": "software",
    "本地软件": "software",
    "外部工具": "software",
    "tools": "software",
    "software": "software",
    "备注": "notes",
    "其他": "notes",
    "notes": "notes",
}

_METADATA = {
    "靶点": "target",
    "靶标": "target",
    "蛋白": "target",
    "target": "target",
    "protein": "target",
    "种属": "species",
    "物种": "species",
    "species": "species",
    "organism": "species",
    "适应症": "indication",
    "疾病": "indication",
    "indication": "indication",
    "disease": "indication",
    "目标": "goal",
    "目的": "goal",
    "goal": "goal",
    "objective": "goal",
    "分子文件": "molecule_file",
    "化合物文件": "molecule_file",
    "molecule_file": "molecule_file",
    "smiles_file": "molecule_file",
    "smilesfile": "molecule_file",
}

_PROPERTIES = {
    "分子量": "MW",
    "mw": "MW",
    "molwt": "MW",
    "molecularweight": "MW",
    "精确分子量": "MW",
    "logp": "LogP",
    "clogp": "LogP",
    "alogp": "LogP",
    "脂水分配系数": "LogP",
    "油水分配系数": "LogP",
    "tpsa": "TPSA",
    "极性表面积": "TPSA",
    "拓扑极性表面积": "TPSA",
    "氢键供体": "HBD",
    "氢键供体数": "HBD",
    "hbd": "HBD",
    "numhdonors": "HBD",
    "氢键受体": "HBA",
    "氢键受体数": "HBA",
    "hba": "HBA",
    "numhacceptors": "HBA",
    "可旋转键": "RotBonds",
    "可旋转键数": "RotBonds",
    "rotatablebonds": "RotBonds",
    "nrot": "RotBonds",
    "qed": "QED",
    "sa": "SA",
    "sascore": "SA",
    "合成可及性": "SA",
    "合成难度": "SA",
    "芳香环": "AromaticRings",
    "芳香环数": "AromaticRings",
    "环数": "Rings",
    "重原子": "HeavyAtoms",
    "重原子数": "HeavyAtoms",
    "fsp3": "Fsp3",
    "csp3": "Fsp3",
    "fractioncsp3": "Fsp3",
    "碳饱和度": "Fsp3",
    "ic50": "IC50",
    "ki": "Ki",
    "kd": "Kd",
    "herg": "hERG",
    "ames": "AMES",
    "对接结合能": "Affinity",
    "结合能": "Affinity",
    "bindingaffinity": "Affinity",
    "affinity": "Affinity",
}

_DEFAULT_OP = {
    "QED": ">=",
    "Fsp3": ">=",
}

_TOOL_KEYS = {
    "name": "name",
    "名称": "name",
    "type": "type",
    "类型": "type",
    "path": "path",
    "路径": "path",
    "命令": "path",
    "可执行文件": "path",
    "url": "url",
    "地址": "url",
    "api": "url",
    "接口": "url",
    "method": "method",
    "方法": "method",
    "args": "args",
    "参数": "args",
    "receptor": "receptor",
    "受体": "receptor",
    "timeout": "timeout",
    "超时": "timeout",
    "header": "header",
    "请求头": "header",
}

_AVOID = {"避免", "avoid", "阴性", "negative", "否", "无", "不需要", "不含"}
_REQUIRE = {"需要", "阳性", "positive", "是", "有"}
_SMILES_RE = re.compile(r"^[A-Za-z0-9@+\-\[\]\(\)=#$\\/%.]+$")
_RANGE_RE = re.compile(
    r"^\s*(-?\d+(?:\.\d+)?)\s*(?:-|~|到|至)\s*(-?\d+(?:\.\d+)?)\s*(.*)$"
)
_CMP_RE = re.compile(
    r"^\s*(<=|>=|<|>|==|=|≤|≥)\s*(-?\d+(?:\.\d+)?)\s*(.*)$"
)
_NUM_RE = re.compile(r"^\s*(-?\d+(?:\.\d+)?)\s*(.*)$")
_KV_RE = re.compile(r"^([^:：]{1,80})[:：]\s*(.*)$")


def parse_request(text: str) -> RequestSpec:
    spec = RequestSpec(raw_text=text or "")
    section = "basic"
    current: ToolSpec | None = None

    for raw_line in (text or "").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        new_section = _match_section(line)
        if new_section:
            if section == "software" and current is not None:
                spec.tools.append(_finalize_tool(current))
                current = None
            section = new_section
            continue

        if section == "software":
            current = _consume_tool_line(line, current, spec)
            continue

        if section == "molecules":
            _consume_molecule_line(line, spec)
            continue

        if section == "notes":
            spec.notes.append(line)
            continue

        matched = _KV_RE.match(line)
        if matched:
            _apply_kv(matched.group(1).strip(), matched.group(2).strip(), spec)
        else:
            spec.notes.append(line)

    if current is not None:
        spec.tools.append(_finalize_tool(current))
    spec.smiles = _dedupe(spec.smiles)
    return spec


def attach_molecule_file(spec: RequestSpec, search_dirs: list[Path]) -> None:
    """把「分子文件」里的 SMILES 追加进来。文件不存在时只记一笔说明。"""
    if not spec.molecule_file:
        return
    path = _find_file(spec.molecule_file, search_dirs)
    if path is None:
        spec.extras.append(f"分子文件不存在: {spec.molecule_file}")
        return
    loaded = _read_molecule_file(path)
    if not loaded:
        spec.extras.append(f"分子文件里没有读到 SMILES: {path}")
        return
    spec.smiles = _dedupe(spec.smiles + loaded)
    spec.extras.append(f"已从分子文件读入 {len(loaded)} 条 SMILES: {path}")


def _match_section(line: str) -> str | None:
    text = re.sub(r"^#+\s*", "", line).strip().strip("[]【】")
    text = re.sub(r"[:：]\s*$", "", text).strip()
    if not text or len(text) > 40 or any(ch in text for ch in ":："):
        return None
    return _SECTIONS.get(text.lower())


def _apply_kv(key: str, value: str, spec: RequestSpec) -> None:
    if not value:
        return
    meta = _METADATA.get(_norm_key(key))
    if meta == "target":
        spec.target = value
        return
    if meta == "species":
        spec.species = value
        return
    if meta == "indication":
        spec.indication = value
        return
    if meta == "goal":
        spec.goal = value
        return
    if meta == "molecule_file":
        spec.molecule_file = value
        return
    if _norm_key(key) in {"smiles", "分子", "化合物"}:
        spec.smiles.extend(_smiles_tokens(value))
        return
    spec.constraints.append(_parse_constraint(key, value))


def _consume_molecule_line(line: str, spec: RequestSpec) -> None:
    matched = _KV_RE.match(line)
    if matched:
        _apply_kv(matched.group(1).strip(), matched.group(2).strip(), spec)
        return
    token = _smiles_tokens(line)
    if token:
        spec.smiles.extend(token)
    else:
        spec.notes.append(line)


def _consume_tool_line(line: str, current: ToolSpec | None, spec: RequestSpec) -> ToolSpec | None:
    matched = _KV_RE.match(line)
    if not matched:
        spec.extras.append(line)
        return current
    key = _TOOL_KEYS.get(_norm_key(matched.group(1)))
    value = matched.group(2).strip()
    if key is None:
        spec.extras.append(line)
        return current
    if key == "name":
        if current is not None:
            spec.tools.append(_finalize_tool(current))
        return ToolSpec(name=value or "tool")
    if current is None:
        current = ToolSpec(name=f"tool_{len(spec.tools) + 1}")
    if key == "type":
        current.type = _normalize_tool_type(value)
    elif key == "timeout":
        numbers = re.search(r"\d+", value)
        current.timeout = int(numbers.group()) if numbers else 300
        current.timeout = min(max(current.timeout, 1), 86400)
    elif key == "method":
        current.method = value.upper() or "POST"
    else:
        setattr(current, key, value)
    return current


def _finalize_tool(tool: ToolSpec) -> ToolSpec:
    if not tool.name:
        tool.name = "tool"
    if tool.type not in {"api", "local"}:
        tool.type = "api" if tool.url else "local"
    if tool.url and not tool.path and tool.type != "api":
        tool.type = "api"
    tool.method = (tool.method or "POST").upper()
    if tool.method not in {"GET", "POST"}:
        tool.method = "POST"
    tool.timeout = min(max(int(tool.timeout or 300), 1), 86400)
    return tool


def _normalize_tool_type(value: str) -> str:
    key = _norm_key(value)
    if key in {"api", "http", "https", "接口", "远程"}:
        return "api"
    return "local"


def _parse_constraint(key: str, value: str) -> Constraint:
    cleaned = value.translate(str.maketrans({"＜": "<", "＞": ">", "＝": "=", "－": "-", "—": "-", "～": "~"}))
    prop = _PROPERTIES.get(_norm_key(key))
    name = prop or key.strip()
    base = Constraint(name=name, kind="text", op="note", raw=value.strip(), source_key=key.strip())
    lowered = cleaned.strip().lower()
    if lowered in _AVOID:
        base.op = "avoid"
        return base
    if lowered in _REQUIRE:
        base.op = "require"
        return base

    ranged = _RANGE_RE.match(cleaned)
    if ranged:
        low = float(ranged.group(1))
        high = float(ranged.group(2))
        if low > high:
            low, high = high, low
        base.kind = "numeric"
        base.op = "between"
        base.value = low
        base.high = high
        base.unit = ranged.group(3).strip()
        return base

    compared = _CMP_RE.match(cleaned)
    if compared:
        op = {"≤": "<=", "≥": ">=", "=": "=="}.get(compared.group(1), compared.group(1))
        base.kind = "numeric"
        base.op = op
        base.value = float(compared.group(2))
        base.unit = compared.group(3).strip()
        return base

    number = _NUM_RE.match(cleaned)
    if number and prop:
        base.kind = "numeric"
        base.op = _DEFAULT_OP.get(prop, "<=")
        base.value = float(number.group(1))
        base.unit = number.group(2).strip()
        return base
    return base


def _smiles_tokens(value: str) -> list[str]:
    tokens = []
    for part in re.split(r"[;；,，]+", value.strip()):
        piece = part.strip()
        if not piece:
            continue
        token = piece.split()[0]
        if _is_smiles_token(token):
            tokens.append(token)
    return tokens


def _is_smiles_token(token: str) -> bool:
    if not token or len(token) > 500:
        return False
    if not _SMILES_RE.fullmatch(token):
        return False
    if not re.search(r"[A-Za-z]", token):
        return False
    if token.isalpha() and len(token) < 2:
        return False
    return True


def _read_molecule_file(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8-sig")
    lines = text.splitlines()
    if not lines or not any(line.strip() for line in lines):
        return []
    first = next(line for line in lines if line.strip())
    if path.suffix.lower() == ".csv" or ("," in first and not first.strip().startswith("#")):
        rows = list(csv.DictReader(lines))
        if rows:
            column = _smiles_column(rows[0])
            if column:
                found = []
                for row in rows:
                    found.extend(_smiles_tokens(row.get(column) or ""))
                return _dedupe(found)
    found = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        found.extend(_smiles_tokens(line))
    return _dedupe(found)


def _smiles_column(row: dict) -> str | None:
    for key in row:
        if key and key.strip().lower() in {"smiles", "smi", "canonical_smiles", "molecule"}:
            return key
    return None


def _find_file(raw: str, search_dirs: list[Path]) -> Path | None:
    candidate = Path(raw)
    if candidate.is_file():
        return candidate.resolve()
    for directory in search_dirs:
        path = (directory / raw).resolve()
        if path.is_file():
            return path
    return None


def _dedupe(items: list[str]) -> list[str]:
    seen = set()
    ordered = []
    for item in items:
        if item in seen:
            continue
        seen.add(item)
        ordered.append(item)
    return ordered


def _norm_key(key: str) -> str:
    return re.sub(r"[\s_\-]+", "", key).lower()
