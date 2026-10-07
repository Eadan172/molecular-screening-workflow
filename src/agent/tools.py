"""按输入文件调用本地计算程序或 HTTP 接口。"""

from __future__ import annotations

import csv
import json
import os
import re
import shlex
import shutil
import subprocess
from pathlib import Path

from src.agent.models import RequestSpec, ToolSpec


def run_tools(spec: RequestSpec, workdir: Path, root: Path, http_request=None) -> list[dict]:
    results = []
    used_names: dict[str, int] = {}
    for tool in spec.tools:
        folder = _unique_name(tool.name, used_names)
        tool_dir = workdir / folder
        tool_dir.mkdir(parents=True, exist_ok=True)
        if tool.type == "api":
            results.append(_run_api(tool, spec, tool_dir, http_request))
        else:
            results.append(_run_local(tool, spec, tool_dir, root))
    return results


def _run_local(tool: ToolSpec, spec: RequestSpec, workdir: Path, root: Path) -> dict:
    executable = _resolve_executable(tool.path, root)
    if executable is None:
        return _result(tool, "skipped", f"找不到可执行文件: {tool.path or '（未填写路径）'}", workdir)

    receptor = ""
    if tool.receptor:
        receptor_path = _resolve_file(tool.receptor, root)
        if receptor_path is None:
            return _result(tool, "skipped", f"找不到受体文件: {tool.receptor}", workdir)
        receptor = str(receptor_path)

    ligand = workdir / "ligands.sdf"
    smiles_file = workdir / "ligands.smi"
    output_file = workdir / "result.out"
    if spec.smiles:
        smiles_file.write_text("\n".join(spec.smiles) + "\n", encoding="utf-8")
        wrote = _write_sdf(spec.smiles, ligand)
        if "{ligand}" in tool.args and not wrote:
            return _result(tool, "skipped", "参数使用了 {ligand}，但当前环境不能把 SMILES 写成 SDF。", workdir)

    mapping = {
        "ligand": str(ligand),
        "ligand_sdf": str(ligand),
        "smiles": str(smiles_file),
        "ligand_smi": str(smiles_file),
        "receptor": receptor,
        "output": str(output_file),
        "outdir": str(workdir),
        "workdir": str(workdir),
    }
    if "{receptor}" in tool.args and not receptor:
        return _result(tool, "skipped", "参数使用了 {receptor}，但没有提供受体文件。", workdir)
    if "{ligand}" in tool.args and not spec.smiles:
        return _result(tool, "skipped", "参数使用了 {ligand}，但输入里没有分子。", workdir)

    try:
        tokens = shlex.split(tool.args or "", posix=(os.name != "nt"))
    except ValueError as exc:
        return _result(tool, "failed", f"参数无法解析: {exc}", workdir)
    command = [str(executable)] + [_substitute(token, mapping) for token in tokens]
    return _execute(tool, command, workdir)


def _run_api(tool: ToolSpec, spec: RequestSpec, workdir: Path, http_request) -> dict:
    url = (tool.url or "").strip()
    if not url.lower().startswith(("http://", "https://")):
        return _result(tool, "skipped", "接口地址需要以 http:// 或 https:// 开头。", workdir)
    header_name, header_value = _split_header(tool.header)
    headers = {"Content-Type": "application/json", "User-Agent": "mol-workflow/1.0"}
    if header_name:
        headers[header_name] = header_value
    payload = {
        "tool": tool.name,
        "target": spec.target,
        "species": spec.species,
        "indication": spec.indication,
        "goal": spec.goal,
        "constraints": [
            {"name": item.name, "op": item.op, "value": item.value, "high": item.high, "unit": item.unit}
            for item in spec.constraints
        ],
        "molecules": [{"smiles": smiles} for smiles in spec.smiles],
    }
    request = http_request or _default_request
    try:
        response = request(tool.method, url, json=payload, headers=headers, timeout=tool.timeout)
    except Exception as exc:
        return _result(tool, "failed", f"接口调用失败: {exc}"[:300], workdir)

    status_code = getattr(response, "status_code", 0)
    content = getattr(response, "content", b"") or b""
    if isinstance(content, str):
        content = content.encode("utf-8")
    content = content[:2_000_000]
    body_path = workdir / "response.json"
    body_path.write_bytes(content)
    parsed = _parse_saved_files(workdir)
    record = _result(tool, "done" if status_code < 400 else "failed", "", workdir)
    record["command"] = [tool.method, url]
    record["exit_code"] = status_code
    record["parsed_rows"] = len(parsed)
    record["tables"] = parsed
    if status_code >= 400:
        record["reason"] = f"接口返回 {status_code}"
        record["stderr_tail"] = _sanitize(content[:500].decode("utf-8", errors="replace"))
    return record


def _execute(tool: ToolSpec, command: list[str], workdir: Path) -> dict:
    env = os.environ.copy()
    env.pop("LLM_API_KEY", None)
    env.pop("OPENAI_API_KEY", None)
    try:
        completed = subprocess.run(
            command,
            cwd=workdir,
            env=env,
            capture_output=True,
            timeout=tool.timeout,
            check=False,
        )
    except subprocess.TimeoutExpired:
        return _result(tool, "failed", f"运行超过 {tool.timeout} 秒", workdir, command)
    except OSError as exc:
        return _result(tool, "failed", f"无法启动程序: {exc}", workdir, command)

    (workdir / "stdout.txt").write_bytes(completed.stdout or b"")
    (workdir / "stderr.txt").write_bytes(completed.stderr or b"")
    parsed = _parse_saved_files(workdir)
    record = _result(
        tool,
        "done" if completed.returncode == 0 else "failed",
        "" if completed.returncode == 0 else f"退出码 {completed.returncode}",
        workdir,
        command,
    )
    record["exit_code"] = completed.returncode
    record["stdout_tail"] = _sanitize((completed.stdout or b"")[:500].decode("utf-8", errors="replace"))
    record["stderr_tail"] = _sanitize((completed.stderr or b"")[:500].decode("utf-8", errors="replace"))
    record["parsed_rows"] = len(parsed)
    record["tables"] = parsed
    return record


def _result(tool: ToolSpec, status: str, reason: str, workdir: Path, command: list[str] | None = None) -> dict:
    return {
        "name": tool.name,
        "type": tool.type,
        "status": status,
        "reason": _sanitize(reason),
        "command": command or [],
        "exit_code": None,
        "stdout_tail": "",
        "stderr_tail": "",
        "output_dir": str(workdir),
        "parsed_rows": 0,
        "tables": [],
    }


def _parse_saved_files(workdir: Path) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(workdir.rglob("*")):
        if not path.is_file():
            continue
        suffix = path.suffix.lower()
        try:
            if suffix == ".csv":
                rows.extend(_read_csv(path))
            elif suffix == ".json":
                rows.extend(_read_json(path))
        except (OSError, UnicodeError, csv.Error, json.JSONDecodeError, ValueError):
            continue
    return rows


def _read_csv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return [dict(row) for row in csv.DictReader(handle)]


def _read_json(path: Path) -> list[dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(data, list):
        return [item for item in data if isinstance(item, dict)]
    if isinstance(data, dict):
        for key in ("molecules", "results", "data"):
            value = data.get(key)
            if isinstance(value, list):
                return [item for item in value if isinstance(item, dict)]
        if any(str(key).lower() in {"smiles", "smi"} for key in data):
            return [data]
    return []


def _write_sdf(smiles_list: list[str], path: Path) -> bool:
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem
    except ImportError:
        return False
    writer = Chem.SDWriter(str(path))
    wrote = False
    try:
        for index, smiles in enumerate(smiles_list, start=1):
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                continue
            mol.SetProp("_Name", f"mol_{index}")
            mol.SetProp("SMILES", smiles)
            AllChem.Compute2DCoords(mol)
            writer.write(mol)
            wrote = True
    finally:
        writer.close()
    return wrote


def _resolve_executable(raw: str, root: Path) -> Path | None:
    if not raw:
        return None
    return _resolve_file(raw, root) or (_as_path(shutil.which(raw)))


def _resolve_file(raw: str, root: Path) -> Path | None:
    candidate = Path(raw).expanduser()
    if candidate.is_file():
        return candidate.resolve()
    rooted = (root / raw).expanduser()
    if rooted.is_file():
        return rooted.resolve()
    return None


def _as_path(value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    return path if path.is_file() else None


def _substitute(token: str, mapping: dict[str, str]) -> str:
    import re

    def replace(match):
        return mapping.get(match.group(1), match.group(0))

    return re.sub(r"\{([A-Za-z0-9_]+)\}", replace, token)


def _sanitize(text: str) -> str:
    cleaned = re.sub(r"(?i)(authorization\s*[:=]\s*bearer\s+)\S+", r"\1[已省略]", text or "")
    cleaned = re.sub(r"(?i)(bearer\s+)[A-Za-z0-9._\-]{12,}", r"\1[已省略]", cleaned)
    cleaned = re.sub(r"(?i)((?:api[_-]?key|token|secret)\s*[=:]\s*)\S+", r"\1[已省略]", cleaned)
    return cleaned


def _split_header(value: str) -> tuple[str, str]:
    if not value or ":" not in value:
        return "", ""
    name, content = value.split(":", 1)
    return name.strip(), content.strip()


def _unique_name(name: str, used: dict[str, int]) -> str:
    import re

    safe = re.sub(r"[^\w\-]+", "_", name).strip("_") or "tool"
    count = used.get(safe, 0) + 1
    used[safe] = count
    return safe if count == 1 else f"{safe}_{count}"


def _default_request(method, url, json=None, headers=None, timeout=None):
    import requests

    return requests.request(method, url, json=json, headers=headers, timeout=timeout)
