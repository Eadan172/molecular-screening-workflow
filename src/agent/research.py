"""用公开的 ChEMBL 接口做靶点检索。失败时不中断计算。"""

from __future__ import annotations

from src.agent.models import RequestSpec


CHEMBL_TARGET = "https://www.ebi.ac.uk/chembl/api/data/target/search.json"
CHEMBL_ACTIVITY = "https://www.ebi.ac.uk/chembl/api/data/activity.json"

_SPECIES = {
    "human": "homo sapiens",
    "人": "homo sapiens",
    "人类": "homo sapiens",
    "mouse": "mus musculus",
    "小鼠": "mus musculus",
    "rat": "rattus norvegicus",
    "大鼠": "rattus norvegicus",
}


def lookup_target(spec: RequestSpec, http_get=None, timeout: int = 20) -> dict:
    result = {
        "query": spec.target,
        "species": spec.species,
        "targets": [],
        "chosen": None,
        "activities": [],
        "activity_stats": {},
        "error": "",
        "summary": "",
    }
    if not spec.target:
        result["summary"] = "输入里没有靶点名称，已跳过公开数据库检索。"
        return result

    getter = http_get or _default_get
    try:
        payload = getter(CHEMBL_TARGET, {"q": spec.target, "limit": 20}, timeout)
        targets = [_target_fields(item) for item in _as_list(payload, "targets")]
        targets = [item for item in targets if item["chembl_id"] or item["pref_name"]]
        result["targets"] = targets
    except Exception as exc:  # 网络或接口异常都不应挡住后续计算
        result["error"] = str(exc)[:300]
        result["summary"] = "公开数据库检索失败，已跳过。本地计算仍会继续。"
        return result

    if not result["targets"]:
        result["summary"] = f"ChEMBL 没有返回与「{spec.target}」匹配的靶点。"
        return result

    chosen = _choose_target(result["targets"], spec.species, spec.target)
    result["chosen"] = chosen
    if chosen.get("chembl_id"):
        try:
            activity_payload = getter(
                CHEMBL_ACTIVITY,
                {
                    "target_chembl_id": chosen["chembl_id"],
                    "standard_type": "IC50",
                    "limit": 20,
                },
                timeout,
            )
            activities = []
            for item in _as_list(activity_payload, "activities")[:20]:
                activities.append(
                    {
                        "molecule_chembl_id": item.get("molecule_chembl_id") or "",
                        "canonical_smiles": item.get("canonical_smiles") or "",
                        "standard_type": item.get("standard_type") or "IC50",
                        "standard_value": _as_float(item.get("standard_value")),
                        "standard_units": item.get("standard_units") or "",
                    }
                )
            result["activities"] = activities
            result["activity_stats"] = _activity_stats(activities)
        except Exception as exc:
            result["error"] = str(exc)[:300]

    result["summary"] = _summarize(spec, result)
    return result


def _default_get(url, params, timeout):
    import requests

    response = requests.get(
        url,
        params=params,
        timeout=timeout,
        headers={"Accept": "application/json", "User-Agent": "mol-workflow/1.0"},
    )
    response.raise_for_status()
    data = response.json()
    if not isinstance(data, dict):
        raise RuntimeError("ChEMBL 返回的不是 JSON 对象")
    return data


def _target_fields(item: dict) -> dict:
    return {
        "chembl_id": item.get("target_chembl_id") or item.get("chembl_id") or "",
        "pref_name": item.get("pref_name") or item.get("name") or "",
        "organism": item.get("organism") or "",
        "target_type": item.get("target_type") or "",
    }


def _as_list(payload, key: str) -> list:
    if not isinstance(payload, dict):
        return []
    value = payload.get(key) or []
    return value if isinstance(value, list) else []


def _choose_target(targets: list[dict], species: str, query: str) -> dict:
    ranked = sorted(
        enumerate(targets),
        key=lambda item: (-_target_score(item[1], species, query), item[0]),
    )
    return ranked[0][1]


def _target_score(item: dict, species: str, query: str) -> int:
    name = (item.get("pref_name") or "").lower()
    query_text = (query or "").strip().lower()
    organism = item.get("organism") or ""
    target_type = (item.get("target_type") or "").upper()
    score = 0
    if species and organism:
        score += 5 if species_matches(species, organism) else -3
    if "SINGLE" in target_type:
        score += 3
    if "COMPLEX" in target_type or "PROTEIN-PROTEIN" in target_type:
        score -= 2
    tokens = set(name.replace("/", " ").replace("-", " ").split())
    if query_text and query_text == name:
        score += 8
    elif query_text and query_text in tokens:
        score += 4
    if "substrate" in name and "substrate" not in query_text:
        score -= 6
    return score


def species_matches(requested: str, organism: str) -> bool:
    if not requested or not organism:
        return True
    left = requested.strip().lower()
    right = organism.strip().lower()
    if left in right or right in left:
        return True
    mapped = _SPECIES.get(left, left)
    return mapped in right or right in mapped


def _activity_stats(activities: list[dict]) -> dict:
    values = [item["standard_value"] for item in activities if item["standard_value"] is not None]
    units = [item["standard_units"] for item in activities if item["standard_units"]]
    stats = {"n": len(activities), "n_with_value": len(values)}
    if values:
        stats["min"] = min(values)
        stats["max"] = max(values)
    if units:
        stats["units"] = max(set(units), key=units.count)
    return stats


def _summarize(spec: RequestSpec, result: dict) -> str:
    chosen = result.get("chosen") or {}
    lines = [
        (
            f"ChEMBL 匹配到 {chosen.get('pref_name') or spec.target}"
            f"（{chosen.get('chembl_id') or '无编号'}），"
            f"种属 {chosen.get('organism') or '未知'}，"
            f"类型 {chosen.get('target_type') or '未知'}。"
        )
    ]
    if spec.species and chosen.get("organism") and not species_matches(spec.species, chosen["organism"]):
        lines.append(
            f"输入种属是 {spec.species}，与这条记录的种属不一致。下面的活性只作背景，不能直接当成该种属的数据。"
        )
    elif spec.species and chosen.get("organism"):
        lines.append(f"记录种属与输入的 {spec.species} 一致。")

    stats = result.get("activity_stats") or {}
    if stats.get("n"):
        unit = stats.get("units") or "未标单位"
        if "min" in stats:
            lines.append(
                f"取回 {stats['n']} 条 IC50 记录，其中 {stats['n_with_value']} 条有数值，"
                f"范围 {stats['min']:g}–{stats['max']:g} {unit}。"
            )
        else:
            lines.append(f"取回 {stats['n']} 条 IC50 记录，但没有可用数值。")
        examples = [
            item.get("canonical_smiles")
            for item in result.get("activities") or []
            if item.get("canonical_smiles")
        ][:3]
        if examples:
            lines.append("示例分子：" + "；".join(examples))
    else:
        lines.append("没有取到 IC50 记录。")
    if result.get("error"):
        lines.append(f"活性检索没有完成：{result['error']}")
    lines.append("这是公开数据库摘要，不是本次分子的活性结论。")
    return "\n".join(lines)


def _as_float(value):
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
