"""把要求分析、拆解、调研、计算、数据分析和报告串起来。"""

from __future__ import annotations

import json
import os
import re
from datetime import datetime
from pathlib import Path

from src.agent.calculate import annotate, calculate_molecules, merge_external
from src.agent.jev import JevClient
from src.agent.llm import LLMClient
from src.agent.parser import attach_molecule_file, parse_request
from src.agent.report import render_report, rule_findings, write_properties
from src.agent.research import lookup_target
from src.agent.settings import Settings
from src.agent.tools import run_tools


_SECRET_LINE = re.compile(
    r"(?i)(authorization|api[_-]?key|token|secret|请求头|header)\s*[:：=]"
)


def run_pipeline(
    request_path: Path,
    settings: Settings,
    output_root: Path,
    project_root: Path | None = None,
    llm=None,
    jev=None,
    http_get=None,
    http_request=None,
    run_id: str | None = None,
) -> dict:
    request_path = Path(request_path)
    if not request_path.is_file():
        raise FileNotFoundError(f"找不到输入文件: {request_path}")

    root = Path(project_root or request_path.resolve().parent)
    text = request_path.read_text(encoding="utf-8-sig")
    spec = parse_request(text)
    attach_molecule_file(spec, [request_path.resolve().parent, root])

    log: list[str] = []

    def say(message: str) -> None:
        print(message)
        log.append(message)

    use_llm = bool(settings.llm_enabled)
    if use_llm and llm is None:
        llm = LLMClient(settings)
    if use_llm and llm is not None and not getattr(llm, "enabled", True):
        use_llm = False
    use_jev = bool(settings.jev_enabled)
    if use_jev and jev is None:
        jev = JevClient(settings)
    if use_jev and jev is not None and not getattr(jev, "enabled", True):
        use_jev = False

    llm_errors: list[str] = []
    jev_errors: list[str] = []
    analysis: dict = {
        "requirement_analysis": {},
        "research_llm": {},
        "findings_llm": {},
        "llm_errors": llm_errors,
    }

    say("1/6 要求分析")
    if use_llm:
        parsed = _llm_json(
            llm,
            "analyze",
            {
                "target": spec.target,
                "species": spec.species,
                "indication": spec.indication,
                "goal": spec.goal,
                "notes": spec.notes,
                "constraints": [item.display() for item in spec.constraints],
                "n_smiles": len(spec.smiles),
                "smiles_preview": spec.smiles[:30],
                "tools": spec.public_dict()["tools"],
                "request_excerpt": _redact(spec.raw_text)[:8000],
            },
            llm_errors,
        )
        if parsed:
            analysis["requirement_analysis"] = parsed
            say("已用大模型理解需求。")
        else:
            say("大模型没有返回需求理解，继续使用规则解析。")
    else:
        say("未调用大模型，按输入文本做规则解析。")

    say("2/6 任务拆解")
    say(f"分子 {len(spec.smiles)} 个，数值/文字要求 {len(spec.constraints)} 条，外部软件 {len(spec.tools)} 个。")

    say("3/6 调研")
    research = lookup_target(spec, http_get=http_get)
    say(research.get("summary") or "调研结束。")
    if use_llm:
        parsed = _llm_json(
            llm,
            "research",
            {
                "target": spec.target,
                "species": spec.species,
                "indication": spec.indication,
                "goal": spec.goal,
                "database_summary": research.get("summary") or "",
                "chosen_target": research.get("chosen"),
                "activity_stats": research.get("activity_stats"),
            },
            llm_errors,
        )
        if parsed:
            analysis["research_llm"] = parsed

    say("4/6 计算")
    calc = calculate_molecules(spec)
    run_dir = _allocate_run_dir(output_root, run_id)
    tool_results = run_tools(spec, run_dir / "tools", root, http_request=http_request)
    tables = []
    for tool in tool_results:
        tables.extend(tool.get("tables") or [])
        label = {"done": "完成", "skipped": "跳过", "failed": "失败"}.get(tool["status"], tool["status"])
        say(f"外部软件 {tool['name']}: {label}" + (f"（{tool['reason']}）" if tool.get("reason") else ""))
    matched = merge_external(calc["molecules"], tables)
    calc["summary"]["external_rows_matched"] = matched
    if matched:
        note = f"已按 SMILES 对齐 {matched} 条外部计算结果。"
        calc["summary"]["notes"].append(note)
        say(note)
    elif tables:
        note = "外部结果里有表格，但没有和输入 SMILES 对上。"
        calc["summary"]["notes"].append(note)
        say(note)
    annotate(spec, calc["molecules"], calc["summary"])

    say("5/6 计算数据分析")
    findings = rule_findings(spec, calc, tool_results)
    jev_decision = {"status": "disabled", "answers": {}, "policy": {}}
    if use_jev:
        try:
            jev_decision = jev.decide(_jev_state(spec, research, calc, tool_results))
            policy_route = (jev_decision.get("policy") or {}).get("route") or "未返回"
            say(f"Jev 决策门控完成，路由: {policy_route}")
        except Exception as exc:
            jev_errors.append(str(exc))
            jev_decision = {
                "status": "failed",
                "error": str(exc),
                "answers": {},
                "policy": {
                    "route": "human_review",
                    "allow_automatic_progress": False,
                    "reasons": ["Jev 调用失败，保守转人工复核。"],
                },
            }
            say("Jev 调用失败，保守转人工复核。")
    if use_llm:
        parsed = _llm_json(
            llm,
            "findings",
            {
                "counts": {
                    "n_input": calc["summary"]["n_input"],
                    "n_valid": calc["summary"]["n_valid"],
                    "n_passed": calc["summary"]["n_passed"],
                    "n_failed": calc["summary"]["n_failed"],
                    "n_incomplete": calc["summary"]["n_incomplete"],
                    "n_invalid": calc["summary"]["n_invalid"],
                },
                "fail_counts": calc["summary"].get("fail_counts"),
                "text_requirements": calc["summary"].get("text_requirements"),
                "not_applied": calc["summary"].get("not_applied"),
                "notes": calc["summary"].get("notes"),
                "tools": [
                    {
                        "name": tool["name"],
                        "type": tool["type"],
                        "status": tool["status"],
                        "reason": tool.get("reason") or "",
                        "parsed_rows": tool.get("parsed_rows", 0),
                    }
                    for tool in tool_results
                ],
                "molecules_preview": [
                    {
                        "smiles": row.get("input_smiles"),
                        "status": row.get("numeric_status"),
                        "failed": row.get("failed"),
                        "pending": row.get("pending"),
                        "MW": row.get("MW"),
                        "LogP": row.get("LogP"),
                        "QED": row.get("QED"),
                    }
                    for row in calc["molecules"][:15]
                ],
            },
            llm_errors,
        )
        if parsed:
            analysis["findings_llm"] = parsed

    mode_parts = []
    if use_llm:
        mode_parts.append("LLM 增强" + ("（部分回退）" if llm_errors else ""))
    else:
        mode_parts.append("仅计算（未配置 LLM API Key）")
    if use_jev:
        mode_parts.append("Jev 决策" + ("（失败转人工）" if jev_errors else ""))
    else:
        mode_parts.append("未启用 Jev")
    mode = " + ".join(mode_parts)

    tasks = _tasks(spec, research, calc, tool_results, jev_decision)
    say("6/6 报告整理")
    report_path = run_dir / "report.md"
    properties_path = run_dir / "properties.csv"
    context = {
        "spec": spec,
        "calc": calc,
        "analysis": analysis,
        "research": research,
        "jev": jev_decision,
        "jev_enabled": use_jev,
        "jev_errors": jev_errors,
        "rule_findings": findings,
        "tools": tool_results,
        "tasks": tasks,
        "mode": mode,
        "run_id": run_dir.name,
        "input_path": str(request_path),
        "report_name": "report.md",
        "properties_name": "properties.csv",
    }
    report_path.write_text(render_report(context), encoding="utf-8")
    write_properties(properties_path, calc["molecules"])
    _write_json(run_dir / "spec.json", spec.public_dict())
    _write_json(run_dir / "research.json", research)
    _write_json(run_dir / "tasks.json", tasks)
    _write_json(
        run_dir / "calculations.json",
        {
            "summary": calc["summary"],
            "molecules": calc["molecules"],
            "tools": [_public_tool(tool) for tool in tool_results],
        },
    )
    _write_json(
        run_dir / "analysis.json",
        {
            "mode": mode,
            "llm_errors": llm_errors,
            "jev_errors": jev_errors,
            "jev": jev_decision,
            "rule_findings": findings,
            **analysis,
        },
    )
    say(f"报告已写入 {report_path}")
    (run_dir / "run.log").write_text("\n".join(log) + "\n", encoding="utf-8")
    return {
        "report": str(report_path),
        "properties": str(properties_path),
        "run_dir": str(run_dir),
        "mode": mode,
    }


def _tasks(spec, research: dict, calc: dict, tools: list[dict], jev: dict) -> list[dict]:
    summary = calc["summary"]
    if not spec.smiles:
        calc_status, calc_detail = "skipped", "输入里没有分子。"
    elif summary.get("rdkit"):
        calc_status, calc_detail = "done", f"已计算 {summary.get('n_valid', 0)} 个有效结构。"
    else:
        calc_status, calc_detail = "skipped", "环境里没有 RDKit。"
    tasks = [
        {"name": "要求分析", "status": "done", "detail": "已提取靶点、种属、阈值和分子。"},
        {
            "name": "任务拆解",
            "status": "done",
            "detail": "执行哪些计算由输入文件决定。模型只能补充建议，不能新增命令。",
        },
        {
            "name": "调研",
            "status": "done" if research.get("summary") else "skipped",
            "detail": (research.get("summary") or "").split("\n")[0][:180],
        },
        {"name": "理化性质计算", "status": calc_status, "detail": calc_detail},
    ]
    for tool in tools:
        tasks.append(
            {
                "name": f"外部计算：{tool['name']}",
                "status": tool["status"],
                "detail": tool.get("reason") or f"解析到 {tool.get('parsed_rows', 0)} 行结果。",
            }
        )
    tasks.append({"name": "计算数据分析", "status": "done", "detail": "已按数值阈值标注通过、未通过和数据不足。"})
    if jev.get("status") == "done":
        tasks.append(
            {
                "name": "Jev 决策门控",
                "status": "done",
                "detail": f"建议路由：{(jev.get('policy') or {}).get('route') or '未返回'}。",
            }
        )
    elif jev.get("status") == "failed":
        tasks.append({"name": "Jev 决策门控", "status": "failed", "detail": "调用失败，已保守转人工复核。"})
    else:
        tasks.append({"name": "Jev 决策门控", "status": "skipped", "detail": "未配置 JEV_API_KEY。"})
    tasks.append({"name": "报告整理", "status": "done", "detail": "已写入 report.md 和 properties.csv。"})
    return tasks


def _jev_state(spec, research: dict, calc: dict, tools: list[dict]) -> dict:
    """只给 Jev 决策所需的摘要，避免密钥、原始请求头和大批分子结构外发。"""
    summary = calc["summary"]
    return {
        "project": {
            "target": spec.target,
            "species": spec.species,
            "indication": spec.indication,
            "goal": spec.goal,
            "constraints": [item.display() for item in spec.constraints],
            "notes": spec.notes,
        },
        "evidence": {
            "database_summary": research.get("summary") or "",
            "chosen_target": research.get("chosen"),
            "activity_stats": research.get("activity_stats") or {},
            "database_error": research.get("error") or "",
        },
        "calculation": {
            "n_input": summary.get("n_input"),
            "n_valid": summary.get("n_valid"),
            "n_passed": summary.get("n_passed"),
            "n_failed": summary.get("n_failed"),
            "n_incomplete": summary.get("n_incomplete"),
            "n_invalid": summary.get("n_invalid"),
            "fail_counts": summary.get("fail_counts") or {},
            "not_applied": summary.get("not_applied") or [],
            "text_requirements": summary.get("text_requirements") or [],
            "external_rows_matched": summary.get("external_rows_matched") or 0,
        },
        "tools": [
            {
                "name": tool.get("name"),
                "type": tool.get("type"),
                "status": tool.get("status"),
                "reason": tool.get("reason") or "",
                "parsed_rows": tool.get("parsed_rows") or 0,
            }
            for tool in tools
        ],
        "decision_scope": (
            "Jev 只做证据充分性、项目模态和下一步路由。"
            "SMILES 有效性、立体化学和数值性质由确定性化学工具判断。"
        ),
    }


def _llm_json(llm, task: str, payload: dict, errors: list[str]) -> dict | None:
    try:
        data = llm.complete_json(task, payload)
    except Exception as exc:
        errors.append(f"{task}: {exc}")
        return None
    if not isinstance(data, dict):
        errors.append(f"{task}: 返回内容不是对象")
        return None
    return data


def _allocate_run_dir(output_root: Path, run_id: str | None) -> Path:
    output_root = Path(output_root)
    output_root.mkdir(parents=True, exist_ok=True)
    stamp = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
    stamp = re.sub(r"[^\w\-]", "_", stamp) or "run"
    path = output_root / stamp
    if path.exists():
        path = output_root / f"{stamp}_{os.getpid()}"
    path.mkdir(parents=True, exist_ok=False)
    return path


def _public_tool(tool: dict) -> dict:
    return {key: value for key, value in tool.items() if key != "tables"}


def _write_json(path: Path, payload) -> None:
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")


def _redact(text: str) -> str:
    lines = []
    for line in (text or "").splitlines():
        if _SECRET_LINE.search(line):
            lines.append("# [已省略可能包含密钥的一行]")
        else:
            lines.append(line)
    return "\n".join(lines)
