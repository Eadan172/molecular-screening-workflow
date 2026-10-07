"""把各步结果整理成可读的 Markdown 报告和 CSV。"""

from __future__ import annotations

import csv
from pathlib import Path


STATUS_LABEL = {
    "pass": "通过",
    "fail": "未通过",
    "incomplete": "数据不足",
    "invalid": "结构无效",
    "done": "完成",
    "skipped": "跳过",
    "failed": "失败",
}

_CSV_COLUMNS = [
    ("input_smiles", "输入SMILES"),
    ("canonical_smiles", "标准SMILES"),
    ("numeric_status", "数值要求"),
    ("failed", "未通过项"),
    ("pending", "缺少数据"),
    ("MW", "分子量"),
    ("LogP", "LogP"),
    ("TPSA", "TPSA"),
    ("HBD", "氢键供体"),
    ("HBA", "氢键受体"),
    ("RotBonds", "可旋转键"),
    ("Fsp3", "Fsp3"),
    ("QED", "QED"),
    ("SA", "SA"),
    ("SAProxy", "SA代理"),
    ("IC50", "IC50"),
    ("Ki", "Ki"),
    ("Kd", "Kd"),
    ("Affinity", "结合能"),
    ("lipinski", "Lipinski参考"),
    ("veber", "Veber参考"),
    ("alerts", "结构警示"),
    ("reason", "说明"),
]


def rule_findings(spec, calc: dict, tools: list[dict]) -> dict:
    summary = calc["summary"]
    molecules = calc["molecules"]
    lines = []
    if not spec.smiles:
        lines.append("没有提供分子，本次没有做性质计算。")
    elif not summary.get("rdkit"):
        lines.append("分子已经读入，但环境里没有 RDKit，理化性质没有算出来。")
    else:
        lines.append(
            f"共 {summary['n_input']} 个分子，有效结构 {summary['n_valid']} 个："
            f"数值要求通过 {summary['n_passed']} 个，未通过 {summary['n_failed']} 个，"
            f"数据不足 {summary['n_incomplete']} 个，结构无效 {summary['n_invalid']} 个。"
        )
        lines.append("“通过”只表示已经算出的数值阈值都满足，不包括尚未计算的项目，也不包括文字要求。")

    fail_counts = summary.get("fail_counts") or {}
    if fail_counts:
        ranked = sorted(fail_counts.items(), key=lambda item: item[1], reverse=True)
        lines.append("未通过次数最多的是：" + "；".join(f"{name}（{count}）" for name, count in ranked[:5]))

    pending = []
    for row in molecules:
        for item in row.get("pending") or []:
            if item not in pending:
                pending.append(item)
    if pending:
        lines.append("这些数值要求还没有数据，相关分子标为数据不足：" + "；".join(pending))
    if summary.get("text_requirements"):
        lines.append(
            "以下文字要求没有自动判断：" + "；".join(summary["text_requirements"])
        )
    if summary.get("not_applied"):
        lines.append("以下要求没有用于过滤：" + "；".join(summary["not_applied"]))

    alert_names = []
    for row in molecules:
        for name in row.get("alerts") or []:
            if name not in alert_names:
                alert_names.append(name)
    if alert_names:
        lines.append("部分结构命中了启发式子结构警示（" + "、".join(alert_names) + "）。这不是毒性结论。")

    next_steps = []
    if any(tool.get("status") == "skipped" for tool in tools):
        next_steps.append("有外部软件被跳过。安装程序或改正输入文件里的路径、受体和接口后再运行。")
    if any(tool.get("status") == "failed" for tool in tools):
        next_steps.append("有外部软件运行失败。查看该工具目录中的 stdout.txt、stderr.txt 或 response.json。")
    if pending:
        next_steps.append("如需自动判断 IC50 或结合能，在计算软件段写上本地程序或接口，并让结果里带 SMILES 列。")
    if summary.get("text_requirements"):
        next_steps.append("hERG、AMES 等文字要求需要预测接口或实验，不能用结构警示代替。")
    if summary.get("n_failed") and not summary.get("n_passed"):
        next_steps.append("没有分子通过已算出的数值阈值。可以放宽最常失败的条件，或更换分子。")
    if not spec.smiles:
        next_steps.append("在 [分子] 段写入 SMILES，或用分子文件给出结构。")
    if not summary.get("rdkit") and spec.smiles:
        next_steps.append("运行 setup.sh 或 setup.bat，装好带 RDKit 的 .venv 后再计算理化性质。")
    if not next_steps:
        next_steps.append("打开 properties.csv，优先看数值要求为“通过”、结构警示较少的分子。")

    reasons = [f"{name}：{count} 个分子未通过" for name, count in fail_counts.items()]
    return {"overview": "\n".join(lines), "failure_reasons": reasons, "next_steps": next_steps}


def render_report(context: dict) -> str:
    spec = context["spec"]
    summary = context["calc"]["summary"]
    molecules = context["calc"]["molecules"]
    analysis = context["analysis"]
    research = context["research"]
    rule = context["rule_findings"]
    lines = [
        "# 分子设计运行报告",
        "",
        f"- 运行编号: {context['run_id']}",
        f"- 模式: {context['mode']}",
        f"- 输入文件: {context['input_path']}",
        "",
        "数值以本报告中的表格和 `properties.csv` 为准。模型文字只作解释，不覆盖计算结果。",
        "",
        "## 1. 要求分析",
        "",
        f"- 靶点: {spec.target or '未填写'}",
        f"- 种属: {spec.species or '未填写'}",
        f"- 适应症: {spec.indication or '未填写'}",
        f"- 目标: {spec.goal or '未填写'}",
        f"- 分子数: {len(spec.smiles)}",
        "",
        "### 已解析的要求",
        "",
    ]
    if spec.constraints:
        for item in spec.constraints:
            lines.append(f"- {item.display()}")
    else:
        lines.append("- 没有解析到阈值或文字要求。")
    if spec.notes:
        lines.extend(["", "### 备注", ""])
        lines.extend(f"- {note}" for note in spec.notes)
    if spec.extras:
        lines.extend(["", "### 读取说明", ""])
        lines.extend(f"- {item}" for item in spec.extras)

    llm_analysis = analysis.get("requirement_analysis") or {}
    if llm_analysis.get("summary"):
        lines.extend(["", "### 模型对需求的理解", "", str(llm_analysis["summary"])])
    for title, key in (("假设", "assumptions"), ("信息缺口", "missing_information"), ("风险", "risks")):
        values = _as_list(llm_analysis.get(key))
        if values:
            lines.extend(["", f"### {title}", ""])
            lines.extend(f"- {item}" for item in values)
    if context["mode"].startswith("仅计算"):
        lines.extend(["", "未配置 LLM API Key，本段由规则解析生成，没有做额外推断。"])
    elif analysis.get("llm_errors"):
        lines.extend(["", "部分模型调用失败，失败步骤已改用规则结果："])
        lines.extend(f"- {item}" for item in analysis["llm_errors"])

    lines.extend(["", "## 2. 任务拆解", ""])
    for task in context["tasks"]:
        status = STATUS_LABEL.get(task["status"], task["status"])
        detail = task.get("detail") or ""
        lines.append(f"- {task['name']}：{status}" + (f"。{detail}" if detail else ""))
    suggested = _as_list(llm_analysis.get("suggested_checks"))
    if suggested:
        lines.extend(["", "模型建议核对以下事项，这些事项不会被自动执行：", ""])
        lines.extend(f"- {item}" for item in suggested)
    lines.append("")
    lines.append("实际执行的外部程序只来自输入文件中的计算软件段，不会采用模型临时写出的命令。")

    lines.extend(["", "## 3. 调研", "", research.get("summary") or "没有调研结果。"])
    llm_research = analysis.get("research_llm") or {}
    if llm_research.get("background"):
        lines.extend(["", "### 模型补充", "", str(llm_research["background"])])
    if llm_research.get("species_note"):
        lines.extend(["", str(llm_research["species_note"])])
    for title, key in (("设计提示", "design_implications"), ("不确定处", "caveats")):
        values = _as_list(llm_research.get(key))
        if values:
            lines.extend(["", f"### {title}", ""])
            lines.extend(f"- {item}" for item in values)

    lines.extend(["", "## 4. 计算", "", "### 内置理化与类药性", ""])
    for note in summary.get("notes") or []:
        lines.append(f"- {note}")
    if not summary.get("notes"):
        lines.append(
            f"- RDKit {'可用' if summary.get('rdkit') else '不可用'}；"
            f"真实 SA Score {'可用' if summary.get('real_sa') else '不可用'}。"
        )
    lines.append("- Lipinski 与 Veber 两列只供参考，除非你在输入里写了对应阈值，否则不参与通过与否。")
    lines.extend(["", _molecule_table(molecules, summary), ""])
    lines.extend(["### 外部软件与接口", ""])
    if not context["tools"]:
        lines.append("输入文件没有启用外部计算软件。需要对接或预测接口时，在 `[计算软件]` 段写上本地路径或 API 地址。")
    else:
        for tool in context["tools"]:
            status = STATUS_LABEL.get(tool["status"], tool["status"])
            lines.append(f"- {tool['name']}（{tool['type']}）：{status}")
            if tool.get("reason"):
                lines.append(f"  - {tool['reason']}")
            if tool.get("command"):
                lines.append("  - 命令: `" + " ".join(str(part) for part in tool["command"]) + "`")
            lines.append(f"  - 解析到的结果行: {tool.get('parsed_rows', 0)}")
            if tool["status"] == "failed" and tool.get("stderr_tail"):
                lines.append("  - 错误输出: " + " ".join(tool["stderr_tail"].split()))

    lines.extend(["", "## 5. 计算数据分析", ""])
    llm_findings = analysis.get("findings_llm") or {}
    if llm_findings.get("overview"):
        lines.extend([str(llm_findings["overview"]), "", "以下统计来自计算表：", ""])
    lines.append(rule["overview"])
    if rule["failure_reasons"]:
        lines.extend(["", "### 未通过项", ""])
        lines.extend(f"- {item}" for item in rule["failure_reasons"])
    extra_reasons = _as_list(llm_findings.get("failure_reasons"))
    if extra_reasons:
        lines.extend(["", "### 模型指出的原因", ""])
        lines.extend(f"- {item}" for item in extra_reasons)

    lines.extend(["", "### 建议的下一步", ""])
    steps = list(rule["next_steps"])
    for item in _as_list(llm_findings.get("next_steps")):
        if item not in steps:
            steps.append(item)
    lines.extend(f"- {item}" for item in steps)

    lines.extend(
        [
            "",
            "## 6. 报告整理",
            "",
            "本次运行写了下面这些文件：",
            "",
            f"- 报告: `{context['report_name']}`",
            f"- 性质表: `{context['properties_name']}`",
            "- 结构化需求: `spec.json`",
            "- 调研: `research.json`",
            "- 任务: `tasks.json`",
            "- 计算与分析: `calculations.json`、`analysis.json`",
            "- 日志: `run.log`",
            "",
            "外部工具的原始输出在各自的子目录里。重新运行会新建一个带时间戳的目录，不会覆盖这一次的结果。",
            "",
        ]
    )
    return "\n".join(lines)


def write_properties(path: Path, molecules: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow([label for _, label in _CSV_COLUMNS])
        for row in molecules:
            writer.writerow([_csv_value(row, key) for key, _ in _CSV_COLUMNS])


def _molecule_table(molecules: list[dict], summary: dict) -> str:
    if not molecules:
        return "没有分子计算结果。"
    sa_key = "SA" if summary.get("real_sa") else "SAProxy"
    sa_label = "SA" if summary.get("real_sa") else "SA代理"
    header = ["输入SMILES", "数值要求", "分子量", "LogP", "QED", sa_label, "未通过项", "缺少数据"]
    rows = [
        "| " + " | ".join(header) + " |",
        "| " + " | ".join("---" for _ in header) + " |",
    ]
    shown = molecules[:30]
    for row in shown:
        rows.append(
            "| "
            + " | ".join(
                [
                    _cell(row.get("input_smiles")),
                    STATUS_LABEL.get(row.get("numeric_status"), ""),
                    _fmt(row.get("MW")),
                    _fmt(row.get("LogP")),
                    _fmt(row.get("QED")),
                    _fmt(row.get(sa_key)),
                    _cell("；".join(row.get("failed") or [])),
                    _cell("；".join(row.get("pending") or [])),
                ]
            )
            + " |"
        )
    if len(molecules) > len(shown):
        rows.append("")
        rows.append(f"表中只显示前 {len(shown)} 行，全部 {len(molecules)} 行在 properties.csv。")
    return "\n".join(rows)


def _csv_value(row: dict, key: str):
    value = row.get(key)
    if key == "numeric_status":
        return STATUS_LABEL.get(value, value or "")
    if key in {"failed", "pending", "alerts"}:
        return "；".join(value or [])
    if isinstance(value, bool):
        return "是" if value else "否"
    if value is None:
        return ""
    if isinstance(value, float):
        return _fmt(value)
    return value


def _fmt(value) -> str:
    if value is None or value == "":
        return ""
    number = float(value)
    if abs(number - round(number)) < 1e-8 and abs(number) >= 1:
        return str(int(round(number)))
    return f"{number:.2f}"


def _cell(value) -> str:
    return str(value or "").replace("|", "\\|").replace("\n", " ")


def _as_list(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        text = value.strip()
        return [text] if text else []
    if isinstance(value, list):
        return [str(item).strip() for item in value if str(item).strip()]
    return [str(value).strip()]
