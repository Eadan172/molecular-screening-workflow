#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""读取需求文本并运行分子设计流程。"""

from __future__ import annotations

import argparse
import sys
import traceback
from pathlib import Path


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="读取分子需求文本，完成要求分析、任务拆解、调研、计算、数据分析和报告。"
    )
    parser.add_argument(
        "input",
        nargs="?",
        default="input/request.txt",
        help="需求文本路径，默认 input/request.txt",
    )
    parser.add_argument(
        "--calc-only",
        action="store_true",
        help="即使填写了 LLM/Jev API Key，也只运行计算和规则报告",
    )
    parser.add_argument(
        "--output",
        default=None,
        help="结果根目录，默认 <项目>/results/runs",
    )
    args = parser.parse_args(argv)

    try:
        from src.agent.pipeline import run_pipeline
        from src.agent.settings import load_settings
    except ImportError as exc:
        print(f"无法加载流程: {exc}")
        print("请先在项目目录运行 bash setup.sh 或 setup.bat。")
        return 1

    request_path = _resolve_input(args.input, ROOT)
    if request_path is None:
        example = ROOT / "input" / "request.txt"
        print(f"找不到输入文件: {args.input}")
        print(f"可以编辑示例文件后重试: {example}")
        return 1

    settings = load_settings(ROOT, calc_only=args.calc_only)
    output_root = Path(args.output) if args.output else ROOT / "results" / "runs"
    if not output_root.is_absolute():
        output_root = ROOT / output_root

    print(f"Python: {sys.executable}")
    if settings.llm_enabled:
        print(f"已检测到 LLM API Key，将调用 {settings.model} 做分析、调研补充和数据解释。")
    else:
        print("未检测到 LLM API Key，直接运行计算，并用规则生成报告。")
    if settings.jev_enabled:
        print(f"已检测到 Jev API Key，将调用 {settings.jev_model} 做证据门控和下一步路由。")
    else:
        print("未检测到 Jev API Key，跳过 Jev 决策门控。")

    try:
        result = run_pipeline(
            request_path=request_path,
            settings=settings,
            output_root=output_root,
            project_root=ROOT,
        )
    except FileNotFoundError as exc:
        print(exc)
        return 1
    except Exception:
        traceback.print_exc()
        return 1

    print(f"模式: {result['mode']}")
    print(f"报告: {result['report']}")
    print(f"性质表: {result['properties']}")
    return 0


def _resolve_input(raw: str, root: Path) -> Path | None:
    candidate = Path(raw)
    if candidate.is_file():
        return candidate.resolve()
    rooted = root / raw
    if rooted.is_file():
        return rooted.resolve()
    return None


if __name__ == "__main__":
    sys.exit(main())
