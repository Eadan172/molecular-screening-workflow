#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""命中之后的可选结构工具。不改筛选、对接或模型预测。"""

from __future__ import annotations

import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.leadopt.cli import main


if __name__ == "__main__":
    raise SystemExit(main())
