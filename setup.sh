#!/usr/bin/env bash
# 在项目目录创建 .venv 并安装运行依赖。
set -euo pipefail
cd "$(dirname "$0")"

if command -v python3 >/dev/null 2>&1; then
  PY=python3
elif command -v python >/dev/null 2>&1; then
  PY=python
else
  echo "未找到 Python。请先安装 Python 3.9 或更高版本。" >&2
  exit 1
fi

"$PY" - <<'PY'
import sys
if sys.version_info < (3, 9):
    raise SystemExit("需要 Python 3.9 或更高版本，当前是 " + sys.version.split()[0])
PY

if ! "$PY" -m venv .venv; then
  echo "创建 .venv 失败。Debian/Ubuntu 请先执行: sudo apt install python3-venv" >&2
  exit 1
fi
.venv/bin/python -m pip install --upgrade pip
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python - <<'PY'
import rdkit
import requests
print("rdkit", rdkit.__version__)
print("requests", requests.__version__)
PY

if [[ ! -f .env ]]; then
  cp .env.example .env
  echo "已创建 .env。"
fi

if [[ -t 0 ]]; then
  current="$(.venv/bin/python - <<'PY'
from pathlib import Path
for line in Path(".env").read_text(encoding="utf-8").splitlines():
    if line.startswith("LLM_API_KEY="):
        print(line.split("=", 1)[1].strip())
        break
PY
)"
  if [[ -z "${current}" ]]; then
    echo "如需大模型分析，请输入 LLM API Key。直接回车则只运行计算。"
    read -r -p "LLM_API_KEY: " key || true
    if [[ -n "${key:-}" ]]; then
      SETUP_LLM_KEY="$key" .venv/bin/python - <<'PY'
import os
from pathlib import Path
key = os.environ["SETUP_LLM_KEY"]
path = Path(".env")
lines = []
found = False
for line in path.read_text(encoding="utf-8").splitlines():
    if line.startswith("LLM_API_KEY="):
        lines.append("LLM_API_KEY=" + key)
        found = True
    else:
        lines.append(line)
if not found:
    lines.append("LLM_API_KEY=" + key)
path.write_text("\n".join(lines) + "\n", encoding="utf-8")
PY
      echo "已写入 .env。"
    fi
  fi
fi

echo
echo "环境已就绪: $(pwd)/.venv"
echo "下一步: 编辑 input/request.txt ，然后运行 ./run.sh"
