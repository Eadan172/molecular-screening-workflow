#!/usr/bin/env bash
# 使用项目目录里的 .venv 运行。还没有环境时会先执行 setup.sh。
set -euo pipefail
cd "$(dirname "$0")"
if [[ ! -x .venv/bin/python ]]; then
  bash setup.sh
fi
exec .venv/bin/python run.py "$@"
