"""读取可选的大模型配置。不把密钥写进 os.environ。"""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


DEFAULT_BASE_URL = "https://api.openai.com/v1"
DEFAULT_MODEL = "gpt-4o-mini"

_PLACEHOLDERS = {
    "",
    "your_api_key",
    "your-api-key",
    "your_key",
    "sk-xxx",
    "sk-...",
    "changeme",
    "none",
    "null",
    "todo",
    "xxxx",
    "请填写",
    "请填写你的key",
    "请填写你的密钥",
}


@dataclass
class Settings:
    api_key: str = ""
    base_url: str = DEFAULT_BASE_URL
    model: str = DEFAULT_MODEL
    calc_only: bool = False

    @property
    def llm_enabled(self) -> bool:
        if self.calc_only:
            return False
        key = (self.api_key or "").strip()
        if key.lower() in _PLACEHOLDERS:
            return False
        if key.startswith("请填写"):
            return False
        return True


def read_env_file(path: Path) -> dict[str, str]:
    data: dict[str, str] = {}
    if not path.is_file():
        return data
    for raw in path.read_text(encoding="utf-8-sig").splitlines():
        line = raw.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        data[key.strip()] = value.strip().strip('"').strip("'")
    return data


def load_settings(root: Path, calc_only: bool = False) -> Settings:
    file_values = read_env_file(Path(root) / ".env")

    def pick(name: str, default: str = "") -> str:
        env_value = os.environ.get(name)
        if env_value is not None and env_value.strip():
            return env_value.strip()
        file_value = file_values.get(name)
        if file_value:
            return file_value.strip()
        return default

    return Settings(
        api_key=pick("LLM_API_KEY"),
        base_url=pick("LLM_BASE_URL", DEFAULT_BASE_URL) or DEFAULT_BASE_URL,
        model=pick("LLM_MODEL", DEFAULT_MODEL) or DEFAULT_MODEL,
        calc_only=calc_only,
    )
