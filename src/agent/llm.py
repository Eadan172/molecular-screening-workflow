"""OpenAI 兼容的大模型客户端。没有密钥时不会被调用。"""

from __future__ import annotations

import json
import re

from src.agent.settings import Settings


class LLMError(RuntimeError):
    pass


_PROMPTS = {
    "analyze": (
        "你是药物化学助手。根据用户的分子设计要求和已经解析出的字段，"
        "输出一个 JSON 对象，键为 summary、assumptions、missing_information、risks、suggested_checks。"
        "summary 用中文，80到150字。其余四个键是字符串数组。"
        "不要改写用户给出的数值阈值，不要编造实验数据或文献编号。只输出 JSON。"
    ),
    "research": (
        "你是药物化学助手。只根据给定的用户要求和公开数据库摘要，"
        "输出一个 JSON 对象，键为 background、species_note、design_implications、caveats。"
        "background 和 species_note 是字符串，后两个键是字符串数组。"
        "不要编造 ChEMBL 编号、文献或活性数值。数据库摘要里没有的内容就写进 caveats。只输出 JSON。"
    ),
    "findings": (
        "你是药物化学助手。只根据给定的计算统计解释结果，"
        "输出一个 JSON 对象，键为 overview、failure_reasons、next_steps。"
        "overview 是中文字符串，后两个键是字符串数组。"
        "不要编造表里没有的分子、分数或实验结论。只输出 JSON。"
    ),
}


def chat_url(base_url: str) -> str:
    base = (base_url or "").strip().rstrip("/")
    if not base:
        base = "https://api.openai.com/v1"
    if base.endswith("/chat/completions"):
        return base
    return base + "/chat/completions"


def extract_json(text: str) -> dict:
    raw = (text or "").strip()
    fenced = re.search(r"```(?:json)?\s*(.*?)```", raw, flags=re.DOTALL | re.IGNORECASE)
    if fenced:
        raw = fenced.group(1).strip()
    start = raw.find("{")
    end = raw.rfind("}")
    if start < 0 or end <= start:
        raise LLMError("模型没有返回 JSON 对象")
    data = json.loads(raw[start : end + 1])
    if not isinstance(data, dict):
        raise LLMError("模型返回的 JSON 不是对象")
    return data


class LLMClient:
    def __init__(self, settings: Settings, http_request=None):
        self.settings = settings
        self._http_request = http_request

    @property
    def enabled(self) -> bool:
        return self.settings.llm_enabled

    def complete_json(self, task: str, payload: dict) -> dict:
        if not self.enabled:
            raise LLMError("未配置 LLM API Key")
        prompt = _PROMPTS.get(task)
        if prompt is None:
            raise LLMError(f"未知模型任务: {task}")
        body = {
            "model": self.settings.model,
            "temperature": 0.2,
            "messages": [
                {"role": "system", "content": prompt},
                {"role": "user", "content": json.dumps(payload, ensure_ascii=False)},
            ],
        }
        request = self._http_request or _default_request
        response = request(
            "POST",
            chat_url(self.settings.base_url),
            json=body,
            headers={
                "Authorization": f"Bearer {self.settings.api_key}",
                "Content-Type": "application/json",
            },
            timeout=90,
        )
        status = getattr(response, "status_code", 0)
        text = getattr(response, "text", "") or ""
        if status >= 400:
            raise LLMError(f"模型接口返回 {status}: {text[:180]}")
        try:
            data = response.json() if hasattr(response, "json") else json.loads(text)
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise LLMError(f"模型接口不是 JSON: {text[:180]}") from exc
        content = _message_text(data)
        return extract_json(content)


def _default_request(method, url, json=None, headers=None, timeout=None):
    import requests

    return requests.request(method, url, json=json, headers=headers, timeout=timeout)


def _message_text(data: dict) -> str:
    try:
        content = data["choices"][0]["message"]["content"]
    except (KeyError, IndexError, TypeError) as exc:
        raise LLMError("模型响应里没有 message.content") from exc
    if isinstance(content, str):
        return content
    if isinstance(content, list):
        parts = []
        for item in content:
            if isinstance(item, str):
                parts.append(item)
            elif isinstance(item, dict):
                parts.append(str(item.get("text") or ""))
        return "".join(parts)
    raise LLMError("无法读取模型文本")
