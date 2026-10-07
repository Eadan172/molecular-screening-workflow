"""流程里流转的数据结构。"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class Constraint:
    """用户写下的一条要求。kind 为 numeric 或 text。"""

    name: str
    kind: str
    op: str
    value: float | None = None
    high: float | None = None
    unit: str = ""
    raw: str = ""
    source_key: str = ""

    def display(self) -> str:
        label = self.source_key or self.name
        if self.kind == "numeric" and self.op == "between":
            return f"{label}: {self.value:g}–{self.high:g} {self.unit}".strip()
        if self.kind == "numeric" and self.value is not None:
            return f"{label}: {self.op} {self.value:g} {self.unit}".strip()
        return f"{label}: {self.raw}".strip()


@dataclass
class ToolSpec:
    """输入文件里声明的本地程序或 HTTP 接口。"""

    name: str
    type: str = "local"
    path: str = ""
    url: str = ""
    method: str = "POST"
    args: str = ""
    receptor: str = ""
    timeout: int = 300
    header: str = ""


@dataclass
class RequestSpec:
    """从需求文本解析出的结构化要求。"""

    raw_text: str = ""
    target: str = ""
    species: str = ""
    indication: str = ""
    goal: str = ""
    notes: list[str] = field(default_factory=list)
    smiles: list[str] = field(default_factory=list)
    molecule_file: str = ""
    constraints: list[Constraint] = field(default_factory=list)
    tools: list[ToolSpec] = field(default_factory=list)
    extras: list[str] = field(default_factory=list)

    def public_dict(self) -> dict:
        """去掉原文和请求头，避免把密钥写进结果目录。"""
        return {
            "target": self.target,
            "species": self.species,
            "indication": self.indication,
            "goal": self.goal,
            "notes": list(self.notes),
            "smiles": list(self.smiles),
            "molecule_file": self.molecule_file,
            "constraints": [
                {
                    "name": item.name,
                    "kind": item.kind,
                    "op": item.op,
                    "value": item.value,
                    "high": item.high,
                    "unit": item.unit,
                    "raw": item.raw,
                    "source_key": item.source_key,
                    "display": item.display(),
                }
                for item in self.constraints
            ],
            "tools": [
                {
                    "name": tool.name,
                    "type": tool.type,
                    "path": tool.path,
                    "url": tool.url,
                    "method": tool.method,
                    "args": tool.args,
                    "receptor": tool.receptor,
                    "timeout": tool.timeout,
                    "has_header": bool(tool.header),
                }
                for tool in self.tools
            ],
            "extras": list(self.extras),
        }
