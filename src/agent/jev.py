"""Jev 决策层：只做语义门控和有界路由，不替代化学计算。"""

from __future__ import annotations

from src.agent.settings import Settings


TREND_OPTIONS = {
    "conventional_small_molecule": "常规小分子：已有可成药口袋，以活性、选择性、ADMET 和合成性推进。",
    "structure_based_3d": "结构驱动 3D 设计：有可靠口袋/复合物结构，适合对接、构象、MD 或口袋条件生成。",
    "covalent_allosteric": "共价或变构设计：存在可验证的反应性位点或变构口袋，需额外选择性与安全性验证。",
    "targeted_protein_degradation": "靶向蛋白降解/分子胶：抑制不足以覆盖生物学机制，且具备 E3、三元复合物或降解证据。",
    "rna_targeting": "RNA 靶向：疾病机制和可操作结构主要位于 RNA，需专门的 RNA 结构、选择性和递送证据。",
    "insufficient_evidence": "当前证据不足以选择模态，应先补齐靶点机制、结构或实验信息。",
}

EVIDENCE_OPTIONS = {
    "sufficient": "现有状态足以支持当前阶段的靶点/种属识别和计算路线；不等于支持疗效结论。",
    "insufficient": "关键信息缺失，现有状态不足以支持项目路线。",
    "contradictory": "输入要求与检索结果之间存在实质冲突，需要人工核对。",
}

ACTION_OPTIONS = {
    "proceed": "继续当前计算路线，同时保留正常人工复核。",
    "add_structure_affinity": "先补充可靠蛋白结构、结合位点、对接/MD 或活性数据。",
    "add_admet_synthesis": "先补充 ADMET、不确定性、适用域、合成路线或原料可得性。",
    "human_review": "暂停自动推进，由药化/生物学专家复核靶点、种属、模态或冲突证据。",
}

RISK_TIERS = [
    "低：信息和计算覆盖较完整，只有常规项目风险。",
    "中：存在一项重要数据缺口，但可通过明确的计算或实验补齐。",
    "高：存在多个关键缺口、超出适用域或模态选择不明确。",
    "很高：靶点/种属/证据冲突，或自动推进可能传播错误。",
]


class JevError(RuntimeError):
    pass


class JevClient:
    def __init__(self, settings: Settings, http_request=None):
        self.settings = settings
        self._http_request = http_request

    @property
    def enabled(self) -> bool:
        return self.settings.jev_enabled

    def decide(self, state: dict) -> dict:
        if not self.enabled:
            raise JevError("未配置 JEV_API_KEY")
        payload = {
            "model": self.settings.jev_model,
            "state": state,
            "questions": {
                "evidence_gate": {
                    "type": "choice",
                    "instructions": (
                        "只根据 state 判断当前证据是否足以支持靶点/种属识别和当前阶段路线。"
                        "不要把数据库背景当成本次分子的活性证据。"
                    ),
                    "criteria": EVIDENCE_OPTIONS,
                },
                "project_modality": {
                    "type": "choice",
                    "instructions": (
                        "选择最适合优先评估的药物发现模态。只有 state 中已有证据时才选择新模态；"
                        "不确定时选 insufficient_evidence。"
                    ),
                    "criteria": TREND_OPTIONS,
                },
                "next_action": {
                    "type": "choice",
                    "instructions": (
                        "选择风险最低、信息价值最高的下一步。缺少活性/结构/ADMET 时不得直接选择 proceed。"
                    ),
                    "criteria": ACTION_OPTIONS,
                },
                "decision_risk": {
                    "type": "score",
                    "instructions": "评估按当前状态继续项目决策的风险，低到高。",
                    "criteria": RISK_TIERS,
                },
                "needs_human_review": {
                    "type": "noul",
                    "instructions": "当前状态是否需要领域专家人工复核后才能继续自动推进？",
                },
            },
        }
        request = self._http_request or _default_request
        response = request(
            "POST",
            decision_url(self.settings.jev_base_url),
            json=payload,
            headers={
                "Authorization": f"Bearer {self.settings.jev_api_key}",
                "Content-Type": "application/json",
            },
            timeout=90,
        )
        status = getattr(response, "status_code", 0)
        text = getattr(response, "text", "") or ""
        if status >= 400:
            raise JevError(f"Jev 接口返回 {status}: {text[:180]}")
        try:
            data = response.json()
        except (TypeError, ValueError) as exc:
            raise JevError(f"Jev 接口不是 JSON: {text[:180]}") from exc
        if not isinstance(data, dict) or not isinstance(data.get("answers"), dict):
            raise JevError("Jev 响应里没有 answers")
        return summarize_decision(data)


def decision_url(base_url: str) -> str:
    base = (base_url or "").strip().rstrip("/")
    if not base:
        base = "https://jev-ai.org/api/v1"
    if base.endswith("/systemone"):
        return base + "/"
    return base + "/systemone/"


def summarize_decision(data: dict) -> dict:
    answers = data.get("answers") or {}
    result = {
        "status": "done",
        "model": data.get("model") or "",
        "model_version": data.get("model_version") or "",
        "request_id": data.get("id") or "",
        "latency_ms": data.get("latency_ms"),
        "usage": data.get("usage") or {},
        "answers": {},
        "policy": {},
    }
    for question_id, answer in answers.items():
        if not isinstance(answer, dict):
            continue
        kind = answer.get("type") or ""
        compact = {"type": kind}
        if kind == "choice":
            compact.update(
                {
                    "choice": answer.get("choice") or "",
                    "confidence": _number(answer.get("confidence")),
                    "probabilities": answer.get("probabilities") or {},
                }
            )
        elif kind == "score":
            compact.update(
                {
                    "score": _number(answer.get("score")),
                    "confidence": _number(answer.get("confidence")),
                    "probabilities": answer.get("probabilities") or {},
                    "legend": answer.get("legend") or {},
                }
            )
        elif kind == "noul":
            compact["noul"] = _number(answer.get("noul"))
        result["answers"][question_id] = compact

    result["policy"] = apply_policy(result["answers"])
    return result


def apply_policy(answers: dict) -> dict:
    evidence = answers.get("evidence_gate") or {}
    action = answers.get("next_action") or {}
    review = answers.get("needs_human_review") or {}
    evidence_choice = evidence.get("choice") or "unknown"
    action_choice = action.get("choice") or "human_review"
    evidence_confidence = evidence.get("confidence")
    action_confidence = action.get("confidence")
    review_probability = review.get("noul")

    reasons = []
    route = action_choice
    if evidence_choice == "contradictory":
        route = "human_review"
        reasons.append("证据门控判为 contradictory。")
    elif evidence_choice != "sufficient":
        route = "add_structure_affinity"
        reasons.append("证据门控没有达到 sufficient。")
    if evidence_confidence is None or evidence_confidence < 0.85:
        route = "human_review"
        reasons.append("证据门控置信度低于 0.85。")
    if action_confidence is None or action_confidence < 0.75:
        route = "human_review"
        reasons.append("下一步路由置信度低于 0.75。")
    if review_probability is not None and review_probability >= 0.70:
        route = "human_review"
        reasons.append("人工复核概率达到 0.70。")
    return {
        "route": route,
        "allow_automatic_progress": route == "proceed",
        "thresholds": {
            "evidence_confidence": 0.85,
            "action_confidence": 0.75,
            "human_review_probability": 0.70,
        },
        "reasons": reasons,
    }


def trend_label(choice: str) -> str:
    return TREND_OPTIONS.get(choice, choice or "未返回")


def action_label(choice: str) -> str:
    return ACTION_OPTIONS.get(choice, choice or "未返回")


def evidence_label(choice: str) -> str:
    return EVIDENCE_OPTIONS.get(choice, choice or "未返回")


def _default_request(method, url, json=None, headers=None, timeout=None):
    import requests

    return requests.request(method, url, json=json, headers=headers, timeout=timeout)


def _number(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
