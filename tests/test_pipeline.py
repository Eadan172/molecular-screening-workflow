import json
import os
import sys
import tempfile
import unittest
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.agent.llm import LLMClient, chat_url, extract_json
from src.agent.parser import attach_molecule_file, parse_request
from src.agent.pipeline import run_pipeline
from src.agent.research import lookup_target
from src.agent.settings import Settings, load_settings


def _no_network(url, params, timeout):
    return {"targets": []}


def _tree_text(root: Path) -> str:
    parts = []
    for path in root.rglob("*"):
        if path.is_file():
            parts.append(path.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(parts)


class ParserTests(unittest.TestCase):
    def test_sections_constraints_and_tools(self):
        spec = parse_request(
            """
            [基本信息]
            靶点: AKT1
            种属: Homo sapiens
            适应症: 肿瘤
            目标: 口服抑制剂
            LogP: 1-3
            [分子]
            CCO aspirin
            c1ccccc1
            [活性]
            IC50: < 100 nM
            hERG: 避免
            [计算软件]
            name: remote
            url: https://example.com/predict
            header: Authorization: Bearer SUPER_SECRET_TOKEN
            name: local_tool
            type: local
            path: /opt/vina
            args: --ligand {ligand}
            """
        )
        self.assertEqual(spec.target, "AKT1")
        self.assertEqual(spec.species, "Homo sapiens")
        self.assertEqual(spec.smiles, ["CCO", "c1ccccc1"])
        logp = next(item for item in spec.constraints if item.name == "LogP")
        self.assertEqual(logp.op, "between")
        self.assertEqual(logp.value, 1)
        self.assertEqual(logp.high, 3)
        ic50 = next(item for item in spec.constraints if item.name == "IC50")
        self.assertEqual(ic50.op, "<")
        self.assertEqual(ic50.value, 100)
        self.assertEqual(ic50.unit, "nM")
        herg = next(item for item in spec.constraints if item.name == "hERG")
        self.assertEqual(herg.op, "avoid")
        self.assertEqual(spec.tools[0].type, "api")
        self.assertEqual(spec.tools[1].type, "local")
        public = json.dumps(spec.public_dict(), ensure_ascii=False)
        self.assertNotIn("SUPER_SECRET_TOKEN", public)
        self.assertTrue(spec.public_dict()["tools"][0]["has_header"])

    def test_example_file_and_molecule_file(self):
        example = (ROOT / "input" / "request.txt").read_text(encoding="utf-8")
        spec = parse_request(example)
        self.assertEqual(spec.target, "AKT1")
        self.assertEqual(len(spec.smiles), 3)
        self.assertEqual(spec.tools, [])
        self.assertTrue(any("晶体" in note for note in spec.notes))

        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            (folder / "mols.csv").write_text("SMILES,name\nCCO,ethanol\n", encoding="utf-8")
            request = folder / "request.txt"
            request.write_text("[基本信息]\n分子文件: mols.csv\n[分子]\nCCC\n", encoding="utf-8")
            loaded = parse_request(request.read_text(encoding="utf-8"))
            attach_molecule_file(loaded, [folder])
            self.assertEqual(loaded.smiles, ["CCC", "CCO"])


class SettingsTests(unittest.TestCase):
    def setUp(self):
        self.saved = {name: os.environ.get(name) for name in ("LLM_API_KEY", "LLM_BASE_URL", "LLM_MODEL")}
        for name in self.saved:
            os.environ.pop(name, None)

    def tearDown(self):
        for name, value in self.saved.items():
            if value is None:
                os.environ.pop(name, None)
            else:
                os.environ[name] = value

    def test_placeholder_and_override(self):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            (root / ".env").write_text("LLM_API_KEY=your_api_key\nLLM_MODEL=file-model\n", encoding="utf-8")
            self.assertFalse(load_settings(root).llm_enabled)
            (root / ".env").write_text(
                "LLM_API_KEY=filekey\nLLM_BASE_URL=https://file.example/v1\nLLM_MODEL=file-model\n",
                encoding="utf-8",
            )
            settings = load_settings(root)
            self.assertTrue(settings.llm_enabled)
            self.assertEqual(settings.model, "file-model")
            self.assertFalse(load_settings(root, calc_only=True).llm_enabled)
            os.environ["LLM_API_KEY"] = "envkey"
            self.assertEqual(load_settings(root).api_key, "envkey")


class ResearchTests(unittest.TestCase):
    def test_prefers_matching_species_and_survives_errors(self):
        spec = parse_request("[基本信息]\n靶点: AKT1\n种属: Homo sapiens\n")

        def http_get(url, params, timeout):
            if "activity" in url:
                return {
                    "activities": [
                        {
                            "molecule_chembl_id": "CHEMBL999",
                            "canonical_smiles": "CCO",
                            "standard_type": "IC50",
                            "standard_value": "10",
                            "standard_units": "nM",
                        }
                    ]
                }
            return {
                "targets": [
                    {
                        "target_chembl_id": "CHEMBL1",
                        "pref_name": "AKT1 mouse",
                        "organism": "Mus musculus",
                        "target_type": "SINGLE PROTEIN",
                    },
                    {
                        "target_chembl_id": "CHEMBL2",
                        "pref_name": "AKT1 human",
                        "organism": "Homo sapiens",
                        "target_type": "SINGLE PROTEIN",
                    },
                ]
            }

        found = lookup_target(spec, http_get=http_get)
        self.assertEqual(found["chosen"]["chembl_id"], "CHEMBL2")
        self.assertIn("一致", found["summary"])

        def substrate_first(url, params, timeout):
            if "activity" in url:
                return {"activities": []}
            return {
                "targets": [
                    {
                        "target_chembl_id": "CHEMBL1255161",
                        "pref_name": "Proline-rich AKT1 substrate 1",
                        "organism": "Homo sapiens",
                        "target_type": "SINGLE PROTEIN",
                    },
                    {
                        "target_chembl_id": "CHEMBL4282",
                        "pref_name": "RAC-alpha serine/threonine-protein kinase",
                        "organism": "Homo sapiens",
                        "target_type": "SINGLE PROTEIN",
                    },
                ]
            }

        kinase = lookup_target(spec, http_get=substrate_first)
        self.assertEqual(kinase["chosen"]["chembl_id"], "CHEMBL4282")

        mouse = parse_request("[基本信息]\n靶点: AKT1\n种属: 小鼠\n")

        def only_human(url, params, timeout):
            if "activity" in url:
                raise RuntimeError("down")
            return {
                "targets": [
                    {
                        "target_chembl_id": "CHEMBL2",
                        "pref_name": "AKT1",
                        "organism": "Homo sapiens",
                        "target_type": "SINGLE PROTEIN",
                    }
                ]
            }

        mismatch = lookup_target(mouse, http_get=only_human)
        self.assertIn("不一致", mismatch["summary"])

        def broken(url, params, timeout):
            raise RuntimeError("offline")

        failed = lookup_target(spec, http_get=broken)
        self.assertIn("跳过", failed["summary"])
        empty = lookup_target(parse_request("目标: 只写了目标\n"), http_get=broken)
        self.assertIn("没有靶点", empty["summary"])


class LlmUnitTests(unittest.TestCase):
    def test_json_and_request_shape(self):
        self.assertEqual(chat_url("https://example.com/v1/"), "https://example.com/v1/chat/completions")
        self.assertEqual(
            chat_url("https://example.com/v1/chat/completions"),
            "https://example.com/v1/chat/completions",
        )
        self.assertEqual(extract_json('说明\n```json\n{"summary": "ok"}\n```'), {"summary": "ok"})

        captured = {}

        def http_request(method, url, json=None, headers=None, timeout=None):
            captured["method"] = method
            captured["url"] = url
            captured["json"] = json
            captured["headers"] = headers

            class Response:
                status_code = 200
                text = ""

                def json(self):
                    return {"choices": [{"message": {"content": '{"summary": "模型摘要"}'}}]}

            return Response()

        client = LLMClient(
            Settings(api_key="sk-SUPERKEY", base_url="https://example.com/v1", model="demo"),
            http_request=http_request,
        )
        self.assertEqual(client.complete_json("analyze", {"target": "AKT1"})["summary"], "模型摘要")
        self.assertEqual(captured["headers"]["Authorization"], "Bearer sk-SUPERKEY")
        self.assertNotIn("sk-SUPERKEY", json.dumps(captured["json"]))


class PipelineTests(unittest.TestCase):
    def _run(self, text, settings, http_request=None, llm=None):
        tmp = tempfile.TemporaryDirectory()
        self.addCleanup(tmp.cleanup)
        folder = Path(tmp.name)
        request = folder / "request.txt"
        request.write_text(text, encoding="utf-8")
        result = run_pipeline(
            request_path=request,
            settings=settings,
            output_root=folder / "out",
            project_root=folder,
            llm=llm,
            http_get=_no_network,
            http_request=http_request,
            run_id="case",
        )
        return result

    def test_calc_only_skips_llm_and_records_missing_tool(self):
        class Boom:
            enabled = True

            def complete_json(self, task, payload):
                raise AssertionError("不应该调用大模型")

        text = f"""
        [基本信息]
        靶点: AKT1
        [分子]
        CCO
        [计算软件]
        name: echo
        type: local
        path: {sys.executable}
        args: -c print(123)
        name: missing
        type: local
        path: /this/does/not/exist/vina
        """
        result = self._run(text, Settings(api_key="sk-SHOULD-NOT-LEAK", calc_only=True), llm=Boom())
        report = Path(result["report"]).read_text(encoding="utf-8")
        for title in ("要求分析", "任务拆解", "调研", "计算", "计算数据分析", "报告整理"):
            self.assertIn(title, report)
        self.assertIn("仅计算", report)
        self.assertIn("找不到可执行文件", report)
        blob = _tree_text(Path(result["run_dir"]))
        self.assertNotIn("sk-SHOULD-NOT-LEAK", blob)
        data = json.loads(Path(result["run_dir"], "calculations.json").read_text(encoding="utf-8"))
        echo = next(item for item in data["tools"] if item["name"] == "echo")
        missing = next(item for item in data["tools"] if item["name"] == "missing")
        self.assertEqual(echo["status"], "done")
        self.assertIn("123", echo["stdout_tail"])
        self.assertEqual(missing["status"], "skipped")

    def test_api_result_merges_and_header_stays_out_of_report(self):
        secret = "SUPER_SECRET_TOKEN"
        seen = {}

        def http_request(method, url, json=None, headers=None, timeout=None):
            seen["headers"] = headers
            seen["json"] = json

            class Response:
                status_code = 200
                content = b'{"molecules":[{"smiles":"CCO","ic50":12}]}'

            return Response()

        text = f"""
        [分子]
        CCO
        [活性]
        IC50: < 100 nM
        [计算软件]
        name: score
        type: api
        url: http://127.0.0.1:9/predict
        header: Authorization: Bearer {secret}
        """
        result = self._run(text, Settings(), http_request=http_request)
        self.assertEqual(seen["headers"]["Authorization"], f"Bearer {secret}")
        self.assertNotIn(secret, json.dumps(seen["json"]))
        self.assertNotIn(secret, _tree_text(Path(result["run_dir"])))
        data = json.loads(Path(result["run_dir"], "calculations.json").read_text(encoding="utf-8"))
        row = data["molecules"][0]
        self.assertEqual(row["IC50"], 12)
        if row["valid"]:
            self.assertEqual(row["numeric_status"], "pass")

    def test_llm_text_is_added_and_errors_fall_back(self):
        class Fake:
            enabled = True

            def __init__(self):
                self.tasks = []

            def complete_json(self, task, payload):
                self.tasks.append(task)
                if task == "analyze":
                    return {
                        "summary": "模型摘要XYZ",
                        "assumptions": ["假设A"],
                        "missing_information": [],
                        "risks": [],
                        "suggested_checks": ["核对晶体结构"],
                    }
                if task == "research":
                    return {
                        "background": "背景ABC",
                        "species_note": "种属说明",
                        "design_implications": ["口袋保守"],
                        "caveats": [],
                    }
                if task == "findings":
                    return {"overview": "分析DEF", "failure_reasons": ["脂溶性"], "next_steps": ["缩短碳链"]}
                raise AssertionError(task)

        fake = Fake()
        result = self._run("[基本信息]\n靶点: AKT1\n", Settings(api_key="sk-test"), llm=fake)
        report = Path(result["report"]).read_text(encoding="utf-8")
        self.assertEqual(fake.tasks, ["analyze", "research", "findings"])
        self.assertIn("模型摘要XYZ", report)
        self.assertIn("背景ABC", report)
        self.assertIn("分析DEF", report)
        self.assertIn("LLM 增强", result["mode"])
        self.assertNotIn("sk-test", _tree_text(Path(result["run_dir"])))

        class Boom:
            enabled = True

            def complete_json(self, task, payload):
                raise RuntimeError("down")

        failed = self._run("[基本信息]\n靶点: AKT1\n", Settings(api_key="sk-test"), llm=Boom())
        self.assertIn("回退", failed["mode"])
        fallback = Path(failed["report"]).read_text(encoding="utf-8")
        self.assertIn("要求分析", fallback)

    def test_shipped_example_runs(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = run_pipeline(
                request_path=ROOT / "input" / "request.txt",
                settings=Settings(calc_only=True),
                output_root=Path(tmp),
                project_root=ROOT,
                http_get=_no_network,
                run_id="example",
            )
            report = Path(result["report"]).read_text(encoding="utf-8")
            self.assertIn("仅计算", report)
            data = json.loads(Path(result["run_dir"], "calculations.json").read_text(encoding="utf-8"))
            self.assertEqual(data["summary"]["n_input"], 3)
            if not data["summary"]["rdkit"]:
                self.assertIn("没有 RDKit", report)
                return
            rows = {row["input_smiles"]: row for row in data["molecules"]}
            self.assertEqual(rows["CC(=O)Oc1ccccc1C(=O)O"]["numeric_status"], "incomplete")
            self.assertEqual(rows["CC(C)Cc1ccc(cc1)C(C)C(=O)O"]["numeric_status"], "incomplete")
            self.assertEqual(rows["CCCCCCCCCCCCCCCCCCCCCCCC"]["numeric_status"], "fail")
            self.assertTrue(rows["CC(=O)Oc1ccccc1C(=O)O"]["MW"] < 500)


if __name__ == "__main__":
    unittest.main()
