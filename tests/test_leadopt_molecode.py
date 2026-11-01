"""MoleCode 只作为命中后的可选步骤。往返和非法结构拒绝在装了依赖时运行。"""

import csv
import io
import subprocess
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.leadopt.cli import main
from src.leadopt.structure import StructureError, mermaid_to_validated, smiles_to_mermaid, export_hit_list


try:
    import molecode  # noqa: F401
    from rdkit import Chem
except ImportError:
    molecode = None
    Chem = None


_HAVE_MOLECODE = molecode is not None and Chem is not None


class ScreeningBoundaryTests(unittest.TestCase):
    def test_screening_entrypoints_do_not_call_molecode(self):
        for relative in ("run.py", "run_workflow.py", "src/agent/pipeline.py", "setup.sh", "setup.bat"):
            text = (ROOT / relative).read_text(encoding="utf-8")
            self.assertNotIn("molecode", text.lower())
            self.assertNotIn("leadopt", text.lower())


@unittest.skipUnless(
    _HAVE_MOLECODE,
    "可选依赖未安装。执行 python -m pip install -r requirements-leadopt.txt 后再跑本组测试",
)
class MoleCodeLeadOptTests(unittest.TestCase):
    def test_round_trip(self):
        smiles = "CC(=O)Oc1ccccc1C(=O)O"
        graph = smiles_to_mermaid(smiles, name="Aspirin")
        self.assertIn("graph", graph)
        self.assertIn("subgraph", graph)
        self.assertNotIn(smiles, graph)

        validated = mermaid_to_validated(graph)
        canonical = Chem.MolToSmiles(Chem.MolFromSmiles(smiles))
        self.assertEqual(validated.smiles, canonical)
        fenced = mermaid_to_validated("```mermaid\n" + graph + "\n```")
        self.assertEqual(fenced.smiles, canonical)
        self.assertEqual(Chem.MolToSmiles(validated.mol), canonical)
        self.assertEqual(validated.mol.GetNumHeavyAtoms(), 13)

        ethanol = mermaid_to_validated(smiles_to_mermaid("CCO", name="Ethanol"))
        self.assertEqual(ethanol.smiles, "CCO")

    def test_invalid_inputs_are_rejected(self):
        with self.assertRaises(StructureError):
            smiles_to_mermaid("not_a_smiles")
        with self.assertRaises(StructureError):
            smiles_to_mermaid("   ")
        with self.assertRaises(StructureError):
            mermaid_to_validated("this is not a molecule")
        with self.assertRaises(StructureError):
            mermaid_to_validated(
                """
                graph TB
                    subgraph Bad["Bad"]
                        Bad_C_1[Unobtainium]
                        Bad_C_2[CH3]
                        Bad_C_1 --- Bad_C_2
                    end
                """
            )

        pentavalent = """
        graph TB
            subgraph Bad["Bad"]
                Bad_C_1[CH3]
                Bad_C_2[CH3]
                Bad_C_3[CH3]
                Bad_C_1 --- Bad_C_2
                Bad_C_1 --- Bad_C_3
            end
        """
        with self.assertRaises(StructureError):
            mermaid_to_validated(pentavalent)

    def test_constrained_edit(self):
        graph = smiles_to_mermaid("CCO", name="Ethanol")
        self.assertEqual(graph.count("[OH]"), 1)
        edited = graph.replace("[OH]", "[CH3]")

        propane = mermaid_to_validated(edited)
        self.assertEqual(propane.smiles, "CCC")

        kept = mermaid_to_validated(graph, require_smarts="[OX2H]")
        self.assertEqual(kept.smiles, "CCO")
        with self.assertRaises(StructureError):
            mermaid_to_validated(edited, require_smarts="[OX2H]")

    def test_export_hits_leaves_the_source_untouched(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            hits = folder / "hits.csv"
            original = "输入SMILES,标准SMILES,备注\nCCO,CCO,ok\nnot_a_smiles,not_a_smiles,bad\nCCO,CCO,again\n"
            hits.write_text(original, encoding="utf-8")
            out = folder / "leadopt"
            out.mkdir()
            (out / "hit_0002.mmd").write_text("stale\n", encoding="utf-8")
            (out / "notes.txt").write_text("keep\n", encoding="utf-8")
            summary = export_hit_list(hits, out)

            self.assertEqual(hits.read_text(encoding="utf-8"), original)
            self.assertEqual(summary["ok"], 1)
            self.assertEqual(summary["rejected"], 1)
            self.assertEqual(summary["duplicate"], 1)
            self.assertFalse((out / "hit_0002.mmd").exists())
            self.assertEqual((out / "notes.txt").read_text(encoding="utf-8"), "keep\n")
            graph = (out / "hit_0001.mmd").read_text(encoding="utf-8")
            self.assertEqual(mermaid_to_validated(graph).smiles, "CCO")
            self.assertFalse((out / "hit_0001.mmd").read_text(encoding="utf-8").startswith("/"))

            with (out / "index.csv").open(encoding="utf-8", newline="") as handle:
                rows = list(csv.DictReader(handle))
            self.assertEqual([row["status"] for row in rows], ["ok", "rejected", "duplicate"])
            self.assertEqual(rows[0]["file"], "hit_0001.mmd")
            self.assertNotIn(str(folder), (out / "index.csv").read_text(encoding="utf-8"))

            screening = folder / "screening.csv"
            screening.write_text("SMILES,pred_IC50\nc1ccccc1,1\n", encoding="utf-8")
            second = export_hit_list(screening, folder / "from-screening")
            self.assertEqual(second["ok"], 1)
            self.assertEqual(
                mermaid_to_validated((folder / "from-screening" / "hit_0001.mmd").read_text(encoding="utf-8")).smiles,
                "c1ccccc1",
            )

    def test_cli_round_trip_and_rejection(self):
        with tempfile.TemporaryDirectory() as tmp:
            folder = Path(tmp)
            graph_path = folder / "ethanol.mmd"
            out = io.StringIO()
            err = io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                code = main(["smiles-to-mermaid", "CCO", "--name", "Ethanol", "--output", str(graph_path)])
            self.assertEqual(code, 0)
            self.assertTrue(graph_path.is_file())

            out = io.StringIO()
            with redirect_stdout(out), redirect_stderr(err):
                code = main(["mermaid-to-smiles", "--input", str(graph_path), "--require-smarts", "[OX2H]"])
            self.assertEqual(code, 0)
            self.assertEqual(out.getvalue().strip(), "CCO")

            broken = folder / "broken.mmd"
            broken.write_text("not mermaid\n", encoding="utf-8")
            err = io.StringIO()
            with redirect_stdout(io.StringIO()), redirect_stderr(err):
                code = main(["mermaid-to-smiles", "--input", str(broken)])
            self.assertEqual(code, 1)
            self.assertIn("无法解析", err.getvalue())


class CommandHelpTests(unittest.TestCase):
    def test_help_does_not_need_molecode_output_paths(self):
        completed = subprocess.run(
            [sys.executable, str(ROOT / "leadopt.py"), "--help"],
            check=False,
            capture_output=True,
            text=True,
        )
        self.assertEqual(completed.returncode, 0)
        self.assertIn("smiles-to-mermaid", completed.stdout)
        self.assertIn("mermaid-to-smiles", completed.stdout)
        self.assertIn("export-hits", completed.stdout)


if __name__ == "__main__":
    unittest.main()
