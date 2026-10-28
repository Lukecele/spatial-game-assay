import json
import tempfile
import unittest
from pathlib import Path

from insitugame.cli import main


class CliTest(unittest.TestCase):
    def test_generates_nf_core_outputs(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            input_path = root / "expression.tsv"
            input_path.write_text(
                "spot\tSensitive\tResistant\tStroma\n"
                "s1\t1\t2\t3\n"
                "s2\t4\t5\t6\n",
                encoding="utf-8",
            )

            self.assertEqual(
                main(["--input", str(input_path), "--output-prefix", str(root / "test")]),
                0,
            )

            result = json.loads(
                (root / "test_results.json").read_text(encoding="utf-8")
            )
            self.assertEqual(result["status"], "completed")
            self.assertEqual(result["n_rows"], 2)
            self.assertTrue((root / "test_tables.csv").is_file())


if __name__ == "__main__":
    unittest.main()
