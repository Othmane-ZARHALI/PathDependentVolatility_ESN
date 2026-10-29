"""Observable extraction, preservation, caching and output ownership checks."""
import argparse
import csv
import io
import json
from pathlib import Path
import shutil
import sys
import tempfile
import unittest
from unittest.mock import patch

import pymupdf

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
import pdfprep


def make_pdf(path):
    doc = pymupdf.open()
    page = doc.new_page(width=600, height=800)
    page.insert_text((40, 40), "Extraction sample", fontsize=20)
    # Each column should stay together, despite sharing vertical coordinates.
    page.insert_textbox((40, 80, 270, 210), "\n".join(f"LEFT{i} independent observation." for i in range(6)), fontsize=10)
    page.insert_textbox((330, 80, 560, 210), "\n".join(f"RIGHT{i} independent observation." for i in range(6)), fontsize=10)
    page.insert_text((40, 240), "Table 1. Prices in EUR; sample size 2.", fontsize=10)
    xs, ys = [40, 200, 360], [260, 290, 320, 350]
    for x in xs:
        page.draw_line((x, ys[0]), (x, ys[-1]))
    for y in ys:
        page.draw_line((xs[0], y), (xs[-1], y))
    for row, values in enumerate([["Asset", "Price"], ["Alpha", "12.50"], ["Beta", "-3.25"]]):
        for col, value in enumerate(values):
            page.insert_text((xs[col] + 8, ys[row] + 20), value, fontsize=11)
    page.insert_text((40, 370), "Note: negative prices are retained literally.", fontsize=10)
    page.draw_rect(pymupdf.Rect(40, 420, 220, 500), color=(0, 0, 1))
    page.draw_line((40, 495), (210, 425), color=(1, 0, 0), width=3)
    page.insert_text((40, 520), "Figure 1. A vector diagram.", fontsize=10)
    # A raster-only scan on a second page.
    pix = page.get_pixmap(dpi=100)
    scan = doc.new_page(width=600, height=800)
    scan.insert_image(scan.rect, stream=pix.tobytes("png"))
    doc.save(path)
    doc.close()


class PreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.source = self.root / "sample.pdf"
        self.dest = self.root / "prepared"
        make_pdf(self.source)
        self.original_hash = pdfprep.digest(self.source)

    def tearDown(self):
        self.temp.cleanup()

    def run_cli(self, *extra, source=None, dest=None):
        stdout, stderr = io.StringIO(), io.StringIO()
        with patch("sys.stdout", stdout), patch("sys.stderr", stderr):
            status = pdfprep.main([str(source or self.source), "-o", str(dest or self.dest), "--ocr", "never", *extra])
        return status, stdout.getvalue(), stderr.getvalue()

    def test_extract_preserves_values_images_and_column_order(self):
        status, out, err = self.run_cli()
        self.assertEqual(status, 0, err)
        data = json.loads((self.dest / "extraction.json").read_text())
        self.assertEqual(pdfprep.digest(self.source), self.original_hash)
        self.assertEqual(pdfprep.digest(self.dest / "original.pdf"), self.original_hash)
        self.assertEqual(data["summary"]["pages"], 2)
        tables = data["pages"][0]["tables"]
        self.assertEqual(len(tables), 1)
        with (self.dest / tables[0]["csv"]).open(newline="") as stream:
            self.assertEqual(list(csv.reader(stream)), [["Asset", "Price"], ["Alpha", "12.50"], ["Beta", "-3.25"]])
        md = (self.dest / "document.md").read_text()
        self.assertIn("sample size 2", md)
        self.assertIn("negative prices", md)
        self.assertLess(md.index("LEFT5"), md.index("RIGHT0"))
        self.assertTrue(any(f["kind"] == "vector" for f in data["pages"][0]["figures"]))
        self.assertTrue(any(f["kind"] == "raster" for f in data["pages"][1]["figures"]))
        self.assertIn("../original.pdf#page=1", (self.dest / "pages/page-0001.md").read_text())
        for name, expected in data["files"].items():
            self.assertEqual(pdfprep.digest(self.dest / name), expected)
        self.assertTrue(data["pages"][1]["needs_review"])

    def test_cache_and_corrupt_output_recovery(self):
        self.assertEqual(self.run_cli()[0], 0)
        manifest = self.dest / "extraction.json"
        before = manifest.stat().st_mtime_ns
        status, out, err = self.run_cli()
        self.assertEqual(status, 0, err)
        self.assertIn("CACHED", out)
        self.assertEqual(manifest.stat().st_mtime_ns, before)
        (self.dest / "document.md").write_text("broken")
        status, out, err = self.run_cli()
        self.assertEqual(status, 0, err)
        self.assertIn("READY", out)
        self.assertIn("Extraction sample", (self.dest / "document.md").read_text())

    def test_source_change_rebuilds(self):
        self.assertEqual(self.run_cli()[0], 0)
        doc = pymupdf.open(self.source)
        doc[0].insert_text((40, 700), "ADDED MATERIAL")
        doc.saveIncr()
        doc.close()
        status, out, err = self.run_cli()
        self.assertEqual(status, 0, err)
        self.assertIn("READY", out)
        self.assertIn("ADDED MATERIAL", (self.dest / "document.md").read_text())

    def test_missing_or_failed_ocr_is_flagged(self):
        with patch.object(pdfprep, "run_ocr", side_effect=RuntimeError("Tesseract unavailable")):
            status, _, err = self.run_cli("--ocr", "auto")
        self.assertEqual(status, 0, err)
        data = json.loads((self.dest / "extraction.json").read_text())
        self.assertIn("Tesseract unavailable", data["pages"][1]["warnings"])
        self.assertEqual(data["pages"][1]["extracted_characters"], 0)

    def test_ocr_success_recorded_without_invented_tables(self):
        with patch.object(pdfprep, "run_ocr", return_value="Recognised scan content"):
            status, _, err = self.run_cli("--ocr", "always")
        self.assertEqual(status, 0, err)
        data = json.loads((self.dest / "extraction.json").read_text())
        self.assertEqual(data["summary"]["ocr_pages"], [1, 2])
        self.assertEqual(data["summary"]["tables"], 0)
        self.assertIn("Recognised scan content", (self.dest / "document.md").read_text())

    @unittest.skipUnless(shutil.which("tesseract"), "Optional Tesseract is not installed")
    def test_real_tesseract_reads_scan(self):
        status, _, err = self.run_cli("--ocr", "auto")
        self.assertEqual(status, 0, err)
        data = json.loads((self.dest / "extraction.json").read_text())
        self.assertEqual(data["pages"][1]["method"], "ocr")
        self.assertIn("Extraction sample", (self.dest / "pages/page-0002.txt").read_text())

    def test_unrelated_output_refused_even_with_force(self):
        self.dest.mkdir()
        note = self.dest / "personal.txt"
        note.write_text("KEEP")
        status, _, _ = self.run_cli("--force")
        self.assertEqual(status, 1)
        self.assertEqual(note.read_text(), "KEEP")

    def test_added_personal_files_refused(self):
        self.assertEqual(self.run_cli()[0], 0)
        note = self.dest / "my-notes.txt"
        note.write_text("KEEP")
        self.assertEqual(self.run_cli("--force")[0], 1)
        self.assertEqual(note.read_text(), "KEEP")

    def test_failure_leaves_existing_package(self):
        self.assertEqual(self.run_cli()[0], 0)
        before = (self.dest / "extraction.json").read_bytes()
        with patch.object(pdfprep, "extract", side_effect=RuntimeError("interrupted")):
            self.assertEqual(self.run_cli("--force")[0], 1)
        self.assertEqual((self.dest / "extraction.json").read_bytes(), before)
        self.assertFalse(list(self.root.glob(".*.building-*")))
        self.assertFalse(list(self.root.glob("*.pdfprep-lock")))

    def test_nested_batch_does_not_ingest_own_outputs(self):
        folder = self.root / "papers"
        (folder / "nested").mkdir(parents=True)
        (folder / "nested/sample.PDF").write_bytes(self.source.read_bytes())
        output = folder / "prepared"
        status, _, err = self.run_cli("--recursive", source=folder, dest=output)
        self.assertEqual(status, 0, err)
        args = argparse.Namespace(inputs=[str(folder)], output=str(output), recursive=True)
        jobs = pdfprep.jobs_for(args)
        self.assertEqual(len(jobs), 1)
        self.assertEqual(jobs[0][1], (output / "nested/sample").resolve())

    def test_rotated_page_and_linked_images(self):
        with pymupdf.open(self.source) as doc:
            doc[0].set_rotation(90)
            doc.saveIncr()
        status, _, err = self.run_cli()
        self.assertEqual(status, 0, err)
        data = json.loads((self.dest / "extraction.json").read_text())
        page = data["pages"][0]
        self.assertEqual(page["rotation"], 90)
        self.assertFalse(any("failed" in w for w in page["warnings"]), page["warnings"])
        self.assertEqual(len(page["tables"]), 1)


if __name__ == "__main__":
    unittest.main()
