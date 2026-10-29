# Local PDF preparation

`pdfprep` prepares local PDFs for reading without an LLM, API key or model calls.
It produces Markdown, exact extracted CSV strings, linked graphic crops, page
previews, the byte-identical original and an extraction report. Results are cached.
The companion `pdf-preprocess` skill calls the same implementation.

## Use from Terminal

```bash
pdfprep ~/Documents/paper.pdf
pdfprep ~/Documents/paper.pdf --output ~/Documents/prepared/paper
pdfprep ~/Documents/papers --recursive --output ~/Documents/prepared
pdfprep --doctor
```

Single-file default: `paper.prepared` beside the PDF. A single-file `--output`
names the package directory itself. A folder's default is a sibling directory
ending in `.prepared`; its relative subdirectories are retained inside the
collection. For multiple explicitly named files, output names get a short path
hash to distinguish matching basenames. Multiple inputs default to `pdfprep-output`
in the current directory. A batch continues after a failed file and exits nonzero
if any extraction failed. Exit zero means packages were produced, not that their
contents were verified; check the report for missing OCR or other warnings.

```text
paper.prepared/
  original.pdf
  document.md
  report.md
  extraction.json
  pages/             # One Markdown, plain-text and PNG file per page
  figures/           # Raster images and vector graphic regions
  tables/            # CSV strings and corresponding PNG crops
```

Read `report.md` first. The JSON manifest stores page labels, one-based PDF page
numbers, native text blocks and bounding boxes, table headers, extraction warnings,
engine/settings identity and hashes of every generated asset. Coordinates are
unrotated PDF page coordinates in points. Per-page Markdown links use `../`;
combined Markdown links are relative to the package root. Keep the whole package
together. In a local Codex session, provide its path; on another machine, transfer
the folder with its assets. Images are opened when needed rather than embedded as
base64 text.

## OCR and limitations

`--ocr auto` (default) attempts local Tesseract OCR when a page has fewer than
40 non-whitespace native characters. If Tesseract is unavailable, the original
page is retained and the failure is recorded. OCR is optional; on macOS it can be
installed separately with `brew install tesseract`. Language data must match
`--language eng` or, for example, `--language eng+ita`.

Use `--ocr never` to skip OCR or `--ocr always` for a corrupt native text layer.
The latter replaces page text with OCR rather than mixing contradictory sources.
OCR pages retain images but do not get reconstructed CSV tables. `--doctor`
reports executable, engine versions and OCR availability. `--dpi 120` is the
default preview resolution; use up to `--dpi 300` for smaller print. OCR itself
uses 300 DPI.

Reading order and heading detection are heuristic. Even two-column layouts may
need visual verification. Equations are **not** converted to verified LaTeX.
Graphic crops may be fragments or include a table. Tiny regions remain in full
page previews. Default table detection uses explicit ruled cells; borderless
tables can be missed. Try `--table-strategy text` for them, but expect false
detections. Empty CSV cells can mean merged cells; values remain strings, and no
financial/numeric interpretation is performed. Captions and footnotes are kept
in page text and previews. Check material values and symbols against the original.
No extraction-quality guarantee or fixed token-saving percentage is implied.

## Cache and preservation

The cache key includes source SHA-256, implementation SHA-256, pinned extraction
engine versions, options and OCR availability/version/languages. Every output
asset is hashed before reuse. Source changes or missing/changed generated assets
trigger a rebuild. `--force` explicitly rebuilds an owned package. Do not edit
generated files as your working notes; a rebuild replaces them.

Unrelated output directories and packages belonging to another source are
refused. Additional personal files inside an existing package are also refused,
so they cannot be discarded during a rebuild. Extraction is staged in a temporary
sibling directory; a failure preserves the previous package. An empty sibling
lock directory prevents concurrent writers to the same package. If a process is
forcibly killed, remove that lock only after confirming no process is using it.
Keep notes outside generated packages.

## Installation and development

Run from the repository with Python 3.10+ and `uv` on PATH:

```bash
python3.12 tools/pdfprep/install.py
```

The installer creates an isolated environment at `~/.local/share/pdfprep/venv`,
the command at `~/.local/bin/pdfprep` and a self-contained skill at
`~/.codex/skills/pdf-preprocess` (or `$CODEX_HOME/skills/pdf-preprocess`). It does
not modify shell configuration. Downloads occur during initial installation;
normal runs require no network. Dependency versions are pinned to avoid changing
extraction results silently. Re-run the installer to update from this checkout.

In Codex, invoke `$pdf-preprocess` with your PDF/folder path, or ask to preprocess
PDFs using `pdfprep`. If the skill is not visible in an existing session, start a
new session. The skill reads reports and selected content; it does not load all
page images automatically.

```bash
~/.local/share/pdfprep/venv/bin/python -m unittest discover -s tools/pdfprep/tests -v
```

Tests cover exact decimal/minus-sign extraction, column order on a fixture,
raster/vector preservation, rotated pages, source immutability, cache reuse and
repair, OCR success/failure reporting, nested batches and output ownership.
OCR integration is tested when Tesseract is installed; the core suite uses a
controlled substitute to exercise reporting without an external dependency.

Dependencies: [PyMuPDF](https://pymupdf.readthedocs.io/) and
[PyMuPDF4LLM](https://pymupdf.readthedocs.io/en/latest/pymupdf4llm/), each under
its upstream licensing terms. This wrapper does not alter those terms.
