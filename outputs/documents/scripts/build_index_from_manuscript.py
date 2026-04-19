from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
manuscript = ROOT / "Manuscript.qmd"
index = ROOT / "index.qmd"

body = manuscript.read_text(encoding="utf-8").rstrip() + "\n"

frontmatter = """---
title: "Manuscript"
bibliography: vitamins.bib
---

"""

index.write_text(frontmatter + body, encoding="utf-8")
