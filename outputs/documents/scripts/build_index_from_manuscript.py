import re
from pathlib import Path
from typing import Optional


ROOT = Path(__file__).resolve().parents[1]
manuscript = ROOT / "Manuscript.qmd"
index = ROOT / "index.qmd"
quarto_config = ROOT / "_quarto.yml"


def read_book_title(config_path: Path) -> Optional[str]:
    text = config_path.read_text(encoding="utf-8")

    in_book = False
    for line in text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        # Track when we've entered/leaving the top-level `book:` section.
        if re.match(r"^[A-Za-z0-9_-]+:\s*$", line):
            in_book = stripped == "book:"
            continue

        if not in_book:
            continue

        match = re.match(r'^\s{2,}title:\s*"(.*)"\s*$', line)
        if match:
            return match.group(1)

        match = re.match(r"^\s{2,}title:\s*'(.*)'\s*$", line)
        if match:
            return match.group(1)

        match = re.match(r"^\s{2,}title:\s*(.+?)\s*$", line)
        if match:
            return match.group(1)

    return None

body = manuscript.read_text(encoding="utf-8").rstrip() + "\n"
title = read_book_title(quarto_config)

frontmatter_lines = ["---"]
if title:
    frontmatter_lines.append(f'title: "{title}"')
frontmatter_lines.append("bibliography: vitamins.bib")
frontmatter_lines.append("---")
frontmatter = "\n".join(frontmatter_lines) + "\n\n"

index.write_text(frontmatter + body, encoding="utf-8")
