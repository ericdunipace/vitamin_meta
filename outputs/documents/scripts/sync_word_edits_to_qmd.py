#!/usr/bin/env python3
"""Merge collaborator Word edits back into Quarto markdown.

The script aligns source `.qmd` paragraphs to:
1. a baseline rendered `.docx` produced from that source, and
2. an edited `.docx` from collaborators.

For matched paragraphs, it preserves inline Quarto expressions such as
`r ...` and cross-references like `@fig-...` while transplanting the edited
Word prose around them.

This is heuristic. It works best when:
- the edited paragraph still resembles the rendered source paragraph,
- rendered values from inline expressions are still present in the edit, and
- paragraphs have not been split/merged dramatically in Word.
"""

from __future__ import annotations

import argparse
import difflib
import re
import sys
import xml.etree.ElementTree as ET
import zipfile
from dataclasses import dataclass
from pathlib import Path


NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}


INLINE_EXPR_RE = re.compile(r"`r\s+[^`]+`")
CROSSREF_RE = re.compile(r"@(?:fig|tbl|sec|eq)-[A-Za-z0-9_.:-]+")


@dataclass
class QmdParagraph:
    start_line: int
    end_line: int
    text: str


@dataclass
class Placeholder:
    token: str
    rendered: str


def normalize(text: str) -> str:
    text = text.replace("\u00A0", " ")
    text = re.sub(r"\s+", " ", text.strip())
    return text


def normalize_for_match(text: str) -> str:
    text = INLINE_EXPR_RE.sub(" ", text)
    text = CROSSREF_RE.sub(" ", text)
    text = re.sub(r"`[^`]+`", " ", text)
    text = text.replace("—", "-")
    text = normalize(text).lower()
    text = re.sub(r"[^a-z0-9%.,;:() -]", " ", text)
    return normalize(text)


def parse_qmd_paragraphs(path: Path) -> list[QmdParagraph]:
    lines = path.read_text(encoding="utf-8").splitlines()
    paragraphs: list[QmdParagraph] = []

    in_fence = False
    fence_delim = None
    current: list[str] = []
    start_line = 1

    def flush(end_line: int) -> None:
        nonlocal current, start_line
        if not current:
            return
        text = "\n".join(current).strip()
        current = []
        if not text:
            return
        if text.startswith("#"):
            return
        if text.startswith(":::") or text.endswith(":::"):
            return
        paragraphs.append(QmdParagraph(start_line=start_line, end_line=end_line, text=text))

    for idx, line in enumerate(lines, start=1):
        stripped = line.strip()

        if stripped.startswith("```") or stripped.startswith("~~~"):
            delim = stripped[:3]
            if not in_fence:
                flush(idx - 1)
                in_fence = True
                fence_delim = delim
            elif delim == fence_delim:
                in_fence = False
                fence_delim = None
            continue

        if in_fence:
            continue

        if stripped == "":
            flush(idx - 1)
            start_line = idx + 1
            continue

        if not current:
            start_line = idx
        current.append(line)

    flush(len(lines))
    return paragraphs


def docx_paragraph_texts(path: Path) -> list[str]:
    with zipfile.ZipFile(path) as zf:
        root = ET.fromstring(zf.read("word/document.xml"))

    out: list[str] = []
    for p in root.findall(".//w:p", NS):
        text = "".join(t.text for t in p.findall(".//w:t", NS) if t.text)
        text = normalize(text)
        if text:
            out.append(text)
    return out


def similarity(a: str, b: str) -> float:
    return difflib.SequenceMatcher(None, normalize_for_match(a), normalize_for_match(b)).ratio()


def find_best_match(
    text: str,
    candidates: list[str],
    threshold: float = 0.55,
    expected_idx: int | None = None,
    window: int = 8,
) -> tuple[int, float] | None:
    if expected_idx is not None:
        lo = max(0, expected_idx - window)
        hi = min(len(candidates), expected_idx + window + 1)
        candidate_indexes = list(range(lo, hi))
    else:
        candidate_indexes = list(range(len(candidates)))

    best_idx = -1
    best_score = -1.0
    for idx in candidate_indexes:
        score = similarity(text, candidates[idx])
        if score > best_score:
            best_idx = idx
            best_score = score

    if expected_idx is not None:
        global_best = find_best_match(text, candidates, threshold=threshold, expected_idx=None, window=window)
        if global_best is not None and global_best[1] > best_score + 0.08:
            best_idx, best_score = global_best

    if best_idx < 0 or best_score < threshold:
        return None
    return best_idx, best_score


def extract_placeholders(source_text: str, rendered_text: str) -> list[Placeholder] | None:
    token_re = re.compile(r"(`r\s+[^`]+`|@(?:fig|tbl|sec|eq)-[A-Za-z0-9_.:-]+)")
    parts = token_re.split(source_text)
    tokens = token_re.findall(source_text)

    if not tokens:
        return []

    pos = 0
    placeholders: list[Placeholder] = []
    for i, token in enumerate(tokens):
        literal_before = normalize(parts[i * 2])
        literal_after = normalize(parts[i * 2 + 2]) if (i * 2 + 2) < len(parts) else ""

        if literal_before:
            found = rendered_text.find(literal_before, pos)
            if found == -1:
                return None
            pos = found + len(literal_before)

        if literal_after:
            next_pos = rendered_text.find(literal_after, pos)
            if next_pos == -1:
                return None
            rendered = rendered_text[pos:next_pos]
            pos = next_pos
        else:
            rendered = rendered_text[pos:]
            pos = len(rendered_text)

        rendered = normalize(rendered)
        if not rendered:
            return None
        placeholders.append(Placeholder(token=token, rendered=rendered))

    return placeholders


def transplant_placeholders(edited_text: str, placeholders: list[Placeholder]) -> tuple[str, int]:
    result = edited_text
    replaced = 0
    search_from = 0

    for ph in placeholders:
        if not ph.rendered:
            continue

        idx = result.find(ph.rendered, search_from)
        if idx == -1:
            idx = result.find(ph.rendered)
        if idx == -1:
            continue

        result = result[:idx] + ph.token + result[idx + len(ph.rendered):]
        search_from = idx + len(ph.token)
        replaced += 1

    return result, replaced


def merge_qmd(
    qmd_path: Path,
    baseline_docx: Path,
    edited_docx: Path,
    write: bool,
    min_score: float,
    min_replaced: int,
) -> int:
    paragraphs = parse_qmd_paragraphs(qmd_path)
    baseline_paragraphs = docx_paragraph_texts(baseline_docx)
    edited_paragraphs = docx_paragraph_texts(edited_docx)

    source_lines = qmd_path.read_text(encoding="utf-8").splitlines()
    replacements: list[tuple[int, int, str]] = []

    for para in paragraphs:
        if not (INLINE_EXPR_RE.search(para.text) or CROSSREF_RE.search(para.text)):
            continue

        base_match = find_best_match(para.text, baseline_paragraphs, threshold=min_score)
        if not base_match:
            continue

        base_idx, base_score = base_match
        baseline_text = baseline_paragraphs[base_idx]
        placeholders = extract_placeholders(para.text, baseline_text)
        if placeholders is None:
            continue

        edited_match = find_best_match(
            baseline_text,
            edited_paragraphs,
            threshold=min_score,
            expected_idx=base_idx,
        )
        if not edited_match:
            continue

        edited_idx, edited_score = edited_match
        edited_text = edited_paragraphs[edited_idx]
        merged_text, replaced = transplant_placeholders(edited_text, placeholders)

        if replaced < max(min_replaced, len([p for p in placeholders if p.token.startswith("`r ")]) // 2):
            continue

        if normalize(merged_text) == normalize(para.text):
            continue

        replacements.append((para.start_line, para.end_line, merged_text))
        print(
            f"matched lines {para.start_line}-{para.end_line} | "
            f"baseline={base_score:.2f} edited={edited_score:.2f} replaced={replaced}",
            file=sys.stderr,
        )

    if not replacements:
        print("No paragraph updates found.", file=sys.stderr)
        return 0

    new_lines = source_lines[:]
    for start_line, end_line, new_text in reversed(replacements):
        new_block = new_text.splitlines()
        new_lines[start_line - 1:end_line] = new_block

    if write:
        qmd_path.write_text("\n".join(new_lines) + "\n", encoding="utf-8")
    else:
        print("\n".join(new_lines))

    print(f"Updated {len(replacements)} paragraph(s).", file=sys.stderr)
    return len(replacements)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qmd", required=True, type=Path)
    parser.add_argument("--baseline-docx", required=True, type=Path)
    parser.add_argument("--edited-docx", required=True, type=Path)
    parser.add_argument("--write", action="store_true", help="Write changes back to the qmd file.")
    parser.add_argument("--min-score", type=float, default=0.55, help="Minimum fuzzy match score.")
    parser.add_argument(
        "--min-replaced",
        type=int,
        default=1,
        help="Minimum number of placeholders that must be reinserted.",
    )
    args = parser.parse_args()

    merge_qmd(
        qmd_path=args.qmd,
        baseline_docx=args.baseline_docx,
        edited_docx=args.edited_docx,
        write=args.write,
        min_score=args.min_score,
        min_replaced=args.min_replaced,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
