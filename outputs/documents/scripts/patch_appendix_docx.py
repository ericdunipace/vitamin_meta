#!/usr/bin/env python3
"""Patch DOCX to apply appendix heading styles based on hidden markers.

Usage:
  patch_appendix_docx.py input.docx
  patch_appendix_docx.py input.docx output.docx

If no args are provided, the script will try to read Quarto's
post-render output list from environment variables and patch all
`.docx` outputs.
"""
from __future__ import annotations

import os
import sys
import zipfile
import tempfile
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET
import re
from typing import Iterable

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
ET.register_namespace("w", NS["w"])

MARKER_TO_STYLE = {
    "APPENDIX_H2": "AppendixHeading1",
    "APPENDIX_H3": "AppendixHeading2",
    "APPENDIX_H4": "AppendixHeading3",
    "APPENDIX_H5": "AppendixHeading4",
}

STYLE_NAME_HINTS = {
    "AppendixHeading1": "Appendix Heading 1",
    "AppendixHeading2": "Appendix Heading 2",
    "AppendixHeading3": "Appendix Heading 3",
    "AppendixHeading4": "Appendix Heading 4",
}

DOCX_HEADING1_STYLES = {"Heading1", "heading 1", "Heading 1"}
APPENDIX_HEADING1_STYLES = {"AppendixHeading1"}
DOCX_HEADING_STYLE_IDS = {f"Heading{i}" for i in range(1, 10)}
CROSSREF_PREFIXES = (
    "sec-",
    "fig-",
    "tbl-",
    "eq-",
    "lst-",
    "thm-",
    "lem-",
    "cor-",
    "prp-",
    "cnj-",
    "def-",
    "exm-",
    "exr-",
    "ch-",
    "apx-",
)
CITATION_PATTERN = re.compile(r"@([A-Za-z0-9_][A-Za-z0-9_:.#$%&+?/\-]*)")


def eprint(*args: object) -> None:
    print(*args, file=sys.stderr)


def _read_quarto_output_files() -> list[Path]:
    files: list[str] = []

    env_list = os.environ.get("QUARTO_PROJECT_OUTPUT_FILES")
    if env_list:
        # Some Quarto versions use newline-separated lists; others may use os.pathsep.
        if "\n" in env_list:
            files.extend([line for line in env_list.splitlines() if line.strip()])
        else:
            files.extend([p for p in env_list.split(os.pathsep) if p.strip()])

    env_file = os.environ.get("QUARTO_USE_FILE_FOR_PROJECT_OUTPUT_FILES")
    if not env_file:
        env_file = os.environ.get("QUARTO_PROJECT_OUTPUT_FILES_FILE")
    if env_file and Path(env_file).exists():
        try:
            contents = Path(env_file).read_text(encoding="utf-8")
            files.extend([line for line in contents.splitlines() if line.strip()])
        except Exception:
            pass

    # Fallback: common output dir
    if not files:
        out_dir = os.environ.get("QUARTO_PROJECT_OUTPUT_DIR")
        if out_dir and Path(out_dir).exists():
            files.extend(str(p) for p in Path(out_dir).rglob("*.docx"))
    # Fallback: common Quarto output locations relative to CWD
    if not files:
        cwd = Path.cwd()
        for candidate in (cwd / "_book", cwd / "_output", cwd):
            if candidate.exists():
                files.extend(str(p) for p in candidate.rglob("*.docx"))

    # De-dup while preserving order
    seen = set()
    out: list[Path] = []
    for f in files:
        p = Path(f)
        if p in seen:
            continue
        seen.add(p)
        out.append(p)
    return out


def _extract_citation_ids(text: str) -> set[str]:
    ids: set[str] = set()
    for match in CITATION_PATTERN.finditer(text):
        key = match.group(1).rstrip(".,;:")
        if any(key.startswith(prefix) for prefix in CROSSREF_PREFIXES):
            continue
        ids.add(key)
    return ids


def _citation_scopes(project_dir: Path) -> dict[str, set[str]]:
    scopes: dict[str, set[str]] = {}
    chapter_files = {
        "main": project_dir / "index.qmd",
        "appendix": project_dir / "Appendix.qmd",
    }
    for scope, path in chapter_files.items():
        if path.exists():
            scopes[scope] = _extract_citation_ids(path.read_text(encoding="utf-8"))
    return scopes


def _style_id_for_paragraph(p: ET.Element) -> str | None:
    ppr = p.find("w:pPr", NS)
    pstyle = ppr.find("w:pStyle", NS) if ppr is not None else None
    return pstyle.get(f"{{{NS['w']}}}val") if pstyle is not None else None


def _top_level_text(elem: ET.Element) -> str:
    texts = []
    for t in elem.findall(".//w:t", NS):
        if t.text:
            texts.append(t.text)
    return "".join(texts)


def _normalize_heading_text(text: str) -> str:
    stripped = text.strip()
    return re.sub(r"^\d+(?:\.\d+)*\.\s*", "", stripped).strip().lower()


def _filter_docx_bibliographies(root: ET.Element, citation_scopes: dict[str, set[str]], debug: bool = False) -> int:
    body = root.find("w:body", NS)
    if body is None or not citation_scopes:
        return 0

    children = list(body)
    sections: list[tuple[str, int, int]] = []
    in_appendix = False

    i = 0
    while i < len(children):
        node = children[i]
        if node.tag != f"{{{NS['w']}}}p":
            i += 1
            continue

        text = _normalize_heading_text(_top_level_text(node))
        style_id = _style_id_for_paragraph(node)

        if style_id in DOCX_HEADING1_STYLES and text == "appendix":
            in_appendix = True
            i += 1
            continue

        is_references = text == "references" and (
            style_id in DOCX_HEADING1_STYLES or style_id in APPENDIX_HEADING1_STYLES
        )
        if not is_references:
            i += 1
            continue

        scope = "appendix" if in_appendix or style_id in APPENDIX_HEADING1_STYLES else "main"
        start = i + 1
        end = start
        seen_bibliography = False
        while end < len(children):
            child = children[end]
            if child.tag == f"{{{NS['w']}}}p":
                child_style = _style_id_for_paragraph(child)
                if child_style == "Bibliography":
                    seen_bibliography = True
                    end += 1
                    continue
                if seen_bibliography:
                    break
                end += 1
                continue
            if child.tag in (f"{{{NS['w']}}}bookmarkStart", f"{{{NS['w']}}}bookmarkEnd"):
                end += 1
                continue
            if seen_bibliography:
                break
            end += 1

        sections.append((scope, start, end))
        i = end

    removed = 0
    for scope, start, end in reversed(sections):
        allowed = citation_scopes.get(scope)
        if not allowed:
            continue

        section_nodes = children[start:end]
        filtered: list[ET.Element] = []
        k = 0
        while k < len(section_nodes):
            node = section_nodes[k]
            if node.tag == f"{{{NS['w']}}}bookmarkStart":
                name = node.get(f"{{{NS['w']}}}name")
                if name == "refs":
                    filtered.append(node)
                    k += 1
                    continue
                if name and name.startswith("ref-"):
                    entry_id = name[4:]
                    entry_nodes = [node]
                    k += 1
                    while k < len(section_nodes):
                        next_node = section_nodes[k]
                        next_name = (
                            next_node.get(f"{{{NS['w']}}}name")
                            if next_node.tag == f"{{{NS['w']}}}bookmarkStart"
                            else None
                        )
                        if next_name and next_name.startswith("ref-"):
                            break
                        entry_nodes.append(next_node)
                        k += 1
                    if entry_id in allowed:
                        filtered.extend(entry_nodes)
                    else:
                        removed += 1
                    continue
            filtered.append(node)
            k += 1

        for node in section_nodes:
            body.remove(node)
        for offset, node in enumerate(filtered):
            body.insert(start + offset, node)
        if debug:
            eprint(f"Filtered {scope} bibliography: kept {len(filtered)} XML nodes")

    return removed


def _strip_main_heading_numbering(styles_xml: bytes, numbering_xml: bytes | None, debug: bool = False) -> tuple[bytes, bytes | None, int]:
    styles_root = ET.fromstring(styles_xml)
    removed = 0

    for style in styles_root.findall("w:style", NS):
        style_id = style.get(f"{{{NS['w']}}}styleId")
        if style_id not in DOCX_HEADING_STYLE_IDS:
            continue
        ppr = style.find("w:pPr", NS)
        if ppr is None:
            continue
        numpr = ppr.find("w:numPr", NS)
        if numpr is not None:
            ppr.remove(numpr)
            removed += 1
            if debug:
                eprint(f"Removed numbering from style {style_id}")

    new_numbering_xml = numbering_xml
    if numbering_xml:
        numbering_root = ET.fromstring(numbering_xml)
        for lvl in numbering_root.findall(".//w:lvl", NS):
            pstyle = lvl.find("w:pStyle", NS)
            if pstyle is None:
                continue
            val = pstyle.get(f"{{{NS['w']}}}val")
            if val in DOCX_HEADING_STYLE_IDS:
                lvl.remove(pstyle)
                removed += 1
                if debug:
                    eprint(f"Removed numbering mapping for style {val}")
        new_numbering_xml = ET.tostring(numbering_root, encoding="utf-8", xml_declaration=True)

    return (
        ET.tostring(styles_root, encoding="utf-8", xml_declaration=True),
        new_numbering_xml,
        removed,
    )


def _patch_document_xml(
    xml_bytes: bytes,
    citation_scopes: dict[str, set[str]] | None = None,
    debug: bool = False,
) -> tuple[bytes, int]:
    root = ET.fromstring(xml_bytes)
    patched = 0
    caption_patched = 0

    def strip_leading_number(p: ET.Element) -> None:
        runs = []
        for r in p.findall("w:r", NS):
            t = r.find("w:t", NS)
            if t is not None and t.text:
                runs.append((r, t))
        if not runs:
            return
        full = "".join(t.text for _, t in runs)
        stripped = full.lstrip()
        if not stripped:
            return
        leading_ws = len(full) - len(stripped)
        parts = stripped.split(None, 1)
        if not parts:
            return
        token = parts[0]
        match = re.fullmatch(r"(?:[A-Z]\.\d+(?:\.\d+)*|\d+\.\d+(?:\.\d+)*)", token)
        if not match:
            return
        # remove token + following whitespace
        after = stripped[len(token):]
        ws = 0
        for ch in after:
            if ch.isspace() or ch == "\u00A0":
                ws += 1
            else:
                break
        remaining = leading_ws + len(token) + ws
        if debug:
            eprint(f"Stripping section prefix '{full[:remaining]}'")
        for r, t in runs:
            if remaining <= 0:
                break
            text = t.text or ""
            if remaining >= len(text):
                remaining -= len(text)
                t.text = ""
            else:
                t.text = text[remaining:]
                remaining = 0
        # Remove empty runs left behind by stripping.
        for r, t in runs:
            if (t.text or "") == "":
                p.remove(r)

    def get_paragraph_text(p: ET.Element) -> str:
        texts = []
        for t in p.findall(".//w:t", NS):
            if t.text:
                texts.append(t.text)
        return "".join(texts)

    def set_paragraph_text(p: ET.Element, new_text: str) -> None:
        runs = []
        for r in p.findall("w:r", NS):
            t = r.find("w:t", NS)
            if t is not None:
                runs.append((r, t))
        if not runs:
            # create a run if none exist
            r = ET.SubElement(p, f"{{{NS['w']}}}r")
            t = ET.SubElement(r, f"{{{NS['w']}}}t")
            t.text = new_text
            return

        remaining = new_text
        for r, t in runs:
            if remaining == "":
                t.text = ""
                continue
            original_len = len(t.text or "")
            if original_len == 0:
                t.text = ""
                continue
            t.text = remaining[:original_len]
            remaining = remaining[original_len:]
        if remaining:
            # append remaining text in a new run
            r = ET.SubElement(p, f"{{{NS['w']}}}r")
            t = ET.SubElement(r, f"{{{NS['w']}}}t")
            t.text = remaining

    appendix_style_ids = set(MARKER_TO_STYLE.values())

    in_appendix = False
    appendix_letter = None
    appendix_index = 0
    lvl2 = lvl3 = lvl4 = 0
    caption_counters = {"Table": 0, "Figure": 0}
    caption_map = {"Table": {}, "Figure": {}}
    section_label_map = {}

    caption_re = re.compile(r"^(Table|Figure)[\u00A0 ]+(\d+(?:\.\d+)*)")

    for p in root.findall(".//w:p", NS):
        found_marker = None
        # Scan runs for markers.
        for r in list(p.findall("w:r", NS)):
            t = r.find("w:t", NS)
            if t is None or t.text is None:
                continue
            marker = t.text.strip()
            if marker in MARKER_TO_STYLE:
                found_marker = marker
                # Remove the entire run containing the marker.
                p.remove(r)

        if not found_marker:
            # If this paragraph already has an appendix style (e.g., rerun),
            # still strip any leading section numbering.
            ppr = p.find("w:pPr", NS)
            if ppr is not None:
                pstyle = ppr.find("w:pStyle", NS)
                if pstyle is not None:
                    current = pstyle.get(f"{{{NS['w']}}}val")
                    if current in appendix_style_ids:
                        strip_leading_number(p)
            continue

        style_id = MARKER_TO_STYLE[found_marker]
        if debug:
            eprint(f"Applying {style_id} to paragraph")

        ppr = p.find("w:pPr", NS)
        if ppr is None:
            ppr = ET.Element(f"{{{NS['w']}}}pPr")
            p.insert(0, ppr)

        pstyle = ppr.find("w:pStyle", NS)
        if pstyle is None:
            pstyle = ET.SubElement(ppr, f"{{{NS['w']}}}pStyle")

        pstyle.set(f"{{{NS['w']}}}val", style_id)
        # Remove any leading section number injected by Quarto numbering.
        strip_leading_number(p)
        patched += 1

    # Second pass: update appendix caption numbering + in-text references
    # and build a map of appendix section bookmarks -> appendix-style labels.
    body = root.find("w:body", NS)
    pending_bookmarks = []

    def capture_bookmarks(elem: ET.Element) -> None:
        if elem.tag == f"{{{NS['w']}}}bookmarkStart":
            name = elem.get(f"{{{NS['w']}}}name")
            if name:
                pending_bookmarks.append(name)
            return
        for b in elem.findall("w:bookmarkStart", NS):
            name = b.get(f"{{{NS['w']}}}name")
            if name:
                pending_bookmarks.append(name)

    if body is None:
        return ET.tostring(root, encoding="utf-8", xml_declaration=True), patched

    for node in list(body):
        capture_bookmarks(node)
        if node.tag != f"{{{NS['w']}}}p":
            continue

        p = node
        ppr = p.find("w:pPr", NS)
        pstyle = ppr.find("w:pStyle", NS) if ppr is not None else None
        style_id = pstyle.get(f"{{{NS['w']}}}val") if pstyle is not None else None

        text = get_paragraph_text(p)
        normalized_text = _normalize_heading_text(text)
        if not text:
            continue

        # Detect appendix boundary (level-1 heading titled 'Appendix')
        if style_id in ("Heading1", "heading 1", "Heading 1") and normalized_text == "appendix":
            in_appendix = True
            appendix_letter = None
            appendix_index = 0
            lvl2 = lvl3 = lvl4 = 0
            caption_counters = {"Table": 0, "Figure": 0}
            caption_map = {"Table": {}, "Figure": {}}
            continue

        if not in_appendix:
            continue

        # Track current appendix letter based on Appendix Heading 1
        if style_id == "AppendixHeading1":
            appendix_index += 1
            appendix_letter = chr(ord("A") + appendix_index - 1)
            lvl2 = lvl3 = lvl4 = 0
            caption_counters = {"Table": 0, "Figure": 0}
            caption_map = {"Table": {}, "Figure": {}}

        if not appendix_letter:
            continue

        # Build section label for appendix headings
        section_label = None
        if style_id == "AppendixHeading1":
            section_label = f"{appendix_letter}"
        elif style_id == "AppendixHeading2":
            lvl2 += 1
            lvl3 = lvl4 = 0
            section_label = f"{appendix_letter}.{lvl2}"
        elif style_id == "AppendixHeading3":
            lvl3 += 1
            lvl4 = 0
            section_label = f"{appendix_letter}.{lvl2}.{lvl3}"
        elif style_id == "AppendixHeading4":
            lvl4 += 1
            section_label = f"{appendix_letter}.{lvl2}.{lvl3}.{lvl4}"

        if section_label:
            # capture any bookmarks on this paragraph, plus any pending ones
            local_bookmarks = [b.get(f"{{{NS['w']}}}name") for b in p.findall("w:bookmarkStart", NS)]
            for name in pending_bookmarks + [b for b in local_bookmarks if b]:
                section_label_map[name] = section_label
            pending_bookmarks = []

        new_text = text

        # Caption line?
        m = caption_re.match(text)
        if m:
            label = m.group(1)
            old_num = m.group(2)
            caption_counters[label] += 1
            new_num = f"{appendix_letter}.{caption_counters[label]}"
            caption_map[label][old_num] = new_num
            if debug:
                eprint(f"Caption match: '{text[:80]}' -> {label} {new_num}")
            new_text = re.sub(
                r"^(%s)[\u00A0 ]+%s" % (label, re.escape(old_num)),
                f"{label} {new_num}",
                new_text,
                count=1,
            )

        # Update in-text references in appendix using known mappings
        for label, mapping in caption_map.items():
            for old_num, new_num in mapping.items():
                new_text = re.sub(
                    r"\\b%s[\u00A0 ]+%s\\b" % (label, re.escape(old_num)),
                    f"{label} {new_num}",
                    new_text,
                )

        if new_text != text:
            set_paragraph_text(p, new_text)
            caption_patched += 1

    # Third pass: update section crossref hyperlink text in appendix
    for h in root.findall(".//w:hyperlink", NS):
        anchor = h.get(f"{{{NS['w']}}}anchor")
        if not anchor or anchor not in section_label_map:
            continue
        label = section_label_map[anchor]
        current = "".join(t.text for t in h.findall(".//w:t", NS) if t.text)
        if not current:
            continue
        if "section" in current.lower():
            new_text = f"Section {label}"
        else:
            new_text = label
        # replace hyperlink text runs
        runs = h.findall(".//w:r", NS)
        if not runs:
            continue
        first = True
        for r in runs:
            t = r.find("w:t", NS)
            if t is None:
                continue
            if first:
                t.text = new_text
                first = False
            else:
                t.text = ""

    # Fallback: update any remaining section references that point to known anchors
    for h in root.findall(".//w:hyperlink", NS):
        anchor = h.get(f"{{{NS['w']}}}anchor")
        if not anchor or anchor not in section_label_map:
            continue
        label = section_label_map[anchor]
        current = "".join(t.text for t in h.findall(".//w:t", NS) if t.text)
        if not current:
            continue
        if "section" not in current.lower():
            continue
        # ensure the numeric part is replaced even if runs are split oddly
        new_text = f"Section {label}"
        runs = h.findall(".//w:r", NS)
        if not runs:
            continue
        first = True
        for r in runs:
            t = r.find("w:t", NS)
            if t is None:
                continue
            if first:
                t.text = new_text
                first = False
            else:
                t.text = ""

    if debug:
        eprint(f"Caption paragraphs updated: {caption_patched}")

    if citation_scopes:
        removed = _filter_docx_bibliographies(root, citation_scopes, debug=debug)
        patched += removed
        if debug:
            eprint(f"Bibliography entries removed: {removed}")

    return ET.tostring(root, encoding="utf-8", xml_declaration=True), patched


def patch_docx(input_path: Path, output_path: Path | None = None, debug: bool = False) -> int:
    if output_path is None:
        output_path = input_path

    if not input_path.exists() or input_path.suffix.lower() != ".docx":
        return 0

    with zipfile.ZipFile(input_path) as z:
        if "word/document.xml" not in z.namelist():
            return 0
        document_xml = z.read("word/document.xml")
        other_files = {name: z.read(name) for name in z.namelist() if name != "word/document.xml"}

    project_dir = Path.cwd()
    if input_path.parent.name == "_book":
        project_dir = input_path.parent.parent
    citation_scopes = _citation_scopes(project_dir)

    new_xml, patched = _patch_document_xml(
        document_xml,
        citation_scopes=citation_scopes,
        debug=debug,
    )

    styles_xml = other_files.get("word/styles.xml")
    numbering_xml = other_files.get("word/numbering.xml")
    if styles_xml:
        new_styles_xml, new_numbering_xml, style_patched = _strip_main_heading_numbering(
            styles_xml,
            numbering_xml,
            debug=debug,
        )
        other_files["word/styles.xml"] = new_styles_xml
        if new_numbering_xml is not None:
            other_files["word/numbering.xml"] = new_numbering_xml
        patched += style_patched

    # Write out new docx
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "patched.docx"
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("word/document.xml", new_xml)
            for name, data in other_files.items():
                z.writestr(name, data)
        shutil.copyfile(tmp, output_path)

    return patched


def main() -> int:
    args = [a for a in sys.argv[1:] if a != "--debug"]
    debug = "--debug" in sys.argv

    if len(args) == 0:
        targets = _read_quarto_output_files()
        if debug:
            eprint("Targets:", targets)
        total = 0
        for p in targets:
            if p.suffix.lower() != ".docx":
                continue
            total += patch_docx(p, debug=debug)
        if debug:
            eprint(f"Patched {total} paragraphs")
        return 0

    if len(args) not in (1, 2):
        eprint("Usage: patch_appendix_docx.py input.docx [output.docx] [--debug]")
        return 2

    inp = Path(args[0])
    out = Path(args[1]) if len(args) == 2 else None
    patched = patch_docx(inp, out, debug=debug)
    if debug:
        eprint(f"Patched {patched} paragraphs")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
