#!/usr/bin/env python3
"""Create a DOCX reference file with appendix heading styles.

Usage:
  build_reference_doc.py base.docx reference/custom-reference.docx
"""
from __future__ import annotations

import sys
import zipfile
import tempfile
import shutil
from pathlib import Path
import xml.etree.ElementTree as ET

NS = {"w": "http://schemas.openxmlformats.org/wordprocessingml/2006/main"}
ET.register_namespace("w", NS["w"])


def build_reference(base_path: Path, output_path: Path) -> None:
    with zipfile.ZipFile(base_path) as z:
        styles_xml = z.read("word/styles.xml")
        numbering_xml = z.read("word/numbering.xml")
        settings_xml = z.read("word/settings.xml") if "word/settings.xml" in z.namelist() else None
        other_files = {name: z.read(name) for name in z.namelist() if name not in ("word/styles.xml", "word/numbering.xml", "word/settings.xml")}

    styles = ET.fromstring(styles_xml)
    numbering = ET.fromstring(numbering_xml)

    # find max ids
    max_abs = -1
    for absnum in numbering.findall("w:abstractNum", NS):
        val = absnum.get(f"{{{NS['w']}}}abstractNumId")
        if val is not None:
            max_abs = max(max_abs, int(val))

    max_num = -1
    for num in numbering.findall("w:num", NS):
        val = num.get(f"{{{NS['w']}}}numId")
        if val is not None:
            max_num = max(max_num, int(val))

    abs_id = max_abs + 1
    num_id = max_num + 1

    # Build appendix numbering (A, A.1, A.1.1)
    abs_el = ET.Element(f"{{{NS['w']}}}abstractNum", {f"{{{NS['w']}}}abstractNumId": str(abs_id)})

    def add_level(ilvl: int, num_fmt: str, lvl_text: str) -> None:
        lvl = ET.SubElement(abs_el, f"{{{NS['w']}}}lvl", {f"{{{NS['w']}}}ilvl": str(ilvl)})
        ET.SubElement(lvl, f"{{{NS['w']}}}start", {f"{{{NS['w']}}}val": "1"})
        ET.SubElement(lvl, f"{{{NS['w']}}}numFmt", {f"{{{NS['w']}}}val": num_fmt})
        ET.SubElement(lvl, f"{{{NS['w']}}}lvlText", {f"{{{NS['w']}}}val": lvl_text})
        ET.SubElement(lvl, f"{{{NS['w']}}}lvlJc", {f"{{{NS['w']}}}val": "left"})

    add_level(0, "upperLetter", "%1")
    add_level(1, "decimal", "%1.%2")
    add_level(2, "decimal", "%1.%2.%3")
    add_level(3, "decimal", "%1.%2.%3.%4")

    numbering.append(abs_el)
    num_el = ET.Element(f"{{{NS['w']}}}num", {f"{{{NS['w']}}}numId": str(num_id)})
    ET.SubElement(num_el, f"{{{NS['w']}}}abstractNumId", {f"{{{NS['w']}}}val": str(abs_id)})
    numbering.append(num_el)

    # Remove existing appendix styles if present
    for st in list(styles.findall("w:style", NS)):
        sid = st.get(f"{{{NS['w']}}}styleId")
        if sid in ("AppendixHeading1", "AppendixHeading2", "AppendixHeading3", "AppendixHeading4"):
            styles.remove(st)

    def add_style(style_id: str, name: str, based_on: str, ilvl: int) -> None:
        st = ET.Element(
            f"{{{NS['w']}}}style",
            {f"{{{NS['w']}}}type": "paragraph", f"{{{NS['w']}}}styleId": style_id, f"{{{NS['w']}}}customStyle": "1"},
        )
        ET.SubElement(st, f"{{{NS['w']}}}name", {f"{{{NS['w']}}}val": name})
        ET.SubElement(st, f"{{{NS['w']}}}basedOn", {f"{{{NS['w']}}}val": based_on})
        ET.SubElement(st, f"{{{NS['w']}}}qFormat")
        ppr = ET.SubElement(st, f"{{{NS['w']}}}pPr")
        numpr = ET.SubElement(ppr, f"{{{NS['w']}}}numPr")
        ET.SubElement(numpr, f"{{{NS['w']}}}ilvl", {f"{{{NS['w']}}}val": str(ilvl)})
        ET.SubElement(numpr, f"{{{NS['w']}}}numId", {f"{{{NS['w']}}}val": str(num_id)})
        styles.append(st)

    add_style("AppendixHeading1", "Appendix Heading 1", "Heading2", 0)
    add_style("AppendixHeading2", "Appendix Heading 2", "Heading3", 1)
    add_style("AppendixHeading3", "Appendix Heading 3", "Heading4", 2)
    add_style("AppendixHeading4", "Appendix Heading 4", "Heading5", 3)

    # Update caption numbering to use appendix chapter style
    if settings_xml is not None:
        settings = ET.fromstring(settings_xml)
        captions = settings.find("w:captions", NS)
        if captions is None:
            captions = ET.SubElement(settings, f"{{{NS['w']}}}captions")

        def upsert_caption(name: str) -> None:
            cap = None
            for c in captions.findall("w:caption", NS):
                if c.get(f"{{{NS['w']}}}name") == name:
                    cap = c
                    break
            if cap is None:
                cap = ET.SubElement(captions, f"{{{NS['w']}}}caption", {f"{{{NS['w']}}}name": name})
            for child in list(cap):
                cap.remove(child)
            ET.SubElement(cap, f"{{{NS['w']}}}autoCapt", {f"{{{NS['w']}}}val": "0"})
            ET.SubElement(cap, f"{{{NS['w']}}}numFmt", {f"{{{NS['w']}}}val": "decimal"})
            ET.SubElement(cap, f"{{{NS['w']}}}chapStyle", {f"{{{NS['w']}}}val": "AppendixHeading1"})
            ET.SubElement(cap, f"{{{NS['w']}}}chapSep", {f"{{{NS['w']}}}val": "period"})
            ET.SubElement(cap, f"{{{NS['w']}}}includeChapNum")

        upsert_caption("Table")
        upsert_caption("Figure")
    else:
        settings = None

    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "custom-reference.docx"
        with zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as z:
            z.writestr("word/styles.xml", ET.tostring(styles, encoding="utf-8", xml_declaration=True))
            z.writestr("word/numbering.xml", ET.tostring(numbering, encoding="utf-8", xml_declaration=True))
            if settings is not None:
                z.writestr("word/settings.xml", ET.tostring(settings, encoding="utf-8", xml_declaration=True))
            for name, data in other_files.items():
                z.writestr(name, data)
        shutil.copyfile(tmp, output_path)


def main() -> int:
    if len(sys.argv) != 3:
        print("Usage: build_reference_doc.py base.docx reference/custom-reference.docx", file=sys.stderr)
        return 2
    base = Path(sys.argv[1])
    out = Path(sys.argv[2])
    if not base.exists():
        print(f"Base DOCX not found: {base}", file=sys.stderr)
        return 1
    out.parent.mkdir(parents=True, exist_ok=True)
    build_reference(base, out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
