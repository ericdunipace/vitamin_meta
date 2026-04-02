# Quarto → Word Appendix Numbering (Automatic)

This project implements an automatic appendix numbering workflow for DOCX output. The pipeline uses:

- Quarto `reference-doc` to define custom Word styles and numbering.
- A minimal Pandoc Lua filter to detect the appendix boundary and insert hidden markers.
- A post-render DOCX patcher that swaps appendix headings to the custom styles.

No manual Word edits are required after each render.

## Source Pattern

Use this heading structure in your `.qmd` files:

```markdown
# Main Section
## Main subsection

# Appendix {.unnumbered}
## Supplementary Tables
### Table derivations
## Additional Analyses
### Sensitivity analyses
```

Rules:
- The first level-1 heading with text exactly `Appendix` (case-insensitive) is the boundary.
- The boundary heading itself remains unnumbered.
- After that, level-2 headings become `Appendix Heading 1`, level-3 become `Appendix Heading 2`.
- Level-4 headings become `Appendix Heading 3`, level-5 headings become `Appendix Heading 4`.

## How It Works

1. `filters/appendix-boundary.lua` detects the Appendix boundary and injects a hidden marker into subsequent appendix headings (DOCX only).
2. `scripts/patch_appendix_docx.py` runs after render, finds those markers in `word/document.xml`, swaps the paragraph style to a custom appendix heading style, and removes the markers.
3. `reference/custom-reference.docx` contains Word styles and numbering that define `A`, `A.1`, `A.1.1`, etc.

## Render the Demo

From this directory (recommended so project hooks run):

```bash
quarto render templates/appendix-demo.qmd --to docx
```

Or from elsewhere using the project root:

```bash
quarto render /Users/eifer/GoogleDrive/Research/vitamin_meta/outputs/documents
```

This will automatically patch the generated DOCX. In Word you should see:

- `Appendix` unnumbered
- `A Supplementary Tables`
- `A.1 Table derivations`
- `B Additional Analyses`
- `B.1 Sensitivity analyses`

## Custom Word Styles

The reference document defines these custom paragraph styles:

- `Appendix Heading 1`
- `Appendix Heading 2`
- `Appendix Heading 3`
- `Appendix Heading 4`

These are linked to a multilevel list:

- `Appendix Heading 1` → `A, B, C ...`
- `Appendix Heading 2` → `A.1, A.2 ...`
- `Appendix Heading 3` → `A.1.1 ...`
- `Appendix Heading 4` → `A.1.1.1 ...`

## Regenerate `custom-reference.docx`

If you need to rebuild the reference doc from a fresh Quarto/Pandoc base:

1. Create a base DOCX (any Quarto render works):

```bash
quarto render templates/appendix-demo.qmd --to docx
```

2. Build the reference doc from that base:

```bash
python scripts/build_reference_doc.py templates/appendix-demo.docx reference/custom-reference.docx
```

## Troubleshooting

- If appendix headings stay numeric: confirm the post-render hook ran and style names match exactly.
- If numbering is missing for appendix styles: the custom styles in `reference/custom-reference.docx` are not linked to a multilevel list.
- If figure/table numbering in the appendix should include chapter letters (A.1, A.2): ensure caption numbering in the reference doc uses `Appendix Heading 1` as the chapter style.
- If nothing happens: ensure the `Appendix` boundary heading exists and is level-1.

## Files

- `_quarto.yml`: DOCX reference doc, Lua filter, and post-render hook.
- `filters/appendix-boundary.lua`: markers for appendix headings.
- `scripts/patch_appendix_docx.py`: DOCX patcher.
- `reference/custom-reference.docx`: Word styles + numbering.
- `templates/appendix-demo.qmd`: minimal demo.
