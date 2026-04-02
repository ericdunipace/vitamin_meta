-- Mark appendix headings for DOCX post-processing.
-- Boundary: first level-1 heading with plain text "Appendix" (case-insensitive).
-- After boundary:
--   level-2 headings get hidden marker APPENDIX_H2
--   level-3 headings get hidden marker APPENDIX_H3
--   level-4 headings get hidden marker APPENDIX_H4 (optional)
-- The markers are inserted as hidden OpenXML runs so they survive into DOCX XML
-- but do not appear to the reader. A post-render script will swap paragraph styles
-- based on these markers and remove them.

local in_appendix = false

local function is_appendix_boundary(header)
  if header.level ~= 1 then return false end
  local text = pandoc.utils.stringify(header):gsub('%s+', ' '):gsub('^%s+', ''):gsub('%s+$', '')
  return text:lower() == 'appendix'
end

local function hidden_marker_run(marker)
  -- Only used for DOCX output.
  local xml = string.format(
    '<w:r><w:rPr><w:vanish/><w:specVanish/></w:rPr><w:t>%s</w:t></w:r>',
    marker
  )
  return pandoc.RawInline('openxml', xml)
end

function Header(header)
  if is_appendix_boundary(header) then
    in_appendix = true
    -- Ensure the boundary heading stays unnumbered.
    header.classes:insert('unnumbered')
    return header
  end

  if not in_appendix then
    return nil
  end

  if not FORMAT:match('docx') then
    return nil
  end

  local marker = nil
  if header.level == 2 then
    marker = 'APPENDIX_H2'
  elseif header.level == 3 then
    marker = 'APPENDIX_H3'
  elseif header.level == 4 then
    marker = 'APPENDIX_H4'
  elseif header.level == 5 then
    marker = 'APPENDIX_H5'
  end

  if marker then
    table.insert(header.content, 1, hidden_marker_run(marker))
    return header
  end

  return nil
end
