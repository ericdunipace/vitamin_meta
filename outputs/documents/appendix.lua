function Header(el)
  if el.level == 1 and el.identifier == "appendix" then
    in_appendix = true
    return el
  end

  if in_appendix and el.level == 1 then
    el.attributes['custom-style'] = 'Appendix Heading 1'
  end

  return el
end