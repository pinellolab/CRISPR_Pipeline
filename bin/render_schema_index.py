#!/usr/bin/env python3
"""Render the pipeline parameter schema as docs/index.html.

This keeps the nf-core-style parameter page on a stable filename so it can be
served directly by GitHub Pages or opened locally without an extra copy step.
"""

import argparse
import html
import json
import re
from pathlib import Path


def slugify(value):
    slug = re.sub(r"[^a-z0-9]+", "-", value.lower()).strip("-")
    return slug or "section"


def fmt_type(schema):
    value = schema.get("type", "-")
    if isinstance(value, list):
        return " | ".join(value)
    return str(value)


def fmt_default(schema):
    if "default" not in schema:
        return "-"
    value = schema["default"]
    if value is None:
        return "<code>null</code>"
    if isinstance(value, bool):
        return f"<code>{str(value).lower()}</code>"
    return f"<code>{html.escape(repr(value))}</code>"


def fmt_options(schema):
    if "enum" in schema:
        return ", ".join(f"<code>{html.escape(str(option))}</code>" for option in schema["enum"])

    schema_format = schema.get("format")
    if schema_format == "file-path":
        return "file path"
    if schema_format == "directory-path":
        return "directory path"
    if schema_format == "uri":
        return "URL"
    if schema.get("type") == "boolean":
        return "<code>true</code>, <code>false</code>"
    return "-"


def schema_groups(schema):
    defs = schema.get("$defs", {})
    refs = schema.get("allOf", [])

    for ref in refs:
        ref_path = ref.get("$ref", "")
        if not ref_path.startswith("#/$defs/"):
            continue
        key = ref_path.rsplit("/", 1)[-1]
        group = defs.get(key)
        if group:
            yield key, group


def render_row(name, param):
    description = html.escape(param.get("description", "-"))
    help_text = param.get("help_text")
    help_html = ""
    if help_text:
        help_html = (
            "<details><summary>More context</summary>"
            f"<small>{html.escape(help_text)}</small></details>"
        )

    return (
        "<tr>"
        f'<td class="param"><code>{html.escape(name)}</code></td>'
        f"<td><p>{description}</p>{help_html}</td>"
        f"<td>{html.escape(fmt_type(param))}</td>"
        f"<td>{fmt_default(param)}</td>"
        f"<td>{fmt_options(param)}</td>"
        "</tr>"
    )


def render_group(group):
    title = group.get("title", "Parameters")
    description = group.get("description", "")
    properties = group.get("properties", {})
    rows = "\n".join(render_row(name, param) for name, param in properties.items())

    return f"""
        <section class="card" id="{slugify(title)}">
          <div class="card-header">
            <div>
              <h2>{html.escape(title)}</h2>
              <p>{html.escape(description)}</p>
            </div>
            <span class="count">{len(properties)} params</span>
          </div>
          <div class="table-wrap">
            <table>
              <thead>
                <tr>
                  <th>Parameter</th>
                  <th>Description</th>
                  <th>Type</th>
                  <th>Default</th>
                  <th>Options</th>
                </tr>
              </thead>
              <tbody>
                {rows}
              </tbody>
            </table>
          </div>
        </section>
        """


def render_page(schema):
    title = schema.get("title", "Pipeline parameters")
    description = schema.get("description", "")
    groups = [group for _, group in schema_groups(schema)]
    nav = "".join(
        f'<a href="#{slugify(group.get("title", "Parameters"))}">{html.escape(group.get("title", "Parameters"))}</a>'
        for group in groups
    )
    sections = "\n".join(render_group(group) for group in groups)

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{html.escape(title)}</title>
  <style>
    :root {{
      --nf-green: #24b064;
      --nf-dark: #0f172a;
      --nf-muted: #64748b;
      --nf-line: #e2e8f0;
      --nf-bg: #f8fafc;
      --nf-card: #ffffff;
      --nf-code: #f1f5f9;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      background: var(--nf-bg);
      color: var(--nf-dark);
      font-family: Inter, ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
      line-height: 1.5;
    }}
    .hero {{
      background: linear-gradient(135deg, #0f172a 0%, #14532d 100%);
      color: white;
      padding: 44px 32px 38px;
    }}
    .hero h1 {{ margin: 0 0 10px; font-size: clamp(32px, 5vw, 54px); letter-spacing: -0.04em; }}
    .hero p {{ max-width: 960px; margin: 0; color: #d1fae5; font-size: 18px; }}
    .shell {{
      display: grid;
      grid-template-columns: 260px minmax(0, 1fr);
      gap: 28px;
      padding: 28px 32px 52px;
      max-width: 1500px;
      margin: 0 auto;
    }}
    aside {{
      position: sticky;
      top: 20px;
      align-self: start;
      background: var(--nf-card);
      border: 1px solid var(--nf-line);
      border-radius: 16px;
      padding: 18px;
      box-shadow: 0 10px 30px rgba(15, 23, 42, 0.06);
    }}
    aside h2 {{ margin: 0 0 12px; font-size: 15px; text-transform: uppercase; color: var(--nf-muted); letter-spacing: 0.08em; }}
    aside a {{
      display: block;
      padding: 9px 10px;
      border-radius: 10px;
      color: #166534;
      text-decoration: none;
      font-weight: 600;
      font-size: 14px;
    }}
    aside a:hover {{ background: #dcfce7; }}
    .toolbar {{ display: flex; gap: 12px; align-items: center; margin-bottom: 18px; }}
    .search {{
      width: 100%;
      padding: 13px 16px;
      border: 1px solid var(--nf-line);
      border-radius: 14px;
      background: white;
      font-size: 15px;
      box-shadow: 0 8px 25px rgba(15, 23, 42, 0.05);
    }}
    .card {{
      background: var(--nf-card);
      border: 1px solid var(--nf-line);
      border-radius: 18px;
      margin-bottom: 24px;
      overflow: hidden;
      box-shadow: 0 14px 35px rgba(15, 23, 42, 0.06);
    }}
    .card-header {{
      display: flex;
      justify-content: space-between;
      gap: 16px;
      padding: 22px 24px;
      border-bottom: 1px solid var(--nf-line);
      background: linear-gradient(180deg, #ffffff 0%, #f8fafc 100%);
    }}
    .card h2 {{ margin: 0; font-size: 24px; letter-spacing: -0.02em; }}
    .card-header p {{ margin: 6px 0 0; color: var(--nf-muted); }}
    .count {{
      align-self: start;
      white-space: nowrap;
      background: #dcfce7;
      color: #166534;
      border-radius: 999px;
      padding: 6px 12px;
      font-size: 13px;
      font-weight: 700;
    }}
    .table-wrap {{ overflow-x: auto; }}
    table {{ width: 100%; border-collapse: collapse; min-width: 980px; }}
    th, td {{ padding: 13px 16px; text-align: left; vertical-align: top; border-bottom: 1px solid var(--nf-line); }}
    th {{ background: #f8fafc; color: #334155; font-size: 12px; text-transform: uppercase; letter-spacing: 0.08em; }}
    tr:hover td {{ background: #f8fafc; }}
    td.param {{ width: 260px; }}
    code {{
      background: var(--nf-code);
      color: #0f172a;
      border: 1px solid #e2e8f0;
      border-radius: 6px;
      padding: 2px 6px;
      font-size: 13px;
      white-space: nowrap;
    }}
    td p {{ margin: 0; }}
    details {{ margin-top: 7px; color: var(--nf-muted); }}
    summary {{ cursor: pointer; color: #166534; font-weight: 700; }}
    .hidden {{ display: none; }}
    @media (max-width: 900px) {{
      .shell {{ grid-template-columns: 1fr; padding: 20px; }}
      aside {{ position: static; }}
      .hero {{ padding: 34px 20px; }}
    }}
  </style>
</head>
<body>
  <header class="hero">
    <h1>{html.escape(title)}</h1>
    <p>{html.escape(description)}</p>
  </header>
  <div class="shell">
    <aside>
      <h2>Parameter groups</h2>
      {nav}
    </aside>
    <main>
      <div class="toolbar">
        <input id="search" class="search" type="search" placeholder="Search parameters, descriptions, defaults, or options...">
      </div>
      {sections}
    </main>
  </div>
  <script>
    const search = document.getElementById('search');
    search.addEventListener('input', () => {{
      const query = search.value.toLowerCase();
      document.querySelectorAll('tbody tr').forEach(row => {{
        row.classList.toggle('hidden', !row.innerText.toLowerCase().includes(query));
      }});
      document.querySelectorAll('.card').forEach(card => {{
        const visibleRows = card.querySelectorAll('tbody tr:not(.hidden)').length;
        card.classList.toggle('hidden', visibleRows === 0);
      }});
    }});
  </script>
</body>
</html>
"""


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--schema", default="nextflow_schema.json", help="Path to nextflow_schema.json")
    parser.add_argument("--output", default="docs/index.html", help="Output HTML path")
    args = parser.parse_args()

    schema_path = Path(args.schema)
    output_path = Path(args.output)
    schema = json.loads(schema_path.read_text())

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(render_page(schema))
    print(f"Wrote {output_path}")


if __name__ == "__main__":
    main()
