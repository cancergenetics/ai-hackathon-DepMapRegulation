"""Wrap docs/site/report.html (an Artifact-style page fragment) into a standalone docs/site/index.html.

report.html deliberately has no <!doctype>, <html>, <head> or <body> tags: it starts with <title>, <meta>,
<link> and <style> elements and then the page content. This script moves that leading block into <head>,
adds a charset and viewport meta, and puts the rest into <body>, so the page opens directly from a checkout.

Usage:  python scripts/figures/build_site.py
"""
from __future__ import annotations

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "docs" / "site" / "report.html"
DST = ROOT / "docs" / "site" / "index.html"

# Leading head-type elements: <title>…</title>, <meta …>, <link …>, <style>…</style>, comments, whitespace.
HEAD_ITEM = re.compile(
    r"\s*(?:<title>.*?</title>|<meta\b[^>]*>|<link\b[^>]*>|<style\b[^>]*>.*?</style>|<!--.*?-->)",
    re.S | re.I,
)


def split_head(page: str) -> tuple[str, str]:
    pos = 0
    while True:
        m = HEAD_ITEM.match(page, pos)
        if not m:
            break
        pos = m.end()
    return page[:pos].strip(), page[pos:].lstrip()


def build() -> Path:
    page = SRC.read_text(encoding="utf-8")
    if re.search(r"<!doctype|<html\b|<body\b", page[:2000], re.I):
        raise SystemExit("report.html already contains a document skeleton; refusing to wrap it twice")
    head, body = split_head(page)
    if "<title>" not in head:
        raise SystemExit("no leading <title> found in report.html")
    out = (
        "<!doctype html>\n"
        '<html lang="en">\n<head>\n'
        '<meta charset="utf-8">\n'
        '<meta name="viewport" content="width=device-width, initial-scale=1, viewport-fit=cover">\n'
        f"{head}\n"
        "<style>:root{padding-top:env(safe-area-inset-top,0px);padding-bottom:env(safe-area-inset-bottom,0px)}</style>\n"
        "</head>\n<body>\n"
        f"{body}\n"
        "</body>\n</html>\n"
    )
    DST.write_text(out, encoding="utf-8")
    return DST


if __name__ == "__main__":
    p = build()
    print(f"wrote {p.relative_to(ROOT)} ({p.stat().st_size:,} bytes)")
