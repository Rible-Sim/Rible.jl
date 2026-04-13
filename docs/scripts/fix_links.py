#!/usr/bin/env python3
"""Prefix Documenter-generated internal links with locale paths.

DocumenterVitepress renders @ref links as absolute paths (e.g. ](/page))
without any locale prefix. This script rewrites them to include /en/ or /cn/:

    ](/get_started)  →  ](/en/get_started)   (in EN files)
    ](/get_started)  →  ](/cn/get_started)   (in CN files)

Usage:
    python3 scripts/fix_links.py <en_dir> <cn_dir>
"""

import os, re, sys

INTERNAL_LINK = re.compile(r"\]\(/(?!en/|cn/|https?)([^)\s]*)")
YAML_LINK = re.compile(r"(link: )/(?!en/|cn/|https?)(\S+)")


def fix_links(locale_dir, locale_prefix):
    for root, _dirs, files in os.walk(locale_dir):
        for fname in files:
            if not fname.endswith(".md"):
                continue
            fpath = os.path.join(root, fname)
            with open(fpath) as f:
                content = f.read()
            rewritten = INTERNAL_LINK.sub(
                lambda m: f"](/{locale_prefix}/{m.group(1)}", content
            )
            rewritten = YAML_LINK.sub(
                lambda m: f"{m.group(1)}/{locale_prefix}/{m.group(2)}", rewritten
            )
            if rewritten != content:
                with open(fpath, "w") as f:
                    f.write(rewritten)


if __name__ == "__main__":
    fix_links(sys.argv[1], "en")
    fix_links(sys.argv[2], "cn")
