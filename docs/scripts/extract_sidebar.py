#!/usr/bin/env python3
"""Extract the sidebar JSON block from a DocumenterVitepress config.mts
and prefix all link paths with the given locale.

Usage:
    python3 scripts/extract_sidebar.py <config.mts> <locale>

Prints the sidebar array (e.g. [{ text: 'Home', link: '/en/index' }, ...])
"""

import re, sys

with open(sys.argv[1]) as f:
    content = f.read()

idx = content.find("sidebar: [")
if idx < 0:
    print("ERROR: sidebar not found", file=sys.stderr)
    sys.exit(1)

start = idx + len("sidebar: ")
depth = 0
end = start
for i in range(start, len(content)):
    if content[i] == "[":
        depth += 1
    elif content[i] == "]":
        depth -= 1
        if depth == 0:
            end = i + 1
            break

block = content[start:end]
locale = sys.argv[2]

block = re.sub(
    r"link: '/(?!en/|cn/|https?)([^']+)'",
    f"link: '/{locale}/\\1'",
    block,
)

print(block)
