#!/usr/bin/env python3
"""Replace locale placeholders in config.mts with extracted sidebar/nav data.

Placeholders replaced:
  nav: 'REPLACE_ME_EN_NAV'       -> nav: <en_sidebar>
  'REPLACE_ME_EN_SIDEBAR'        -> <en_sidebar>
  nav: 'REPLACE_ME_CN_NAV'       -> nav: <cn_sidebar>
  'REPLACE_ME_CN_SIDEBAR'        -> <cn_sidebar>

Root themeConfig uses marker comments to locate nav/sidebar for replacement
(DocumenterVitepress populates them with un-prefixed links which we override):
  // REPLACE_ME_ROOT_NAV_START ... // REPLACE_ME_ROOT_NAV_END
  The first sidebar: [...] after the marker pair.

Usage:
    python3 scripts/inject_config.py <config.mts> <en_sidebar> <cn_sidebar>

Note: sidebar and nav are identical in DocumenterVitepress output.
"""

import re
import sys


def find_bracket_end(text, start):
    """Find the index after the closing ] that matches the [ at start."""
    depth = 0
    for i in range(start, len(text)):
        if text[i] == "[":
            depth += 1
        elif text[i] == "]":
            depth -= 1
            if depth == 0:
                return i + 1
    return len(text)


config_file = sys.argv[1]
en_sidebar = sys.argv[2]
cn_sidebar = sys.argv[3]

with open(config_file) as f:
    config = f.read()

# Locale nav/sidebar
config = config.replace("nav: 'REPLACE_ME_EN_NAV'", f"nav: {en_sidebar}")
config = config.replace("'REPLACE_ME_EN_SIDEBAR'", en_sidebar)
config = config.replace("nav: 'REPLACE_ME_CN_NAV'", f"nav: {cn_sidebar}")
config = config.replace("'REPLACE_ME_CN_SIDEBAR'", cn_sidebar)

# Root themeConfig — replace nav between markers with EN-prefixed data
config = re.sub(
    r"// REPLACE_ME_ROOT_NAV_START\n.*?// REPLACE_ME_ROOT_NAV_END",
    f"// REPLACE_ME_ROOT_NAV_START\n    nav: {en_sidebar},\n    // REPLACE_ME_ROOT_NAV_END",
    config,
    count=1,
    flags=re.DOTALL,
)

# Root themeConfig — replace the sidebar array after the nav marker
nav_end = config.find("REPLACE_ME_ROOT_NAV_END")
if nav_end >= 0:
    rest = config[nav_end:]
    m = re.search(r"sidebar: \[", rest)
    if m:
        arr_start = nav_end + m.end() - 1  # index of '['
        arr_end = find_bracket_end(config, arr_start)
        config = (
            config[: nav_end + m.start()] + f"sidebar: {en_sidebar}" + config[arr_end:]
        )

with open(config_file, "w") as f:
    f.write(config)
