#!/usr/bin/env bash
# ==============================================================================
# make_combined.sh — Build a single multilingual VitePress site (EN + CN)
# ==============================================================================
#
#   EN at /en/  ·  CN at /cn/  ·  / redirects to /en/
#
# Usage:
#   bash docs/make_combined.sh
#   DOCUMENTER_DRAFT=true bash docs/make_combined.sh
#   python3 -m http.server 8080 --directory docs/build_combined/1/ 
#.  #or
#.  
# Prerequisites:
#   julia --project=docs -e 'import Pkg; Pkg.instantiate()'
#   node + npm available on PATH
# ==============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRIPTS="$SCRIPT_DIR/scripts"
DRAFT="${DOCUMENTER_DRAFT:-false}"
COMBINED="$SCRIPT_DIR/build_combined"

# ── Phase 0: Copy assets into source dirs ──────────────────────────────────
# DocumenterVitepress looks for assets at joinpath(sourcedir, "assets").
# With source = "src/en" / "src/cn", it expects docs/src/{en,cn}/assets/.
# The real files live in docs/src/assets/ — copy them before building.

SRC_ASSETS="$SCRIPT_DIR/src/assets"
for locale in en cn; do
    LOCALE_ASSETS="$SCRIPT_DIR/src/$locale/assets"
    mkdir -p "$LOCALE_ASSETS"
    cp "$SRC_ASSETS"/logo.svg "$SRC_ASSETS"/favicon.ico "$LOCALE_ASSETS/"
done

# ── Phase 1+2: Build EN and CN markdown ──────────────────────────────────

echo "=== Phase 1: Building EN markdown ==="
SKIP_VITEPRESS=true DOCUMENTER_DRAFT="$DRAFT" \
    julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/make.jl"

echo "=== Phase 2: Building CN markdown ==="
SKIP_VITEPRESS=true DOCUMENTER_DRAFT="$DRAFT" \
    julia --project="$SCRIPT_DIR" "$SCRIPT_DIR/make_cn.jl"

# ── Phase 3: Merge into build_combined/ ───────────────────────────────────

echo "=== Phase 3: Merging into build_combined/ ==="

EN="$SCRIPT_DIR/build_en/.documenter"
CN="$SCRIPT_DIR/build_cn/.documenter"
MD="$COMBINED/.documenter"

rm -rf "$COMBINED"
mkdir -p "$MD"

# 3a. Copy EN build as the base (has .vitepress config, public/, components/)
cp -a "$EN/." "$MD/"

# 3b. Move EN content into en/
mkdir -p "$MD/en" "$MD/cn"
for item in "$MD"/*; do
    base="$(basename "$item")"
    [[ "$base" == .vitepress || "$base" == public || "$base" == components \
       || "$base" == en || "$base" == cn || "$base" == assets ]] && continue
    mv "$item" "$MD/en/"
done

# 3c. Copy CN markdown into cn/ (skip .vitepress/ and components/)
cd "$CN"
find . -not -path './.vitepress*' -not -path './components*' -not -name '.' \
    -print0 | while IFS= read -r -d '' f; do
    if [ -d "$CN/$f" ]; then
        mkdir -p "$MD/cn/$f"
    else
        mkdir -p "$MD/cn/$(dirname "$f")"
        cp "$CN/$f" "$MD/cn/$f"
    fi
done

# 3d. Root redirect page
mkdir -p "$MD/public"
cp "$SCRIPTS/redirect.html" "$MD/public/index.html"

# ── Phase 3.5: Fix internal links in merged markdown ──────────────────────

echo "=== Phase 3.5: Fixing internal links ==="
python3 "$SCRIPTS/fix_links.py" "$MD/en" "$MD/cn"

# ── Phase 4: Inject locale data into config.mts ───────────────────────────

echo "=== Phase 4: Injecting locale data ==="

CONFIG="$MD/.vitepress/config.mts"
CN_CONFIG="$CN/.vitepress/config.mts"

EN_SIDEBAR="$(python3 "$SCRIPTS/extract_sidebar.py" "$CONFIG" en)"
CN_SIDEBAR="$(python3 "$SCRIPTS/extract_sidebar.py" "$CN_CONFIG" cn)"

python3 "$SCRIPTS/inject_config.py" "$CONFIG" "$EN_SIDEBAR" "$CN_SIDEBAR"
echo "Config injected."

# ── Phase 5: Run VitePress build ──────────────────────────────────────────

echo "=== Phase 5: Running VitePress build ==="

cd "$SCRIPT_DIR"

DV_PKG="$(julia --project="$SCRIPT_DIR" -e 'using DocumenterVitepress; println(dirname(dirname(pathof(DocumenterVitepress))))')"
DV_TEMPLATE="$DV_PKG/template/package.json"

COPIED_PKG_JSON=false
if [ ! -f package.json ]; then
    cp "$DV_TEMPLATE" package.json
    COPIED_PKG_JSON=true
fi
npm install
npx vitepress build "$MD"

if [ "$COPIED_PKG_JSON" = true ]; then
    rm -f package.json package-lock.json
fi

OUTDIR="$COMBINED/1"
if [ -d "$OUTDIR" ] && [ "$(ls -A "$OUTDIR" 2>/dev/null)" ]; then
    echo 'var DOCUMENTER_CURRENT_VERSION = "dev";' > "$OUTDIR/siteinfo.js"
    echo "=== Build complete ==="
    echo "Output: $OUTDIR"
    echo "EN at /en/  ·  CN at /cn/  ·  Redirect at /"
else
    echo "WARNING: VitePress output not found at $OUTDIR"
    exit 1
fi
