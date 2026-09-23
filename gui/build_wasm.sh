#!/bin/bash
# Build script for RustFEA GUI WASM target
# Requires: wasm-bindgen-cli, wasm-opt (optional)

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Building RustFEA GUI for WASM..."

# Build in release mode for smaller output
cargo build --target wasm32-unknown-unknown --features wasm --no-default-features --release

echo "Running wasm-bindgen..."

# Generate JS bindings directly into web/ folder (alongside index.html)
wasm-bindgen \
    --target web \
    --out-dir web \
    --no-typescript \
    ../target/wasm32-unknown-unknown/release/rust_fea_gui.wasm

# Optional: Optimize WASM (requires wasm-opt from binaryen)
if command -v wasm-opt &> /dev/null; then
    echo "Optimizing WASM with wasm-opt..."
    wasm-opt -Oz web/rust_fea_gui_bg.wasm -o web/rust_fea_gui_bg.wasm
fi

echo ""
echo "Build complete! Output in web/"
echo ""
echo "To serve locally, run:"
echo "  cd web && python3 -m http.server 8080"
echo "Then open http://localhost:8080"
