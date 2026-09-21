#!/bin/bash
# Build script for RustFEA GUI WASM target
# Requires: wasm-bindgen-cli, wasm-opt (optional)

set -e

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd "$SCRIPT_DIR"

echo "Building RustFEA GUI for WASM..."

# Build in release mode for smaller output
cargo build --target wasm32-unknown-unknown --features web --no-default-features --release

echo "Running wasm-bindgen..."

# Create output directory
mkdir -p web/dist

# Generate JS bindings
wasm-bindgen \
    --target web \
    --out-dir web/dist \
    ../target/wasm32-unknown-unknown/release/rust_fea_gui.wasm

# Optional: Optimize WASM (requires wasm-opt from binaryen)
if command -v wasm-opt &> /dev/null; then
    echo "Optimizing WASM with wasm-opt..."
    wasm-opt -Oz web/dist/rust_fea_gui_bg.wasm -o web/dist/rust_fea_gui_bg.wasm
fi

# Copy index.html
cp web/index.html web/dist/

echo ""
echo "Build complete! Output in web/dist/"
echo ""
echo "To serve locally, run:"
echo "  cd web/dist && python3 -m http.server 8080"
echo "Then open http://localhost:8080"
