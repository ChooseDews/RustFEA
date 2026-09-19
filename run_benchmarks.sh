#!/bin/bash
# Run FEA benchmark suite and generate report
# Usage: ./run_benchmarks.sh [options]

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# Default values
OUTPUT_DIR="examples/output"
OUTPUT_FORMAT="markdown"
OUTPUT_FILE=""
BENCHMARKS=""
VERBOSE=""

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        -o|--output)
            OUTPUT_FILE="$2"
            shift 2
            ;;
        -f|--format)
            OUTPUT_FORMAT="$2"
            shift 2
            ;;
        -b|--benchmarks)
            BENCHMARKS="$2"
            shift 2
            ;;
        -v|--verbose)
            VERBOSE="1"
            shift
            ;;
        -h|--help)
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  -o, --output FILE    Output file path (default: stdout)"
            echo "  -f, --format FORMAT  Output format: markdown, json, csv (default: markdown)"
            echo "  -b, --benchmarks     Comma-separated list of benchmarks to run"
            echo "  -v, --verbose        Enable verbose output"
            echo "  -h, --help           Show this help message"
            echo ""
            echo "Examples:"
            echo "  $0                                    # Run all benchmarks"
            echo "  $0 -o report.md                       # Save to file"
            echo "  $0 -b uniaxial,shear -f json          # Run specific benchmarks"
            echo "  $0 -o report.csv -f csv               # Export as CSV"
            exit 0
            ;;
        *)
            echo "Unknown option: $1"
            exit 1
            ;;
    esac
done

# Build if needed
echo "Building benchmark runner..."
cargo build --bin run_benchmarks --quiet

# Construct command
CMD="cargo run --bin run_benchmarks --quiet --"
CMD="$CMD --format $OUTPUT_FORMAT"

if [ -n "$OUTPUT_FILE" ]; then
    CMD="$CMD --output $OUTPUT_FILE"
fi

if [ -n "$BENCHMARKS" ]; then
    CMD="$CMD --benchmarks $BENCHMARKS"
fi

if [ -n "$VERBOSE" ]; then
    CMD="$CMD --verbose"
    export RUST_LOG=info
fi

# Run benchmarks
echo "Running FEA benchmarks..."
echo ""
eval $CMD

if [ -n "$OUTPUT_FILE" ]; then
    echo ""
    echo "Report saved to: $OUTPUT_FILE"
fi
