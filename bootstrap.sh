#!/bin/sh
set -e  # Exit immediately if a command fails

echo "Running autoreconf to regenerate build system files..."
autoreconf -fi

echo "Bootstrap complete. Now run:"
echo "  ./configure [options]"
echo "  make"
