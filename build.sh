#!/bin/bash
# Build script for Azure deployment
# This builds the Next.js background before deploying Streamlit

echo "Building Next.js background..."
cd my-app
npm install
npm run build
cd ..

echo "✓ Build complete! Ready for deployment."
echo "Next.js static files are in: my-app/embed/"