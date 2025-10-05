#!/bin/bash
# Azure Web App startup script - runs both Next.js and Streamlit

set -e  # Exit on error

echo "=== Starting Impactor-2025 deployment ==="

# Install Python dependencies
echo "Installing Python dependencies..."
pip install --upgrade pip
pip install -r requirements.txt

# Install Node.js dependencies and build Next.js
echo "Installing Node.js dependencies..."
cd my-app
npm install --legacy-peer-deps --production=false
npm run build

# Start Next.js on Azure's port (must be $PORT)
echo "Starting Next.js server on port ${PORT:-8000}..."
PORT=${PORT:-8000} npm run start &
NEXTJS_PID=$!
cd ..

# Wait for Next.js to be ready
echo "Waiting for Next.js to start..."
for i in {1..30}; do
    if curl -s http://localhost:${PORT:-8000} > /dev/null; then
        echo "Next.js is ready!"
        break
    fi
    echo "Waiting... ($i/30)"
    sleep 2
done

# Start Streamlit in background on a different port (not exposed)
echo "Starting Streamlit on port 8501..."
streamlit run main.py \
    --server.port=8501 \
    --server.address=0.0.0.0 \
    --server.headless=true \
    --server.enableCORS=false \
    --server.enableXsrfProtection=false &

# Keep both running
wait -n