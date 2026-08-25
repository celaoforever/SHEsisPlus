#!/usr/bin/env bash
set -e

# 1. Install Docker & Docker Compose plugin if missing
if ! command -v docker &>/dev/null; then
    echo "[*] Docker not found. Installing Docker..."
    curl -fsSL https://get.docker.com | sh
    systemctl enable --now docker
else
    echo "[*] Docker is already installed."
fi

# Ensure docker compose command is available
if ! docker compose version &>/dev/null; then
    echo "[*] Installing docker-compose-plugin..."
    apt-get update && apt-get install -y docker-compose-plugin
fi

# 2. Create project directory
PROJECT_DIR="/opt/shesisplus"
mkdir -p "$PROJECT_DIR"
cd "$PROJECT_DIR"

# 3. Write docker-compose.yml
echo "[*] Generating docker-compose.yml..."
cat << 'EOF' > docker-compose.yml
version: "3.8"

services:
  shesis-mongo:
    image: mongo:4.4
    container_name: shesis-mongo
    restart: unless-stopped
    volumes:
      - mongo_data:/data/db

  shesis-redis:
    image: redis:6-alpine
    container_name: shesis-redis
    restart: unless-stopped

  shesis-web:
    image: realsamy/shesisplus:latest
    container_name: shesis-web
    restart: unless-stopped
    ports:
      - "5903:5903"
    environment:
      - REDIS_HOST=shesis-redis
      - REDIS_PORT=6379
      - MONGO_HOST=shesis-mongo
      - MONGO_PORT=27017
    depends_on:
      - shesis-mongo
      - shesis-redis
    volumes:
      - shesis_results:/app/SHEsisWebServer/public/tmp

volumes:
  mongo_data:
  shesis_results:
EOF

# 4. Pull latest images and run
echo "[*] Pulling images and starting services..."
docker compose pull
docker compose up -d

echo "=========================================="
echo " SHEsisPlus is running at: http://$(curl -s ifconfig.me || echo 'localhost'):5903"
echo " Project location: $PROJECT_DIR"
echo " View logs: cd $PROJECT_DIR && docker compose logs -f"
echo "=========================================="
