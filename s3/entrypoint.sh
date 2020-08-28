#!/bin/sh

echo "Preparing bucket..."
mkdir -p /data/${MINIO_BUCKET}

echo "Starting minio server..."
env MINIO_ACCESS_KEY=${MINIO_ACCESS_KEY} MINIO_SECRET_KEY=${MINIO_SECRET_KEY} minio server /data &
PID=$!
sleep 5

echo "Setting up client..."
mc config host add s3 http://localhost:9000 ${MINIO_ACCESS_KEY} ${MINIO_SECRET_KEY}

echo "Setting bucket policy..."
mc policy set public s3/${MINIO_BUCKET}

echo "Ready."
wait ${PID}
