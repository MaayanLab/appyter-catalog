#!/bin/sh

echo "Preparing bucket..."
mkdir -p /data/${MINIO_BUCKET}

echo "Starting minio server..."
env MINIO_ACCESS_KEY=${MINIO_ACCESS_KEY} MINIO_SECRET_KEY=${MINIO_SECRET_KEY} minio server /data &
PID=$!

let SUCCESS=1
let RETRIES=0
while [ ${SUCCESS} -ne 0 ]; do
  if [ "${RETRIES}" -ge 5 ]; then
    echo "Maximum retries to set bucket policy, quitting"
    exit 1
  fi
  sleep 5

  echo "Setting up client..."
  mc alias set s3 http://localhost:9000 ${MINIO_ACCESS_KEY} ${MINIO_SECRET_KEY}
  let SUCCESS=$?

  echo "Setting bucket policy..."
  mc policy set download s3/${MINIO_BUCKET}
  let SUCCESS=$?

  let RETRIES="${RETRIES}+1"
done

echo "Ready."
wait ${PID}
