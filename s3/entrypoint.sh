#!/bin/sh

echo "Preparing bucket..."
mkdir -p /data/${MINIO_BUCKET}

echo "Preparing policy..."
cat > public-get.json << EOF
{
  "Version": "2012-10-17",
  "Statement": [
    {
      "Action": [
        "s3:GetObject"
      ],
      "Effect": "Allow",
      "Principal": {
        "AWS": [
          "*"
        ]
      },
      "Resource": [
        "arn:aws:s3:::${MINIO_BUCKET}/*"
      ],
      "Sid": ""
    }
  ]
}
EOF

echo "Starting minio server..."
minio server /data --console-address ":9001" &
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
  mc alias set s3 http://127.0.0.1:9000 "${MINIO_ROOT_USER}" "${MINIO_ROOT_PASSWORD}"
  let SUCCESS=$?

  echo "Setting bucket policy..."
  mc policy set-json public-get.json "s3/${MINIO_BUCKET}"
  let SUCCESS=$?

  let RETRIES="${RETRIES}+1"
done

echo "Ready."
wait ${PID}
