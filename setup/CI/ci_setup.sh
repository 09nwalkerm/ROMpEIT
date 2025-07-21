#!/bin/bash

echo "Starting jenkins server..."
sudo systemctl start jenkins
echo "Sleeping for 10"
sleep 10
echo "Loading git secret key..."
secret=`head -1 tokens/GIT.SECRET`
echo "Starting ngrok tunnel on port 8081..."
nohup ngrok http 8081 --verify-webhook github --verify-webhook-secret $secret> /dev/null 2>&1 &
echo "Sleeping for 10"
sleep 10
echo "Setting up github webhook..."
python3 webhook.py
echo "Done."
