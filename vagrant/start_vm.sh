#!/bin/bash

mkdir -p txbr
cp -r ../* txbr/

echo ""
echo "---- Starting VM This may take a few minutes -----"
echo ""
vagrant up
echo ""
echo "Connect to vm by invoking: vagrant ssh"
echo ""
echo "or for ssh with X redirected"
echo "vagrant ssh -- -X"
echo ""

