#!/usr/bin/bash

echo "PostScript beginning"
#root
hadd -f Result_$1.root *.root
#.q
echo "PostScript done"
