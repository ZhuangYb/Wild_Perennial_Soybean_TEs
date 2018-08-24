#!/usr/bin/bash
`cat temp1.txt|perl -a -F'_' -lne 'print join("\t",@F)' |sort|uniq|sort -k1,1 -k6n,6 >temp2.txt`