#!/bin/sh 
touch rpx.log init.ss 
rpx << EOF

y
../ceres
0
../imp
4

-568 543.8162402007455 0
5

710 0 0
0


init.ss
y
EOF