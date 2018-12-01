#!/bin/sh 
touch rpx.log init.ss 
rpx << EOF


../ceres
0
../imp
4

-568 461.1826009347023 0
5

710 0 0
0


init.ss
y
EOF