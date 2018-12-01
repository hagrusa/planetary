#!/bin/sh 
touch rpx.log init.ss 
rpx << EOF


../ceres
0
../imp
4

-568 237.93408136001378 0
5

710 0 0
0


init.ss
y
EOF