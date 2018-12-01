#!/bin/sh 
touch rpx.log init.ss 
rpx << EOF


../ceres
0
../imp
4

-568 487.57230233063893 0
5

710 0 0
0


init.ss
y
EOF