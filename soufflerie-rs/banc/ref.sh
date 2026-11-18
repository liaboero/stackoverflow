#!/bin/sh
# Références du banc (version d'origine : contrôles wgpu rétablis).
B=target/release/examples/banc3d.exe
export SOUFFLERIE_CONTROLES=1
$B sphere --sortie banc/ref/sphere.txt --champs banc/ref/sphere.bin
for i in 1 2 3; do
  $B ar4 --sortie banc/ref/ar4_$i.txt --champs banc/ref/ar4_$i.bin
  $B ar4b --sortie banc/ref/ar4b_$i.txt --champs banc/ref/ar4b_$i.bin
done
