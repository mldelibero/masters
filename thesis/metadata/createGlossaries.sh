#! /bin/bash
# Run this file to build the index. Need to integrate into TexShop later.
makeindex -s thesis.ist -t thesis.glg -o thesis.gls thesis.glo
makeindex -s thesis.ist -t thesis.alg -o thesis.acr thesis.acn
makeindex -s thesis.ist -t thesis.slg -o thesis.sym thesis.sbl

