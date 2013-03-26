export GPRODMODE=1

gevgen_hadron -p 211 -t 1000060120 -n 500000 -r 100 -k .100 -m hA  
gntpc -f ginuke -i gntp.inuke.100.ghep.root -o Jan_31_13_pip_C_100_v271_hA.ginuke.root 
gevgen_hadron -p 211 -t 1000060120 -n 500000 -r 101 -k .100 -m hN  
gntpc -f ginuke -i gntp.inuke.101.ghep.root -o Jan_31_13_pip_C_100_v271_hN.ginuke.root 

gevgen_hadron -p 211 -t 1000060120 -n 500000 -r 102 -k .160 -m hA  
gntpc -f ginuke -i gntp.inuke.102.ghep.root -o Jan_31_13_pip_C_160_v271_hA.ginuke.root 
gevgen_hadron -p 211 -t 1000060120 -n 500000 -r 103 -k .160 -m hN  
gntpc -f ginuke -i gntp.inuke.103.ghep.root -o Jan_31_13_pip_C_160_v271_hN.ginuke.root 

gevgen_hadron -p 211 -t 1000060120 -n 500000 -r 106 -k .300 -m hA  
gntpc -f ginuke -i gntp.inuke.106.ghep.root -o Jan_31_13_pip_C_300_v271_hA.ginuke.root 
gevgen_hadron -p 211 -t 1000060120 -n 500000 -r 107 -k .300 -m hN  
gntpc -f ginuke -i gntp.inuke.107.ghep.root -o Jan_31_13_pip_C_300_v271_hN.ginuke.root