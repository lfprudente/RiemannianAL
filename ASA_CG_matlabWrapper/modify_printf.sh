#!/usr/bin/env bash

# the idea: replace all "printf" statements with "Printf"
#   Then, if we compile with the -DMEXPRINTF option,
#   "Printf" is linked to "mexPrintf" so that output
#   appears on the Matlab prompt.

# Stephen Becker, March 22 2012; Jan 22 2018

#file='ASA_CG-2.2/asa_cg.c' # change appropriately
file='ASA_CG-3.0/asa_cg.c'


cat > tempFile << EOF
#ifdef MEXPRINTF
#define Printf mexPrintf
#else
#define Printf printf
#endif
EOF

sed 's/printf/Printf/g' $file >> tempFile

mv tempFile $file 
