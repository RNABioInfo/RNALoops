#!/bin/sh

#Shell script for algorithm compilation ease of use by Marius Sebeke
#Written to easier automate the compilation of RNAmotiFold algorithms, including updating the hexdump file used to compile them (with this they would be immobile and couldn't be moved after compilation which is a little annoying)
#VARIABLES: GAPC = location of the installed gapcM, ALG = what to call the output algorithm file (no need for .cc here) ARGS = Arguments (including instance)
BASEDIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd) #Set absolute path from home dir to Misc/Application/RNAmotiFold

#Read arguments from the commandline
for ARGUMENT in "$@"
do
    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=)
    case "$KEY" in
            GAPC)    GAPC=${VALUE} ;;
            ALG)     ALG=${VALUE} ;;
            ARGS)    ARGS=${VALUE} ;;
            FILE)    FILE=${VALUE} ;;
            PERL)    PERL=${VALUE} ;;
            *)   
    esac
done

echo $ARGS

#Check if any of the required arugments are empty (ARGS may be empty, so I'll leave it out here though it should never be)
if [ -z "$ALG" ]; then
    printf '%s\n' "ALG is required" >&2
    exit 2
fi

if [ -z "$FILE" ]; then
    printf '%s\n' "FILE is required" >&2
    exit 2
fi

if [ -z "$PERL" ]; then
    PERL="$(command -v perl)"
    if [ -z "$PERL" ]; then
        echo "No perl interpreter given or found via which"
        exit 1
    else
        echo "No perl interpreter given, using $PERL"
    fi
fi
if [ -z "$GAPC" ]; then
    GAPC="$(command -v gapc)"
    if [ -z "$GAPC" ]; then
        echo "No gapc compiler given or found via which"
        exit 1
    else
        echo "No gapc given, using: $GAPC"
    fi

fi
#Move to RNALoops base dir (3 directories under this files location, change this if you ever move this script!)
cd "$BASEDIR/../../.." || {
    printf '%s\n' "Failed to change to project root" >&2 
    exit 1
}
echo $BASEDIR
#Compile C++ code from gapc file for the given algorithm, using the given arguments and the right algorithm file
"$GAPC"\
    -o "$ALG.cc" \
    -i "$ALG" \
    $ARGS \
    "$FILE" || exit 1


#Add RNA Options because all my algorithms need them
"$PERL" Misc/Applications/addRNAoptions.pl \
    "${ALG}.mf" 0 || exit 1

#Build the algorithm
make -f "${ALG}".mf || exit 1

#Remove all the compilation files generated during the process
rm -fv \
     "${ALG}.mf"\
     "${ALG}.o" \
     "${ALG}.d" \
     "${ALG}.cc" \
     "${ALG}.hh" \
     "${ALG}_main.cc" \
     "${ALG}_main.o" \
     "${ALG}_main.d"