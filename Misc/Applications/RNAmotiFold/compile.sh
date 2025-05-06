#!/bin/sh

#Shell script for algorithm compilation ease of use by Marius Sebeke
#Written to easier automate the compilation of RNAmotiFold algorithms, including updating the hexdump file used to compile them (with this they would be immobile and couldn't be moved after compilation which is a little annoying)
#VARIABLES: GAPC = location of the installed gapcM, ALG = what to call the output algorithm file (no need for .cc here) ARGS = Arguments (including instance)
BASEDIR=$(dirname "$0")

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

if [ ! ${PERL:+1} ]
then
    PERL="$(which perl)"
    if [ "x$PERL" = "x" ]
    then 
    echo "No perl interpreter given or found via which"
    exit 1
    else
    echo "No perl interpreter given, using $PERL"
    fi
fi

if [ ! ${GAPC:+1} ]
then
    GAPC="$(which gapc)"
    if [ "x$GAPC" = "x" ]
    then 
    echo "No gapc compiler given or found via which"
    exit 1
    else
    echo "No gapc given, using: $GAPC"
    fi

fi

cd "$BASEDIR/../../.."
$GAPC -o $ALG.cc -i $ALG $ARGS $FILE
$PERL Misc/Applications/addRNAoptions.pl $ALG.mf 0
make -f "$ALG".mf
rm -fv "$ALG".mf "$ALG".o "$ALG".d "$ALG".cc "$ALG".hh "$ALG"_main.cc "$ALG"_main.o "$ALG"_main.o "$ALG"_main.d
