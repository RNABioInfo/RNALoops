FILEDIR=$(dirname $(realpath "$0"))

HEXDUMP_FILE="$FILEDIR/../../../../Extensions/mot_header.hh"
for ARGUMENT in "$@"
do
    KEY=$(echo $ARGUMENT | cut -f1 -d=)
    VALUE=$(echo $ARGUMENT | cut -f2 -d=)
    case "$KEY" in
            VERSION)    VERSION=${VALUE} ;;
            *)   
    esac
done
rm -fv $HEXDUMP_FILE
MOTIF_FILES="$FILEDIR"/versions/"$VERSION"/*
RFAM_FILES="$FILEDIR"/versions/rfam/*
printf "//RNA 3D Motif Atlas Version: $VERSION\n" >> $HEXDUMP_FILE

for f in $MOTIF_FILES
do
    filename=$(basename "$f")
    extension="${filename##*.}"
    if [ "$extension" = "csv" ]
        then
            grep -v '^#' $f | xxd -i -n $(basename "$f" .csv) | sed 's/\<unsigned\>/static unsigned/g' >> $HEXDUMP_FILE
    fi
done

for f in $RFAM_FILES
do
    filename=$(basename "$f")
    extension="${filename##*.}"
    if [ "$extension" = "csv" ]
        then
            grep -v '^#' $f | xxd -i -n $(basename "$f" .csv) | sed 's/\<unsigned\>/static unsigned/g' >> $HEXDUMP_FILE
    fi
done