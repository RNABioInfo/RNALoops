FILEDIR=$(dirname $(realpath "$0"))

HEXDUMP_FILE="$FILEDIR/../../../../Extensions/mot_header.hh"

rm -fv $HEXDUMP_FILE

echo "#include <unordered_map>" >>$HEXDUMP_FILE
echo "#include <tuple>" >>$HEXDUMP_FILE
echo "#include <string>" >>$HEXDUMP_FILE

RFAM_FILES="$FILEDIR"/versions/rfam/*

declare -a version_array=()

for dir in "$FILEDIR"/versions/combined/*/
do
    FILE=${dir%*/}
    VERSION="${FILE##*/}"
    MOTIF_FILES="$FILEDIR"/versions/combined/"$VERSION"/*

    version_array+=(${VERSION})
    

    printf "//RNA 3D Motif Atlas Version: $VERSION\n" >> $HEXDUMP_FILE
    for f in $MOTIF_FILES
    do
        filename=$(basename "$f")
        varname=$(basename "$f" .csv)
        var_version="${varname}"_"${VERSION}"
        extension="${filename##*.}"
        if [ "$extension" = "csv" ]
            then #chad_gpt generated sed string for replacing "unsigned [filepath]" with "static unsigned [filename]",because xxd on ubuntu 22.04 doesn't have -n parameter
                xxd -i $f | sed -E 's/(unsigned )?(char|int) _[a-zA-Z0-9_]*_([a-zA-Z0-9]+_[a-zA-Z0-9]+)_csv/static \1\2 \3/g' | sed "s/$varname/$var_version/g" >> $HEXDUMP_FILE
        fi
    done
done

for f in $RFAM_FILES
do
    filename=$(basename "$f")
    extension="${filename##*.}"
    if [ "$extension" = "csv" ]
        then
            xxd -i $f | sed -E 's/(unsigned )?(char|int) _[a-zA-Z0-9_]*_([a-zA-Z0-9]+_[a-zA-Z0-9]+)_csv/static \1\2 \3/g' >> $HEXDUMP_FILE
    fi
done

echo "static const std::unordered_map<std::string,std::tuple<unsigned char*, unsigned int, unsigned char*, unsigned int, unsigned char*, unsigned int>> rna3d_versions = {" >> $HEXDUMP_FILE

for version in "${version_array[@]}"
do
    echo '{"'"$version"'",{rna3d_hairpins_'"$version"',rna3d_hairpins_'"$version"'_len,rna3d_internals_'"$version"',rna3d_internals_'"$version"'_len,rna3d_bulges_'"$version"',rna3d_bulges_'"$version"'_len}},' >> $HEXDUMP_FILE
done

echo '};' >> $HEXDUMP_FILE