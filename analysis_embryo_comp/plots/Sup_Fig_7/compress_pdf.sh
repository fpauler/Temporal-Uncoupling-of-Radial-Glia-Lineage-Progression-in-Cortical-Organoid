for f in ./*.pdf; do
    convert -compress jpeg -quality 70 "$f" "compressed_${f##*/}"
done
