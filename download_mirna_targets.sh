#Author Anel Ordabayeva

#!/bin/bash

CSV="DEMs_log2fc1_check.csv"
OUTDIR="./MIRNA_1_TARGETS"
mkdir -p "$OUTDIR"

tail -n +2 "$CSV" | cut -d',' -f1 | while read -r mirna; do
    outname=$(echo "$mirna" | tr '/' '_' | tr -d '[:space:]').txt
    outfile="${OUTDIR}/${outname}"
    url="https://rnasysu.com/encori/moduleDownload.php?source=agoClipRNA&type=txt&value=hg38;lncRNA;${mirna};1;0;0;1;None;all"

    echo "Downloading $mirna..."

    wget --timeout=15 --tries=2 --no-check-certificate \
        --header="User-Agent: Mozilla/5.0 (X11; Linux x86_64)" \
        -q -O "$outfile" "$url"

    if [[ ! -s "$outfile" ]]; then
        echo "⚠️  Download failed or empty for: $mirna"
        rm -f "$outfile"
    else
        echo "✅ Downloaded: $mirna"
    fi
done

echo "✅ All downloads attempted."

