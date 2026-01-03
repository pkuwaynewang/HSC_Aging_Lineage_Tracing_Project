awk -F, '{print NF}' AlleleColonies.txt > UMI_number.txt
paste AlleleAnnotations.txt UMI_number.txt > Barcode_UMI_number.txt
