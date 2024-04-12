#!/bin/bash

# Check if the correct number of arguments are provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 [input_bed_file] [output_bed_file]"
    exit 1
fi

input_bed_file=$1
output_bed_file=$2

# Use awk to process the file
awk -F'\t' '{
    if ($2 > $3) {
        temp = $2;
        $2 = $3;
        $3 = temp;
    }
    print $0;
}' "$input_bed_file" | awk 'BEGIN{OFS="\t"} {print $1,$2,$3,$4}' > "$output_bed_file"
