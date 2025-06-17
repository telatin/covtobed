#!/bin/bash

# Simple test data generator using samtools
# This script generates synthetic BAM files for testing

set -eou pipefail

# Function to create a simple SAM file
create_sam() {
    local output_file="$1"
    local ref_name="$2" 
    local ref_length="$3"
    local pattern_type="$4"
    
    # SAM header
    cat > "$output_file" << EOF
@HD	VN:1.4	SO:coordinate
@SQ	SN:$ref_name	LN:$ref_length
EOF

    case "$pattern_type" in
        "coverage")
            local start="$5"
            local end="$6"
            local depth="$7"
            
            # Generate reads with specific coverage
            local read_length=100
            for ((i=1; i<=depth; i++)); do
                local pos=$((start + (i-1) * (end - start - read_length) / (depth > 1 ? depth - 1 : 1)))
                echo -e "read_${i}\t0\t$ref_name\t$((pos+1))\t30\t${read_length}M\t*\t0\t0\t$(printf 'A%.0s' $(seq 1 $read_length))\t$(printf '!%.0s' $(seq 1 $read_length))"
            done >> "$output_file"
            ;;
            
        "stranded")
            # Forward strand reads
            for ((i=1; i<=3; i++)); do
                local pos=$((100 + i*10))
                echo -e "forward_${i}\t0\tchr1\t$((pos+1))\t30\t50M\t*\t0\t0\t$(printf 'A%.0s' $(seq 1 50))\t$(printf '!%.0s' $(seq 1 50))"
            done >> "$output_file"
            
            # Reverse strand reads  
            for ((i=1; i<=2; i++)); do
                local pos=$((150 + i*10))
                echo -e "reverse_${i}\t16\tchr1\t$((pos+1))\t30\t50M\t*\t0\t0\t$(printf 'T%.0s' $(seq 1 50))\t$(printf '!%.0s' $(seq 1 50))"
            done >> "$output_file"
            ;;
            
        "single")
            echo -e "single_read\t0\t$ref_name\t1001\t30\t100M\t*\t0\t0\t$(printf 'A%.0s' $(seq 1 100))\t$(printf '!%.0s' $(seq 1 100))" >> "$output_file"
            ;;
    esac
}

# Function to convert SAM to sorted BAM
sam_to_bam() {
    local sam_file="$1"
    local bam_file="$2"
    
    if command -v samtools >/dev/null 2>&1; then
        samtools view -Sb "$sam_file" | samtools sort -o "$bam_file" 2>/dev/null
        # Only try to index if it's a proper BAM file
        if [ -f "$bam_file" ]; then
            samtools index "$bam_file" 2>/dev/null || true
        fi
    else
        echo "Warning: samtools not available, cannot convert to BAM"
        return 1
    fi
}

# Main function
main() {
    if [ $# -lt 2 ]; then
        echo "Usage: $0 <output_file> <pattern_type> [args...]"
        echo "Pattern types:"
        echo "  coverage <start> <end> <depth> - Generate specific coverage pattern"
        echo "  stranded - Generate stranded coverage data"
        echo "  single - Generate single read"
        exit 1
    fi
    
    local output_file="$1"
    local pattern_type="$2"
    shift 2
    
    local sam_file="${output_file%.bam}.sam"
    local ref_name="chr1"
    local ref_length=10000
    
    case "$pattern_type" in
        "coverage")
            if [ $# -lt 3 ]; then
                echo "Coverage pattern requires: start end depth"
                exit 1
            fi
            create_sam "$sam_file" "$ref_name" "$ref_length" "coverage" "$1" "$2" "$3"
            ;;
        "stranded")
            create_sam "$sam_file" "$ref_name" "$ref_length" "stranded"
            ;;
        "single")
            create_sam "$sam_file" "$ref_name" "$ref_length" "single"
            ;;
        *)
            echo "Unknown pattern type: $pattern_type"
            exit 1
            ;;
    esac
    
    if sam_to_bam "$sam_file" "$output_file"; then
        rm -f "$sam_file"
        echo "Generated BAM file: $output_file"
    else
        echo "Generated SAM file: $sam_file (samtools not available for BAM conversion)"
    fi
}

main "$@"