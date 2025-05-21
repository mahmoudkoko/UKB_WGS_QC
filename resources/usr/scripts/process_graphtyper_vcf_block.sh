#####################################
# QC function with integrated cleanup
#####################################

process_graphtyper_vcf_block() {

    # Local cleanup tracking for this function
    local temp_files_local=()
    local temp_dirs_local=()
    
    # Function-specific cleanup
    cleanup_vcf_processing() {
        local exit_code=$?
        
        log_message "INFO: Cleaning up working directory"        
        
        # Clean up tracked local files
        for file in "${temp_files_local[@]}"; do
            [[ -f "$file" ]] && rm -f "$file" 2>/dev/null || true
        done
        
        # Clean up tracked local directories
        for dir in "${temp_dirs_local[@]}"; do
            [[ -d "$dir" ]] && rm -rf "$dir" 2>/dev/null || true
        done

        # Always clean up the working directory at the end
        if [[ -n "${vcf_dir:-}" && -d "$vcf_dir" ]]; then
            rm -rf "$vcf_dir" 2>/dev/null || true
        fi
        
        return $exit_code
    }
    
    # Set up function-specific cleanup trap
    trap cleanup_vcf_processing RETURN
    
    # Stage 0: Function initialization
    
    # Read input, remove any quotes and initiate other variables.
    local vcf_input="$(echo $1 | sed -e "s/project-.*://" -e 's/^"//' -e 's/"$//')"
    local vcf_desc=""
    local vcf_hash=""
    local vcf_name=""
    local vcf_dir=""
    local chr=""
    local blk=""
    local file_name_prefix=""
    local output_files_list_hash=""
    local output_files_array=()
    
    log_message "INFO: Processing Input: $1"        

    # Determine if input is file hash or path
    vcf_desc=$(dx describe --json "${DX_PROJECT_CONTEXT_ID}:$vcf_input" 2>/dev/null )
    vcf_hash="$(echo $vcf_desc | jq -r .id )"
    vcf_name="$(echo $vcf_desc | jq -r .name )"

    # Extract chromosome and block
    chr=$(echo "$vcf_name" | cut -f2 -d'_' | tr -d 'c')
    blk=$(echo "$vcf_name" | cut -f3 -d'_' | tr -d 'b')
    
    # Add to output array
    output_files_array+=("${chr}")
    output_files_array+=("${blk}")
    output_files_array+=("${vcf_hash}")


    if [[ -z "$vcf_name" ]] || [[ -z "$chr" ]] || [[ -z "$blk" ]]; then
            log_message "ERROR: Cannot resolve vcf file ID or path"
            return 1
    
    # Validate VCF name
    elif [[ ! "$vcf_name" =~ ^ukb23374_c ]] || [[ ! "$vcf_name" =~ "_v1.vcf.gz$" ]] || [[ ! "$chr" =~ ^([1-9]|1[0-9]|2[0-2]|X|x)$ ]] || [[ ! "$blk" =~ ^-?[0-9]+$ ]]; then
        log_message "ERROR: $vcf_name does not match expected pattern (ukb23374_cx_bx_v1.vcf.gz)"
        return 1
  
    # Check for existing outputs (this prevents the same file from running twice by mistake)
    elif dx find data --name "${vcf_name%.vcf.gz}.qc_outputs.txt" --folder "${DX_PROJECT_CONTEXT_ID}:${OUTPUT_DIR}/logs/" --brief | grep -q . ; then
        log_message "INFO: Found existing outputs for $vcf_name "
        if ! output_files_list_hash=$(dx describe --json "${DX_PROJECT_CONTEXT_ID}:${OUTPUT_DIR}/logs/${vcf_name%.vcf.gz}.qc_outputs.txt" 2> /dev/null); then
            log_message "ERROR: could not pull the hash of the previous output list"
            return 1
        else
            log_message "INFO: Skipping $vcf_name and reporting previous output"
            echo "$vcf_name,$output_files_list_hash"
            return 0
        fi

    else
        # Continue if no errors
        log_message "INFO: Processing vcf file: $vcf_name"
    fi
    
    # Create working directory
    vcf_dir="${vcf_hash#file-}"
    mkdir -p "$vcf_dir"
    temp_dirs_local+=("$vcf_dir")  # Track for cleanup

    track_temp_dir "$vcf_dir"      # Track globally too

    # Create output prefix to write directly in this folder
    file_name_prefix="${vcf_dir}/${vcf_name%".vcf.gz"}"
    
    # Stage 1: Download VCF

    log_message "INFO: Downloading block $blk of chromosome $chr ..."
    
    if ! dx download -f "${DX_PROJECT_CONTEXT_ID}:${vcf_hash}" --output "$vcf_dir/" --no-progress; then
        log_message "ERROR: Failed to download $vcf_hash"
        return 1
    fi
    
    # Track temp files locally and globally
    temp_files_local+=("${file_name_prefix}.vcf.gz")

    track_temp_file "${file_name_prefix}.vcf.gz"
    
    # Check if VCF has variants
    if bcftools view -H "${file_name_prefix}.vcf.gz" | head -n 1 | grep -q . ; then
 
        log_message "INFO: Applying genotype filters ..."

    else

        log_message "INFO: No variants in this VCF ..."

        output_files_array+=("NA" "NA" "NA" "NA" "NA" "NA" "NA" "NA")
        
        # Upload the final output list
        if ! output_files_list_hash=$(echo "${output_files_array[@]}" |\
            dx upload - \
                --path "${OUTPUT_DIR}/logs/${vcf_name%.vcf.gz}.qc_outputs.txt" \
                --tag "wgs_qc_log" \
                --property "$(printf 'chromosome=%s' "$chr")" \
                --property "$(printf 'block=%s' "$blk")" \
                --brief); then
            log_message "ERROR: Failed to upload final output list"
            return 1
        fi
    
        log_message "INFO: Successfully completed processing for $vcf_name"
        echo "$vcf_name,$output_files_list_hash"
        return 0

    fi
    


    # Stage 2: BCF processing

    
    # BCFtools QC pipeline
    if ! {
        bcftools view --no-version -Ou --threads 2 "${file_name_prefix}.vcf.gz" |\
        bcftools annotate --no-version -Ou -x ID,^INFO/MQ,^INFO/AAScore,^INFO/QD,^INFO/QDalt,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL |\
        bcftools norm --no-version -Ou -m-any --keep-sum AD |\
        bcftools view --no-version -Ou  -m2 -M2 |\
        bcftools +setGT --no-version -Ou -- -n u -t a |\
        bcftools +setGT --no-version -Ou -- -n . -t ./x |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' (GT="0/1" | GT="1/1") & ( sSUM(FMT/AD) / FMT/DP ) < 0.75 ' |\
        bcftools +fill-tags --no-version -Ou -- -t 'VAF' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="0/0" & (FMT/AD[:0] < 5 | FMT/VAF > 0.3 )' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="1/1" & (FMT/AD[:1] < 5 | FMT/VAF < 0.7 )' |\
        bcftools +fill-tags --no-version -Ou -- -t 'FMT/pAB:1=binom(FMT/AD)' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="0/1" & (FMT/AD[:0] < 3 | FMT/AD[:1] < 3 | FMT/VAF < 0.1 | FMT/VAF > 0.9 | FMT/pAB < 0.00001)' |\
        bcftools +fill-tags --no-version -Ou -- -t 'FMT/plGQ:1=int(sMEDIAN(FMT/PL)-sMIN(FMT/PL))' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GQ < 20 & plGQ < 20 ' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="0/0" & ( (PL[:0] >= PL[:1]) | (PL[:0] >= PL[:2]) ) ' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="0/1" & ( (PL[:1] >= PL[:0]) | (PL[:1] >= PL[:2]) ) ' |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="1/1" & ( (PL[:2] >= PL[:0]) | (PL[:2] >= PL[:1]) ) ' |\
        bcftools +fill-tags --no-version -Ou -- -t 'NS,AC,AN,AC_Het,AC_Hom,AC_Hemi,ExcHet' |\
        bcftools +tag2tag --no-version -Ou -- --PL-to-GL --replace |\
        bcftools +fill-tags --no-version -Ou -- -t 'INFO/LEN_REF:1=STRLEN(REF),INFO/LEN_ALT:1=STRLEN(ALT),INFO/QAD= ( 0 - ( ( SUM(FMT/GL) - SUM(FMT/GL[:0]) ) / ( SUM(FMT/AD) - SUM(FMT/AD[:0]) ) ) ),INFO/DP_AVG=MEAN(FMT/DP),INFO/GQ_AVG=MEAN(FMT/GQ),INFO/DP10_REF=( COUNT(GT="0/0" & FMT/AD[:0] < 10 ) / COUNT(GT="0/0") ),INFO/DP10_HOM=( COUNT(GT="1/1" & FMT/AD[:1] < 10 ) / COUNT(GT="1/1") ),INFO/DP10_HET=( COUNT(GT="0/1" & SUM(FMT/AD) < 10 ) / COUNT(GT="0/1") ),INFO/GQ40_REF=( COUNT(GT="0/0" & FMT/GQ < 40 & FMT/plGQ < 40 ) /  COUNT(GT="0/0") ),INFO/GQ40_HOM=( COUNT(GT="1/1" & FMT/GQ < 40 & FMT/plGQ < 40 ) /  COUNT(GT="1/1") ),INFO/GQ40_HET=( COUNT(GT="0/1" & FMT/GQ < 40 & FMT/plGQ < 40 ) /  COUNT(GT="0/1") ),INFO/VAF40_HET=( COUNT(GT="0/1" & (VAF < 40 | VAF > 60) ) / COUNT(GT="0/1") )' |\
        bcftools annotate --no-version -Ou -x FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/GL,FORMAT/plGQ,FORMAT/VAF,FORMAT/pAB |\
        bcftools filter --no-version -Ob -s 'AC0' -m+ -i 'AC=0' -l 2 -o "${file_name_prefix}.bcf" --write-index=csi
    }; then
        log_message "ERROR: BCFtools pipeline failed"
        return 1
    fi
    
    # Track temp files
    temp_files_local+=("${file_name_prefix}.bcf" "${file_name_prefix}.bcf.csi")

    track_temp_file "${file_name_prefix}.bcf"
    track_temp_file "${file_name_prefix}.bcf.csi"
    
    # Clean up input VCF to save space
    rm "${file_name_prefix}.vcf.gz"
    
    # Check if there are variants after QC
    if ! bcftools view -H "${file_name_prefix}.bcf" | head -n 1 | grep -q . ; then
        log_message "ERROR: No variants remaining after QC although original file had variants"
        return 1
    fi
    
    # Stage 3: Plink2 processing

    log_message "INFO: Converting BCF to Plink2 files ..."
    
    if ! plink2 --make-pgen 'vzs' \
        --bcf "${file_name_prefix}.bcf" \
        --vcf-half-call m \
        --threads 2 \
        --memory 3000 \
        --out "${file_name_prefix}.qc"; then
        log_message "ERROR: plink2 conversion failed"
        return 1
    fi
    
    # Track temp files
    temp_files_local+=("${file_name_prefix}.qc.pgen" "${file_name_prefix}.qc.psam" "${file_name_prefix}.qc.pvar.zst" "${file_name_prefix}.qc.log")
    
    # Clean up BCF files
    rm "${file_name_prefix}.bcf" "${file_name_prefix}.bcf.csi"
    
    # Stage 4: Create sites-only VCF and collect metrics

    log_message "INFO: Creating sites-only VCF ..."
    
    # Create sites-only VCF (with error checking)
    mv "${file_name_prefix}.qc.pvar.zst" "${file_name_prefix}.qc.pvar.zst.bk"
    temp_files_local+=("${file_name_prefix}.qc.pvar.zst.bk")
    
    # Create VCF header
    if ! {
        zstdcat "${file_name_prefix}.qc.pvar.zst.bk" |\
        awk -F"\t" 'BEGIN{print "##fileformat=VCFv4.2"}!/^#/{exit}1' OFS="\t" |\
        bgzip > "${file_name_prefix}.qc.sites.vcf.gz"
    }; then
        log_message "ERROR: Failed to create VCF header"
        return 1
    fi
    
    # Append variants
    if ! {
        zstdcat "${file_name_prefix}.qc.pvar.zst.bk" |\
        awk '!/^#/' |\
        awk -F'\t' -v blk="${blk}" -v bed_file="${file_name_prefix}.qc.bed" \
        'NR==1{ start_position=$2 }{$3=$1":"blk":"NR; print "chr"$0}END{ end_position=$2; print "chr"$1,start_position,end_position,"b"blk,NR > bed_file }' OFS="\t" |\
        bgzip >> "${file_name_prefix}.qc.sites.vcf.gz"
    }; then
        log_message "ERROR: Failed to append variants to VCF"
        return 1
    fi
    
    # Index the VCF
    if ! bcftools index "${file_name_prefix}.qc.sites.vcf.gz"; then
        log_message "ERROR: Failed to index sites VCF"
        return 1
    fi
    
    temp_files_local+=("${file_name_prefix}.qc.sites.vcf.gz" "${file_name_prefix}.qc.sites.vcf.gz.csi" "${file_name_prefix}.qc.bed")
    
    # Create new pvar file with custom IDs and no INFO fields
    if ! {
        zcat "${file_name_prefix}.qc.sites.vcf.gz" |\
        awk -F"\t" '/^##fileformat/ || /^##FILTER/ || /^##INFO/{next;}/^##/{print}/^#CHROM/{$7=$8;NF=7;print}!/^#/{$1=gsub("^chr","",$1);$6=".";$7=".";NF=7;print}' OFS="\t" |\
        zstd --no-progress -o "${file_name_prefix}.qc.pvar.zst"
    }; then
        log_message "ERROR: Failed to create new pvar file"
        return 1
    fi
    
    # Collect quality metrics
    log_message "INFO: Collecting variant-level QC metrics ..."
    if ! {
        bcftools query -f '%ID %TYPE %LEN_REF %LEN_ALT %NS %AN %AC %AC_Het %AC_Hom %AC_Hemi %FILTER %QUAL %MQ %AAScore %QD %QDalt %QAD %ExcHet %DP_AVG %GQ_AVG %DP10_REF %DP10_HOM %DP10_HET %GQ40_REF %GQ40_HOM %GQ40_HET %VAF40_HET' "${file_name_prefix}.qc.sites.vcf.gz" |\
        tr ' ' '\t' |\
        gzip > "${file_name_prefix}.qc.quality_scores.tsv.gz"
    }; then
        log_message "ERROR: Failed to collect quality metrics"
        return 1
    fi
    
    # Collect sample metrics
    log_message "INFO: Collecting sample QC metrics ..."
    if ! plink2 --pfile "${file_name_prefix}.qc" 'vzs' \
        --sample-counts 'zs' 'cols=homref,homalt,homaltsnp,het,ts,tv,single,missing' \
        --threads 2 \
        --memory 3000 \
        --out "${file_name_prefix}.qc"; then
        log_message "ERROR: Failed to collect sample metrics"
        return 1
    fi
    
    # Reformat sample counts
    if ! {
        zstdcat "${file_name_prefix}.qc.scount.zst" |\
        awk '!/^#/' |\
        gzip > "${file_name_prefix}.qc.sample_stats.tsv.gz"
    }; then
        log_message "ERROR: Failed to reformat sample counts"
        return 1
    fi
    

    # Compress the psam file
    gzip -c "${file_name_prefix}.qc.psam" > "${file_name_prefix}.qc.psam.gz"

    # Track temp files
    temp_files_local+=("${file_name_prefix}.qc.psam.gz" "${file_name_prefix}.qc.quality_scores.tsv.gz" "${file_name_prefix}.qc.sample_stats.tsv.gz")
    
    # Clean up intermediate files
    rm -f "${file_name_prefix}.qc.scount.zst" "${file_name_prefix}.qc.pvar.zst.bk" "${file_name_prefix}.qc.log" "${file_name_prefix}.qc.psam"

    # Organize files for upload
    log_message "INFO: Organizing files for upload ..."
    mkdir -p "${vcf_dir}/pfiles" "${vcf_dir}/stats" "${vcf_dir}/sites" "${vcf_dir}/logs"
    
    mv "${file_name_prefix}.qc.pgen" "${file_name_prefix}.qc.psam.gz" "${file_name_prefix}.qc.pvar.zst" "${vcf_dir}/pfiles/"
    mv "${file_name_prefix}.qc.bed" "${file_name_prefix}.qc.sites.vcf.gz" "${file_name_prefix}.qc.sites.vcf.gz.csi" "${vcf_dir}/sites/"
    mv "${file_name_prefix}.qc.quality_scores.tsv.gz" "${file_name_prefix}.qc.sample_stats.tsv.gz" "${vcf_dir}/stats/"
    
    
    # Upload files
    log_message "INFO: Uploading outputs to destination directory ..."
    local dx_upload_ids
    if ! mapfile -t dx_upload_ids < <(dx upload --recursive "$vcf_dir/" --tag "wgs_qc_output" --path ${OUTPUT_DIR} --brief); then
        log_message "ERROR: Failed to upload files"
        return 1
    fi
    
    # Add to output array
    output_files_array+=("${dx_upload_ids[@]}")
    
    # Upload the final output list
    if ! output_files_list_hash=$(echo "${output_files_array[@]}" |\
        dx upload - \
            --path "${OUTPUT_DIR}/logs/${vcf_name%.vcf.gz}.qc_outputs.txt" \
            --tag "wgs_qc_log" \
            --property "$(printf 'chromosome=%s' "$chr")" \
            --property "$(printf 'block=%s' "$blk")" \
            --brief); then
        log_message "ERROR: Failed to upload final output list"
        return 1
    fi
    

    log_message "INFO: Successfully completed processing for $vcf_name"
    echo "$vcf_name,$output_files_list_hash"
    return 0
}