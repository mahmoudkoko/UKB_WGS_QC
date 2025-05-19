#!/bin/bash

#####
# ENV
#####

set -euo pipefail
# Uncomment for debugging: set -x

# BCFtools is pre-installed in this applet. This will export the path to the required libraries to run it or its plugins
export LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH
export BCFTOOLS_PLUGINS=/usr/lib

##################
# Configurations
##################

# This serves as a reminder of the input parameters in dxapp.js. Going forward, CAPs will indicate variables not specific to a function.
PROJECT_ID="${project_id}"
VCF_LIST="${project_id}:${ukb23374_vcf_list}"
OUTPUT_DIR="${project_id}:${qc_output_dir}"
N_JOBS="${parallel_qc_jobs}"  # FIXED: Correct variable name


#########################
# Global cleanup function
#########################

# Global cleanup tracking
declare -a GLOBAL_TEMP_DIRS=()
declare -a GLOBAL_TEMP_FILES=()

global_cleanup() {
    local exit_code=$?
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting cleanup..." >&2
    
    # Kill any remaining parallel jobs
    pkill -P $$ 2>/dev/null || true
    
    # Clean up any remaining temporary directories
    for dir in "${GLOBAL_TEMP_DIRS[@]}"; do
        if [[ -d "$dir" ]]; then
            echo "Cleaning up directory: $dir" >&2
            rm -rf "$dir" 2>/dev/null || true
        fi
    done
    
    # Clean up any remaining temporary files
    for file in "${GLOBAL_TEMP_FILES[@]}"; do
        if [[ -f "$file" ]]; then
            echo "Cleaning up file: $file" >&2
            rm -f "$file" 2>/dev/null || true
        fi
    done
    
    # Log final status
    if [[ $exit_code -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup completed successfully" >&2
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup completed after error (exit code: $exit_code)" >&2
    fi
    
    exit $exit_code
}

# Set up global cleanup trap
trap global_cleanup EXIT INT TERM


###########################################
# Function to check for ongoing applet runs
###########################################
#!/bin/bash

# IMPROVED: Use safe shell settings with proper error handling
set -euo pipefail
# Uncomment for debugging: set -x

#####
# ENV
#####

# This will export the required libraries to run bcftools and its plugins
export LD_LIBRARY_PATH=/usr/lib:$LD_LIBRARY_PATH
export BCFTOOLS_PLUGINS=/usr/lib

##################
# Configurations
##################

PROJECT_ID="$project_id"
VCF_LIST="$ukb23374_vcf_list"
OUTPUT_DIR="$qc_output_dir"
N_JOBS="$parallel_qc_jobs"  # FIXED: Correct variable name

# Global cleanup tracking
declare -a GLOBAL_TEMP_DIRS=()
declare -a GLOBAL_TEMP_FILES=()

#############
# Global cleanup function
#############

global_cleanup() {
    local exit_code=$?
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting cleanup..." >&2
    
    # Kill any remaining parallel jobs
    pkill -P $$ 2>/dev/null || true
    
    # Clean up any remaining temporary directories
    for dir in "${GLOBAL_TEMP_DIRS[@]}"; do
        if [[ -d "$dir" ]]; then
            echo "Cleaning up directory: $dir" >&2
            rm -rf "$dir" 2>/dev/null || true
        fi
    done
    
    # Clean up any remaining temporary files
    for file in "${GLOBAL_TEMP_FILES[@]}"; do
        if [[ -f "$file" ]]; then
            echo "Cleaning up file: $file" >&2
            rm -f "$file" 2>/dev/null || true
        fi
    done
    
    # Log final status
    if [[ $exit_code -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup completed successfully" >&2
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] Cleanup completed after error (exit code: $exit_code)" >&2
    fi
    
    exit $exit_code
}

# Set up global cleanup trap
trap global_cleanup EXIT INT TERM

#############
# Utility functions
#############

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Function to safely add directories to cleanup list
track_temp_dir() {
    local dir="$1"
    GLOBAL_TEMP_DIRS+=("$dir")
    log_message "Tracking temp directory: $dir"
}

# Function to safely add files to cleanup list
track_temp_file() {
    local file="$1"
    GLOBAL_TEMP_FILES+=("$file")
    log_message "Tracking temp file: $file"
}

#############
# QC function with integrated cleanup
#############

process_graphtyper_vcf_block() {
    # Local cleanup tracking for this function
    local temp_files_local=()
    local temp_dirs_local=()
    local cleanup_stage=0
    
    # Function-specific cleanup
    cleanup_vcf_processing() {
        local exit_code=$?
        
        log_message "Cleaning up VCF processing (stage: $cleanup_stage, exit: $exit_code)"
        
        # Progressive cleanup based on stage
        case $cleanup_stage in
            0) 
                log_message "No cleanup needed - function initialization"
                ;;
            1) 
                log_message "Cleaning up after download failure..."
                [[ -f "${file_name_prefix}.vcf.gz" ]] && rm -f "${file_name_prefix}.vcf.gz" 2>/dev/null || true
                ;;
            2) 
                log_message "Cleaning up after BCF processing failure..."
                rm -f "${file_name_prefix}".{vcf.gz,bcf,bcf.csi} 2>/dev/null || true
                ;;
            3) 
                log_message "Cleaning up after plink2 processing failure..."
                rm -f "${file_name_prefix}".{vcf.gz,bcf,bcf.csi} 2>/dev/null || true
                rm -f "${file_name_prefix}.qc".{pgen,psam,pvar.zst,log,scount.zst} 2>/dev/null || true
                ;;
            4)
                log_message "Cleaning up working directory after processing failure..."
                ;;
        esac
        
        # Always clean up the working directory at the end
        if [[ -n "${vcf_dir:-}" && -d "$vcf_dir" ]]; then
            log_message "Removing working directory: $vcf_dir"
            rm -rf "$vcf_dir" 2>/dev/null || true
        fi
        
        # Clean up tracked local files
        for file in "${temp_files_local[@]}"; do
            [[ -f "$file" ]] && rm -f "$file" 2>/dev/null || true
        done
        
        # Clean up tracked local directories  
        for dir in "${temp_dirs_local[@]}"; do
            [[ -d "$dir" ]] && rm -rf "$dir" 2>/dev/null || true
        done
        
        return $exit_code
    }
    
    # Set up function-specific cleanup trap
    trap cleanup_vcf_processing RETURN
    
    # Stage 0: Function initialization
    cleanup_stage=0
    
    # Read input and validate
    local query_input="$1"
    local vcf_hash=""
    local vcf_name=""
    local chr=""
    local blk=""
    local vcf_dir=""
    local file_name_prefix=""
    local output_files_array=()  # FIXED: Properly declare array
    
    # Determine if input is file hash or path
    if [[ "$query_input" =~ ^file- ]]; then
        vcf_hash="$query_input"
    else
        # Try to resolve as path
        vcf_hash=$(dx describe "$query_input" 2>/dev/null | grep "^ID" | tr -s ' ' | cut -d' ' -f2) || {
            log_message "ERROR: Cannot resolve $query_input to file ID"
            return 2
        }
    fi
    
    # Get VCF name and validate
    vcf_name=$(dx describe "$vcf_hash" | grep "^Name" | tr -s ' ' | cut -d' ' -f2)
    
    if [[ "${vcf_name:0:10}" != "ukb23374_c" ]] || [[ "${vcf_name: -10}" != "_v1.vcf.gz" ]]; then
        log_message "ERROR: $vcf_name does not match expected pattern (ukb23374_cx_bx_v1.vcf.gz)"
        return 2
    fi
    
    # Check for existing outputs
    if dx find data --name "${vcf_name%.vcf.gz}.qc_outputs.txt" --folder "${OUTPUT_DIR}/logs/" --brief | grep -q . ; then
        log_message "Found existing outputs for $vcf_name. Skipping."
        return 2
    fi
    
    # Extract chromosome and block
    chr=$(echo "$vcf_name" | cut -f2 -d'_' | tr -d 'c')
    blk=$(echo "$vcf_name" | cut -f3 -d'_' | tr -d 'b')
    
    # Create working directory
    vcf_dir="${vcf_hash#file-}"
    mkdir -p "$vcf_dir"
    temp_dirs_local+=("$vcf_dir")  # Track for cleanup
    track_temp_dir "$vcf_dir"      # Track globally too
    
    file_name_prefix="${vcf_dir}/${vcf_name%".vcf.gz"}"
    
    # Stage 1: Download
    cleanup_stage=1
    log_message "Downloading block $blk of chromosome $chr"
    
    if ! dx download -f "${PROJECT_ID}:${vcf_hash}" --output "$vcf_dir/" --no-progress; then
        log_message "ERROR: Failed to download $vcf_hash"
        return 2
    fi
    
    temp_files_local+=("${file_name_prefix}.vcf.gz")
    track_temp_file "${file_name_prefix}.vcf.gz"
    
    # Check if VCF has variants
    if ! bcftools view -H "${file_name_prefix}.vcf.gz" | head -n 1 | grep -q . ; then
        log_message "No variants in this VCF. Skipping."
        return 1
    fi
    
    # Stage 2: BCF processing
    cleanup_stage=2
    log_message "Applying genotype filters..."
    
    # The long bcftools pipeline (with error handling)
    if ! {
        bcftools view --no-version -Ou --threads 2 "${file_name_prefix}.vcf.gz" |\
        bcftools annotate --no-version -Ou -x ID,^INFO/MQ,^INFO/AAScore,^INFO/QD,^INFO/QDalt,^FORMAT/GT,^FORMAT/AD,^FORMAT/DP,^FORMAT/GQ,^FORMAT/PL |\
        bcftools norm --no-version -Ou -m-any --keep-sum AD |\
        bcftools view --no-version -Ou  -m2 -M2 |\
        bcftools +setGT --no-version -Ou -- -n u -t a |\
        bcftools +setGT --no-version -Ou -- -n . -t ./x |\
        bcftools +setGT --no-version -Ou -- -n . -t q -i ' GT="alt" & ( sSUM(FMT/AD) / FMT/DP ) < 0.75 ' |\
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
        bcftools +fill-tags --no-version -Ou -- -t 'INFO/QAD= ( 0 - ( ( SUM(FMT/GL) - SUM(FMT/GL[:0]) ) / ( SUM(FMT/AD) - SUM(FMT/AD[:0]) ) ) )' |\
        bcftools +fill-tags --no-version -Ou -- -t 'INFO/LEN_REF:1=STRLEN(REF),INFO/LEN_ALT:1=STRLEN(ALT),INFO/DP_AVG=MEAN(FMT/DP),INFO/GQ_AVG=MEAN(FMT/GQ),INFO/DP10_REF=( COUNT(GT="0/0" & FMT/AD[:0] < 10 ) / COUNT(GT="0/0") ),INFO/DP10_HOM=( COUNT(GT="1/1" & FMT/AD[:1] < 10 ) / COUNT(GT="1/1") ),INFO/DP10_HET=( COUNT(GT="0/1" & SUM(FMT/AD) < 10 ) / COUNT(GT="0/1") ),INFO/GQ40_REF=( COUNT(GT="0/0" & FMT/GQ < 40 & FMT/plGQ < 40 ) /  COUNT(GT="0/0") ),INFO/GQ40_HOM=( COUNT(GT="1/1" & FMT/GQ < 40 & FMT/plGQ < 40 ) /  COUNT(GT="1/1") ),INFO/GQ40_HET=( COUNT(GT="0/1" & FMT/GQ < 40 & FMT/plGQ < 40 ) /  COUNT(GT="0/1") ),INFO/VAF40_HET=( COUNT(GT="0/1" & (VAF < 40 | VAF > 60) ) / COUNT(GT="0/1") )' |\
        bcftools annotate --no-version -Ou -x FORMAT/AD,FORMAT/DP,FORMAT/GQ,FORMAT/GL,FORMAT/plGQ,FORMAT/VAF,FORMAT/pAB |\
        bcftools view --no-version -Ob -i 'AC>0' -l 2 -o "${file_name_prefix}.bcf" --write-index=csi
    }; then
        log_message "ERROR: BCFtools pipeline failed"
        return 2
    fi
    
    temp_files_local+=("${file_name_prefix}.bcf" "${file_name_prefix}.bcf.csi")
    track_temp_file "${file_name_prefix}.bcf"
    track_temp_file "${file_name_prefix}.bcf.csi"
    
    # Clean up input VCF to save space
    rm "${file_name_prefix}.vcf.gz"
    
    # Check if there are variants after QC
    if ! bcftools view -H "${file_name_prefix}.bcf" | head -n 1 | grep -q . ; then
        log_message "No variants remaining after QC. Skipping."
        return 1
    fi
    
    # Stage 3: Plink2 processing
    cleanup_stage=3
    log_message "Converting BCF to Plink2 files..."
    
    if ! plink2 --make-pgen 'vzs' \
        --bcf "${file_name_prefix}.bcf" \
        --vcf-half-call m \
        --threads 2 \
        --memory 3000 \
        --out "${file_name_prefix}.qc"; then
        log_message "ERROR: plink2 conversion failed"
        return 2
    fi
    
    temp_files_local+=("${file_name_prefix}.qc.pgen" "${file_name_prefix}.qc.psam" "${file_name_prefix}.qc.pvar.zst" "${file_name_prefix}.qc.log")
    
    # Clean up BCF files
    rm "${file_name_prefix}.bcf" "${file_name_prefix}.bcf.csi"
    
    # Stage 4: Create sites-only VCF and collect metrics
    cleanup_stage=4
    log_message "Creating sites-only VCF and collecting metrics..."
    
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
        return 2
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
        return 2
    fi
    
    # Index the VCF
    if ! bcftools index "${file_name_prefix}.qc.sites.vcf.gz"; then
        log_message "ERROR: Failed to index sites VCF"
        return 2
    fi
    
    temp_files_local+=("${file_name_prefix}.qc.sites.vcf.gz" "${file_name_prefix}.qc.sites.vcf.gz.csi" "${file_name_prefix}.qc.bed")
    
    # Create new pvar file
    if ! {
        zcat "${file_name_prefix}.qc.sites.vcf.gz" |\
        awk -F"\t" '/^##fileformat/ || /^##FILTER/ || /^##INFO/{next;}/^##/{print}/^#CHROM/{$7=$8;NF=7;print}!/^#/{$1=gsub("^chr","",$1);$6=".";$7=".";NF=7;print}' OFS="\t" |\
        zstd --no-progress -o "${file_name_prefix}.qc.pvar.zst"
    }; then
        log_message "ERROR: Failed to create new pvar file"
        return 2
    fi
    
    # Collect quality metrics
    log_message "Collecting variant-level QC metrics..."
    if ! {
        bcftools query -f '%ID %TYPE %LEN_REF %LEN_ALT %NS %AN %AC %AC_Het %AC_Hom %AC_Hemi %FILTER %QUAL %MQ %AAScore %QD %QDalt %QAD %ExcHet %DP_AVG %GQ_AVG %DP10_REF %DP10_HOM %DP10_HET %GQ40_REF %GQ40_HOM %GQ40_HET %VAF40_HET' "${file_name_prefix}.qc.sites.vcf.gz" |\
        tr ' ' '\t' |\
        gzip > "${file_name_prefix}.qc.quality_scores.tsv.gz"
    }; then
        log_message "ERROR: Failed to collect quality metrics"
        return 2
    fi
    
    # Collect sample metrics
    log_message "Collecting sample QC metrics..."
    if ! plink2 --pfile "${file_name_prefix}.qc" 'vzs' \
        --sample-counts 'zs' 'cols=homref,homalt,homaltsnp,het,ts,tv,single,missing' \
        --threads 2 \
        --memory 3000 \
        --out "${file_name_prefix}.qc"; then
        log_message "ERROR: Failed to collect sample metrics"
        return 2
    fi
    
    # Reformat sample counts
    if ! {
        zstdcat "${file_name_prefix}.qc.scount.zst" |\
        awk '!/^#/' |\
        gzip > "${file_name_prefix}.qc.sample_stats.tsv.gz"
    }; then
        log_message "ERROR: Failed to reformat sample counts"
        return 2
    fi
    
    temp_files_local+=("${file_name_prefix}.qc.quality_scores.tsv.gz" "${file_name_prefix}.qc.sample_stats.tsv.gz")
    
    # Organize files for upload
    log_message "Organizing files for upload..."
    mkdir -p "${vcf_dir}/pfiles" "${vcf_dir}/stats" "${vcf_dir}/sites" "${vcf_dir}/logs"
    
    mv "${file_name_prefix}.qc.pgen" "${file_name_prefix}.qc.psam" "${file_name_prefix}.qc.pvar.zst" "${vcf_dir}/pfiles/"
    mv "${file_name_prefix}.qc.bed" "${file_name_prefix}.qc.sites.vcf.gz" "${file_name_prefix}.qc.sites.vcf.gz.csi" "${vcf_dir}/sites/"
    mv "${file_name_prefix}.qc.quality_scores.tsv.gz" "${file_name_prefix}.qc.sample_stats.tsv.gz" "${vcf_dir}/stats/"
    
    # Clean up intermediate files
    rm -f "${file_name_prefix}.qc.scount.zst" "${file_name_prefix}.qc.pvar.zst.bk" "${file_name_prefix}.qc.log"
    
    # Upload files
    log_message "Uploading outputs to destination directory..."
    local dx_upload_ids
    if ! mapfile -t dx_upload_ids < <(dx upload --recursive "$vcf_dir/" --tag "wgs_qc_block" --brief); then
        log_message "ERROR: Failed to upload files"
        return 2
    fi
    
    # Add to output array
    output_files_array+=("${chr}")
    output_files_array+=("${blk}")
    output_files_array+=("${vcf_hash}")
    output_files_array+=("${dx_upload_ids[@]}")
    
    # Upload the final output list
    local output_files_list_hash
    if ! output_files_list_hash=$(echo "${output_files_array[@]}" |\
        dx upload - \
            --path "${OUTPUT_DIR}/logs/${vcf_name%.vcf.gz}.qc_outputs.txt" \
            --tag "wgs_qc_block" \
            --property "$(printf 'chromosome=%s' "$chr")" \
            --property "$(printf 'block=%s' "$blk")" \
            --brief); then
        log_message "ERROR: Failed to upload output list"
        return 2
    fi
    
    # Success! Mark cleanup stage as complete so directory won't be deleted
    cleanup_stage=5
    
    log_message "Successfully completed processing for $vcf_name"
    echo "$vcf_name,$output_files_list_hash"
    return 0
}

###################################
# Function to log execution context
###################################

log_execution_context() {
    log_message "=== Execution Context ==="
    log_message "Job ID: ${DX_JOB_ID:-'Not available'}"
    log_message "Applet ID: ${DX_APPLET_ID:-'Not available'}"
    log_message "Executable ID: ${DX_EXECUTABLE_ID:-'Not available'}"
    log_message "Project: ${DX_PROJECT_CONTEXT_ID:-'Not available'}"
    log_message "User: ${DX_USER_ID:-'Not available'}"
    log_message "========================="
}

###########################################
# Function to check for ongoing applet runs
###########################################

check_ongoing_applet_runs() {
    local input_vcf_list="$1"
    
    # Get current applet and job IDs from environment
    local current_applet_id="${DX_APPLET_ID:-}"
    local current_job_id="${DX_JOB_ID:-}"
    
    if [[ -z "$current_applet_id" ]]; then
        log_message "WARNING: Could not determine current applet ID from environment"
        # Fallback to hardcoded name
        current_applet_id="ukb23374_wgs_qc_applet"
    else
        log_message "Using applet ID from environment: $current_applet_id"
    fi
    
    log_message "Checking for ongoing applet runs with same input ..."
    
    # Find running or runnable jobs for this applet (using actual applet ID)
    local ongoing_jobs
    ongoing_jobs=$(dx find jobs \
        --state running,runnable,waiting \
        --tag "wgs_qc_batch" \
        --executable "$current_applet_id" \
        --brief 2>/dev/null || true)
    
    if [[ -z "$ongoing_jobs" ]]; then
        log_message "No ongoing applet runs found"
        return 0
    fi
    
    local job_count=$(echo "$ongoing_jobs" | wc -l)
    log_message "Found $job_count ongoing job(s) for this applet"
    
    # Check each ongoing job for matching input
    local job_id
    while read -r job_id; do
        if [[ -n "$job_id" ]]; then
            # Skip checking our own job
            if [[ "$job_id" == "$current_job_id" ]]; then
                log_message "Skipping self (job: $job_id)"
                continue
            fi
            
            # Get job details without jq - parse the describe output directly
            local job_describe
            job_describe=$(dx describe "$job_id" 2>/dev/null || echo "")
            
            if [[ -z "$job_describe" ]]; then
                log_message "WARNING: Could not describe job $job_id"
                continue
            fi
            
            # Extract input file using grep and awk (works for both JSON and text output)
            local job_input=""
            if echo "$job_describe" | grep -q "ukb23374_vcf_list"; then
                # Try to extract from text output format first
                job_input=$(echo "$job_describe" | grep "ukb23374_vcf_list" | head -1 | awk '{print $NF}' | tr -d '",' || echo "")
                
                # If that didn't work, try JSON format parsing
                if [[ -z "$job_input" || "$job_input" == "ukb23374_vcf_list" ]]; then
                    job_input=$(echo "$job_describe" | grep -A1 "ukb23374_vcf_list" | tail -1 | grep "file-" | awk '{print $1}' | tr -d '",' || echo "")
                fi
            fi
            
            # Get job status using simple grep
            local job_status
            job_status=$(echo "$job_describe" | grep -i "state\|status" | head -1 | awk '{print $NF}' | tr -d '",' || echo "unknown")
            
            # Get creation time
            local job_created
            job_created=$(echo "$job_describe" | grep -i "created" | head -1 | cut -d'"' -f4 || echo "unknown")
            
            if [[ "$job_input" == "$input_vcf_list" ]]; then
                log_message "ERROR: Found ongoing applet run with same input:"
                log_message "  Job ID: $job_id"
                log_message "  Input: $job_input"
                log_message "  Status: $job_status"
                log_message "  Created: $job_created"
                log_message "Aborting to avoid duplicate processing"
                return 1
            else
                log_message "Found job $job_id with different input: $job_input (status: $job_status)"
            fi
        fi
    done <<< "$ongoing_jobs"
    
    log_message "No ongoing runs found with same input - proceeding"
    return 0
}


###################
# Utility functions
###################

log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}

# Function to safely add directories to cleanup list
track_temp_dir() {
    local dir="$1"
    GLOBAL_TEMP_DIRS+=("$dir")
    log_message "Tracking temp directory: $dir"
}

# Function to safely add files to cleanup list
track_temp_file() {
    local file="$1"
    GLOBAL_TEMP_FILES+=("$file")
    log_message "Tracking temp file: $file"
}

#####################################
# QC function with integrated cleanup
#####################################

process_graphtyper_vcf_block() {
    # Local cleanup tracking for this function
    local temp_files_local=()
    local temp_dirs_local=()
    local cleanup_stage=0
    
    # Function-specific cleanup
    cleanup_vcf_processing() {
        local exit_code=$?
        
        log_message "Cleaning up VCF processing (stage: $cleanup_stage, exit: $exit_code)"
        
        # Progressive cleanup based on stage
        case $cleanup_stage in
            0) 
                log_message "No cleanup needed - function initialization"
                ;;
            1) 
                log_message "Cleaning up after download failure..."
                [[ -f "${file_name_prefix}.vcf.gz" ]] && rm -f "${file_name_prefix}.vcf.gz" 2>/dev/null || true
                ;;
            2) 
                log_message "Cleaning up after BCF processing failure..."
                rm -f "${file_name_prefix}".{vcf.gz,bcf,bcf.csi} 2>/dev/null || true
                ;;
            3) 
                log_message "Cleaning up after plink2 processing failure..."
                rm -f "${file_name_prefix}".{vcf.gz,bcf,bcf.csi} 2>/dev/null || true
                rm -f "${file_name_prefix}.qc".{pgen,psam.gz,pvar.zst,log,scount.zst} 2>/dev/null || true
                ;;
            4)
                log_message "Cleaning up working directory after processing failure..."
                ;;
        esac
        
        # Always clean up the working directory at the end
        if [[ -n "${vcf_dir:-}" && -d "$vcf_dir" ]]; then
            log_message "Removing working directory: $vcf_dir"
            rm -rf "$vcf_dir" 2>/dev/null || true
        fi
        
        # Clean up tracked local files
        for file in "${temp_files_local[@]}"; do
            [[ -f "$file" ]] && rm -f "$file" 2>/dev/null || true
        done
        
        # Clean up tracked local directories
        for dir in "${temp_dirs_local[@]}"; do
            [[ -d "$dir" ]] && rm -rf "$dir" 2>/dev/null || true
        done
        
        return $exit_code
    }
    
    # Set up function-specific cleanup trap
    trap cleanup_vcf_processing RETURN
    
    # Stage 0: Function initialization
    cleanup_stage=0
    
    # Read input and validate
    local query_input="$1"
    local vcf_hash=""
    local vcf_name=""
    local chr=""
    local blk=""
    local vcf_dir=""
    local file_name_prefix=""
    local output_files_array=()
    

    # Determine if input is file hash or path
    if [[ "$query_input" =~ ^file- ]]; then
        vcf_hash="$query_input"
    else
        # Try to resolve as path
        vcf_hash=$(dx describe "$query_input" 2>/dev/null | grep "^ID" | tr -s ' ' | cut -d' ' -f2) || {
            log_message "ERROR: Cannot resolve $query_input to file ID"
            return 2
            # Exit code 2 reserved for errors.
        }
    fi
    
    # Get VCF name and validate
    vcf_name=$(dx describe "$vcf_hash" | grep "^Name" | tr -s ' ' | cut -d' ' -f2)
    
    if [[ "${vcf_name:0:10}" != "ukb23374_c" ]] || [[ "${vcf_name: -10}" != "_v1.vcf.gz" ]]; then
        log_message "ERROR: $vcf_name does not match expected pattern (ukb23374_cx_bx_v1.vcf.gz)"
        return 2
    fi
    
    # Check for existing outputs (this prevents the same file from running twice)
    if dx find data --name "${vcf_name%.vcf.gz}.qc_outputs.txt" --folder "${OUTPUT_DIR}/logs/" --brief | grep -q . ; then
        log_message "Found existing outputs for $vcf_name. Skipping."
        return 1
    fi
    
    # Extract chromosome and block
    chr=$(echo "$vcf_name" | cut -f2 -d'_' | tr -d 'c')
    blk=$(echo "$vcf_name" | cut -f3 -d'_' | tr -d 'b')
    
    # Create working directory
    vcf_dir="${vcf_hash#file-}"
    mkdir -p "$vcf_dir"
    temp_dirs_local+=("$vcf_dir")  # Track for cleanup
    track_temp_dir "$vcf_dir"      # Track globally too

    # Create output prefix to write directly in this folder
    file_name_prefix="${vcf_dir}/${vcf_name%".vcf.gz"}"
    
    # Stage 1: Download VCF
    cleanup_stage=1
    log_message "Downloading block $blk of chromosome $chr ..."
    
    if ! dx download -f "${PROJECT_ID}:${vcf_hash}" --output "$vcf_dir/" --no-progress; then
        log_message "ERROR: Failed to download $vcf_hash"
        return 2
    fi
    
    # Track temp files locally and globally
    temp_files_local+=("${file_name_prefix}.vcf.gz")
    track_temp_file "${file_name_prefix}.vcf.gz"
    
    # Check if VCF has variants
    if ! bcftools view -H "${file_name_prefix}.vcf.gz" | head -n 1 | grep -q . ; then
        log_message "No variants in this VCF. Skipping."
        return 1
    fi
    
    # Stage 2: BCF processing
    cleanup_stage=2
    log_message "Applying genotype filters ..."
    
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
        bcftools view --no-version -Ob -i 'AC>0' -l 2 -o "${file_name_prefix}.bcf" --write-index=csi
    }; then
        log_message "ERROR: BCFtools pipeline failed"
        return 2
    fi
    
    # Track temp files
    temp_files_local+=("${file_name_prefix}.bcf" "${file_name_prefix}.bcf.csi")
    track_temp_file "${file_name_prefix}.bcf"
    track_temp_file "${file_name_prefix}.bcf.csi"
    
    # Clean up input VCF to save space
    rm "${file_name_prefix}.vcf.gz"
    
    # Check if there are variants after QC
    if ! bcftools view -H "${file_name_prefix}.bcf" | head -n 1 | grep -q . ; then
        log_message "No variants remaining after QC. Skipping."
        return 1
    fi
    
    # Stage 3: Plink2 processing
    cleanup_stage=3
    log_message "Converting BCF to Plink2 files ..."
    
    if ! plink2 --make-pgen 'vzs' \
        --bcf "${file_name_prefix}.bcf" \
        --vcf-half-call m \
        --threads 2 \
        --memory 3000 \
        --out "${file_name_prefix}.qc"; then
        log_message "ERROR: plink2 conversion failed"
        return 2
    fi
    
    # Track temp files
    temp_files_local+=("${file_name_prefix}.qc.pgen" "${file_name_prefix}.qc.psam" "${file_name_prefix}.qc.pvar.zst" "${file_name_prefix}.qc.log")
    
    # Clean up BCF files
    rm "${file_name_prefix}.bcf" "${file_name_prefix}.bcf.csi"
    
    # Stage 4: Create sites-only VCF and collect metrics
    cleanup_stage=4
    log_message "Creating sites-only VCF ..."
    
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
        return 2
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
        return 2
    fi
    
    # Index the VCF
    if ! bcftools index "${file_name_prefix}.qc.sites.vcf.gz"; then
        log_message "ERROR: Failed to index sites VCF"
        return 2
    fi
    
    temp_files_local+=("${file_name_prefix}.qc.sites.vcf.gz" "${file_name_prefix}.qc.sites.vcf.gz.csi" "${file_name_prefix}.qc.bed")
    
    # Create new pvar file with custom IDs and no INFO fields
    if ! {
        zcat "${file_name_prefix}.qc.sites.vcf.gz" |\
        awk -F"\t" '/^##fileformat/ || /^##FILTER/ || /^##INFO/{next;}/^##/{print}/^#CHROM/{$7=$8;NF=7;print}!/^#/{$1=gsub("^chr","",$1);$6=".";$7=".";NF=7;print}' OFS="\t" |\
        zstd --no-progress -o "${file_name_prefix}.qc.pvar.zst"
    }; then
        log_message "ERROR: Failed to create new pvar file"
        return 2
    fi
    
    # Collect quality metrics
    log_message "Collecting variant-level QC metrics ..."
    if ! {
        bcftools query -f '%ID %TYPE %LEN_REF %LEN_ALT %NS %AN %AC %AC_Het %AC_Hom %AC_Hemi %FILTER %QUAL %MQ %AAScore %QD %QDalt %QAD %ExcHet %DP_AVG %GQ_AVG %DP10_REF %DP10_HOM %DP10_HET %GQ40_REF %GQ40_HOM %GQ40_HET %VAF40_HET' "${file_name_prefix}.qc.sites.vcf.gz" |\
        tr ' ' '\t' |\
        gzip > "${file_name_prefix}.qc.quality_scores.tsv.gz"
    }; then
        log_message "ERROR: Failed to collect quality metrics"
        return 2
    fi
    
    # Collect sample metrics
    log_message "Collecting sample QC metrics ..."
    if ! plink2 --pfile "${file_name_prefix}.qc" 'vzs' \
        --sample-counts 'zs' 'cols=homref,homalt,homaltsnp,het,ts,tv,single,missing' \
        --threads 2 \
        --memory 3000 \
        --out "${file_name_prefix}.qc"; then
        log_message "ERROR: Failed to collect sample metrics"
        return 2
    fi
    
    # Reformat sample counts
    if ! {
        zstdcat "${file_name_prefix}.qc.scount.zst" |\
        awk '!/^#/' |\
        gzip > "${file_name_prefix}.qc.sample_stats.tsv.gz"
    }; then
        log_message "ERROR: Failed to reformat sample counts"
        return 2
    fi
    
    # Compress the psam file
    gzip -c "${file_name_prefix}.qc.psam" > "${file_name_prefix}.qc.psam.gz"

    # Track temp files
    temp_files_local+=("${file_name_prefix}.qc.psam.gz" "${file_name_prefix}.qc.quality_scores.tsv.gz" "${file_name_prefix}.qc.sample_stats.tsv.gz")
    
    # Organize files for upload
    log_message "Organizing files for upload ..."
    mkdir -p "${vcf_dir}/pfiles" "${vcf_dir}/stats" "${vcf_dir}/sites" "${vcf_dir}/logs"
    
    mv "${file_name_prefix}.qc.pgen" "${file_name_prefix}.qc.psam.gz" "${file_name_prefix}.qc.pvar.zst" "${vcf_dir}/pfiles/"
    mv "${file_name_prefix}.qc.bed" "${file_name_prefix}.qc.sites.vcf.gz" "${file_name_prefix}.qc.sites.vcf.gz.csi" "${vcf_dir}/sites/"
    mv "${file_name_prefix}.qc.quality_scores.tsv.gz" "${file_name_prefix}.qc.sample_stats.tsv.gz" "${vcf_dir}/stats/"
    
    # Clean up intermediate files
    rm -f "${file_name_prefix}.qc.scount.zst" "${file_name_prefix}.qc.pvar.zst.bk" "${file_name_prefix}.qc.log" "${file_name_prefix}.qc.psam"
    
    # Upload files
    log_message "Uploading outputs to destination directory ..."
    local dx_upload_ids
    if ! mapfile -t dx_upload_ids < <(dx upload --recursive "$vcf_dir/" --tag "wgs_qc_output" --path ${OUTPUT_DIR} --brief); then
        log_message "ERROR: Failed to upload files"
        return 2
    fi
    
    # Add to output array
    output_files_array+=("${chr}")
    output_files_array+=("${blk}")
    output_files_array+=("${vcf_hash}")
    output_files_array+=("${dx_upload_ids[@]}")
    
    # Upload the final output list
    local output_files_list_hash
    if ! output_files_list_hash=$(echo "${output_files_array[@]}" |\
        dx upload - \
            --path "${OUTPUT_DIR}/logs/${vcf_name%.vcf.gz}.qc_outputs.txt" \
            --tag "wgs_qc_log" \
            --property "$(printf 'chromosome=%s' "$chr")" \
            --property "$(printf 'block=%s' "$blk")" \
            --brief); then
        log_message "ERROR: Failed to upload output list"
        return 2
    fi
    
    # Success! Mark cleanup stage as complete so directory won't be deleted
    cleanup_stage=5
    
    log_message "Successfully completed processing for $vcf_name"
    echo "$vcf_name,$output_files_list_hash"
    return 0
}


############################
# Main function with cleanup
############################

main() {
    # Execution context

    log_execution_context

    # Log files
    local log_dir="./LOG_DIR"
    local vcf_list_file=""
    local joblog=""
    local jobsummary=""
    
   # Check for ongoing applet runs with same input
    if ! check_ongoing_applet_runs "$VCF_LIST"; then
        log_message "FATAL: Duplicate applet run detected - exiting"
        exit 1
    fi
    

    # Track main temp directory
    track_temp_dir "$log_dir"
    
    # Create log directory
    mkdir -p "$log_dir"
    
    # Download the input file
    vcf_list_file=$(dx describe "$VCF_LIST" | grep "^Name" | tr -s ' ' | cut -d' ' -f2)
    

    if [[ "${vcf_list_file}" =~ ^vcf_list_ ]] || [[ "${vcf_list_file: -10}" != ".input.txt" ]]; then
        log_message "ERROR: $vcf_list_file does not match expected pattern (vcf_list_n_of_N.input.txt)"
        return 1
    elif ! dx download "${VCF_LIST}" -o "${log_dir}/${vcf_list_file}"; then
        log_message "ERROR: Failed to download VCF list file"
        exit 1
    fi
    
    joblog="${log_dir}/${vcf_list_file%.input.txt}.log"
    jobsummary="${log_dir}/${vcf_list_file%.input.txt}.summary.txt.gz"
    
    track_temp_file "$joblog"
    track_temp_file "$jobsummary"
    
    # Export the function for parallel
    export -f process_graphtyper_vcf_block
    export -f log_message
    export -f track_temp_dir
    export -f track_temp_file
    
    # Run parallel jobs with timeout and retry
    log_message "Starting parallel processing with $N_JOBS jobs ..."
    
    parallel \
        --jobs "$N_JOBS" \
        --results "$log_dir" \
        --joblog "$joblog" \
        --timeout 7200 \
        --retries 1 \
        process_graphtyper_vcf_block :::: "${log_dir}/${vcf_list_file}"
    
    # Generate summary
    log_message "Generating summary ..."
    echo "----- Summary -----"
    echo "Input: $vcf_list_file"
    echo "RECORDS: $(awk 'NR>1' "$joblog" | wc -l)"
    echo "COMPLETED: $(awk -F"\t" 'NR > 1 && $7 == 0' "$joblog" | wc -l)"
    echo "SKIPPED: $(awk -F"\t" 'NR > 1 && $7 == 1' "$joblog" | wc -l)"
    echo "FAILED: $(awk -F"\t" 'NR > 1 && $7 > 1' "$joblog" | wc -l)"
    
    # Combine all logs
    log_message "Combining logs ..."
    if ! {
        find "$log_dir" -name stderr -exec cat {} + | gzip > "${jobsummary}"
        cat "$joblog" | gzip >> "${jobsummary}"
    }; then
        log_message "WARNING: Failed to combine some logs"
    fi
    

    # Upload final summary and capture the file ID for output
    log_message "Uploading final summary ..."
    local qc_log_file_id
    if ! qc_log_file_id=$(dx upload "${jobsummary}" --brief); then
        log_message "WARNING: Failed to upload final summary"
        exit 1
    fi
    
    # Set the output for the applet
    dx-jobutil-add-output qc_log_file "$qc_log_file_id"
    
    log_message "Main function completed successfully"
    log_message "Output log file ID: $qc_log_file_id"

}