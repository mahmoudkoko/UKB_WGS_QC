#!/bin/bash

set -euo pipefail

source /usr/scripts/export_global_vars.sh
source /usr/scripts/check_ongoing_applet_runs.sh
source /usr/scripts/process_graphtyper_vcf_block.sh
source /usr/scripts/global_cleanup.sh
source /usr/scripts/utilities.sh


# Export the functions used within the main entry point below for parallel runs
export -f process_graphtyper_vcf_block
export -f log_message
export -f track_temp_dir
export -f track_temp_file


# Set global cleanup trap
trap global_cleanup EXIT INT TERM


# Set global vars
export_global_vars


# Print execution context
log_execution_context


##################
# Main entry point
##################

main() {


    # Initiate dirs and files
    local log_dir="/home/dnanexus/LOG_DIR"
    local log_file="${VCF_LIST_NAME%.input.txt}.log"
    local job_summary="${VCF_LIST_NAME%.input.txt}.summary.txt.gz"
    local job_summary_dx_id=""

    # Track main temp directory
    track_temp_dir "$log_dir"
    track_temp_file "$log_file"
    track_temp_file "$job_summary"    


    # Create log directory

    mkdir -p "$log_dir"
    
    # Process the input file

    log_message "INFO: Getting the VCF list ..."

   # Check for ongoing applet runs with same input
    if ! check_ongoing_applet_runs ; then
        log_message "ERROR: Duplicate applet run detected - exiting"
        exit 1

    elif [[ ! "${VCF_LIST_NAME}" =~ ^vcf_list_ ]] || [[ ! "${VCF_LIST_NAME}" =~ .input.txt$ ]] ; then

        # If the name isn't right, issue an error and exit
        log_message "ERROR: ${VCF_LIST_NAME} does not match expected pattern (vcf_list_n_of_N.input.txt) - exiting"
        exit 1

        # Otherwise download the file
    elif ! dx download "${DX_PROJECT_CONTEXT_ID}:${VCF_LIST_HASH}" -o "${log_dir}/${VCF_LIST_NAME}"; then

        # If the download fails, issue an error and exit
        log_message "ERROR: Failed to download VCF list file - exiting"
        exit 1
    else

        log_message "INFO: List contains $(cat ${log_dir}/${VCF_LIST_NAME} | wc -l) elements"
    fi
    

    # Run parallel jobs with timeout (~ 3 hours) and a single retry
    log_message "INFO: Starting parallel processing with ${N_JOBS} jobs ..."
    
    parallel \
        --jobs "$N_JOBS" \
        --results "${log_dir}" \
        --joblog "${log_dir}/${log_file}" \
        --timeout 12000 \
        --retries 1 \
        process_graphtyper_vcf_block :::: "${log_dir}/${VCF_LIST_NAME}"
    
    log_message "INFO: Generating summary ..."

    # Generate summary (Standard output - will appear on the platform log this way)
    echo "----- Summary -----"
    echo "Input: ${VCF_LIST_NAME} ${VCF_LIST_HASH}"
    echo "RECORDS: $(awk 'NR>1' "${log_dir}/${log_file}" | wc -l)"
    echo "COMPLETED: $(awk -F"\t" 'NR > 1 && $7 == 0' "${log_dir}/${log_file}" | wc -l)"
    echo "SKIPPED: $(awk -F"\t" 'NR > 1 && $7 == 1' "${log_dir}/${log_file}" | wc -l)"
    echo "FAILED: $(awk -F"\t" 'NR > 1 && $7 > 1' "${log_dir}/${log_file}" | wc -l)"
    
    # Combine all logs and upload
    log_message "INFO: Combining logs ..."

    if ! {
        find "${log_dir}" -name stderr -exec cat {} + | gzip > "${log_dir}/${job_summary}"
        cat "${log_dir}/${log_file}" | gzip >> "${log_dir}/${job_summary}"
    }; then
        # If compessing all log files fails, issue an error and exit.
        log_message "ERROR: Failed to collect logs"
        log_message "INFO: QC completed with errors"
        exit 1
    else
        # If compressed successfully, upload and capture the file ID for output
        log_message "INFO: Uploading final summary ..."

        if ! job_summary_dx_id=$(dx upload "${log_dir}/${job_summary}" --path "${DX_PROJECT_CONTEXT_ID}:${VCF_LIST_DIR}/${job_summary}" --brief); then

            # If upload failed, issue an error and exit
            log_message "ERROR: Failed to upload final summary"
            log_message "INFO: QC completed with errors"
            exit 1

        else
            # If uploaded, set the output for the applet

            if ! dx-jobutil-add-output qc_log_file "${job_summary_dx_id}" ; then

                # If setting outputs failed, issue an error and exit
                log_message "ERROR: Failed to capture the final output"
                log_message "INFO: QC completed with errors"
                exit 1

            fi

         
        fi

    fi

    log_message "INFO: QC completed successfully"
    log_message "INFO: Summaries uploaded as ${job_summary_dx_id}"

    # Exit gracefully
    echo "LOG: ${job_summary_dx_id}"
}