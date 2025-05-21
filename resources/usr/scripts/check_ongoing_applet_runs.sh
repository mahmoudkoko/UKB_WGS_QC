###########################################
# Function to check for ongoing applet runs
###########################################

check_ongoing_applet_runs() {
    
    if [[ -z "$DX_APPLET_ID" ]] || [[ "$DX_APPLET_ID" == null ]]; then
        log_message "ERROR: Could not determine current applet ID from environment. Using default: ukb23374_wgs_qc_applet"
        return 1
    else
        log_message "INFO: Using applet ID from environment: $DX_APPLET_ID"
    fi
    
    log_message "INFO: Checking for ongoing applet runs with same input $VCF_LIST_HASH ..."
    
    # Find running or runnable jobs for this applet (using actual applet ID)
    local ongoing_jobs
    ongoing_jobs=$(dx find jobs \
        --project ${DX_PROJECT_CONTEXT_ID} \
        --executable "$DX_APPLET_ID" \
        --tag "wgs_qc_batch" \
        --brief |\
        sed 's/project-.*://' 2>/dev/null || true)
    
    if [[ -z "$ongoing_jobs" ]]; then
        log_message "ERROR: Failed to properly identify running jobs (should at least identify self)"
        return 1
    fi

    local job_count=$(echo "$ongoing_jobs" | wc -l)
    log_message "INFO: Found $job_count ongoing job(s) for this applet"
    
    # Check each ongoing job for matching input
    local job_id
    while read -r job_id; do
        if [[ -n "$job_id" ]]; then

            # Skip checking our own job
            if [[ "$job_id" == "$DX_JOB_ID" ]]; then
                continue
            fi
            
            # Get job details with jq - parse the describe output directly to get the input file
            local job_describe
            job_describe=$(dx describe --json "$job_id"  | sed 's/project-.*://' 2>/dev/null || echo "")
            
            # Skip of no description
            if [[ -z "$job_describe" ]]; then
                log_message "ERROR: Could not describe job $job_id"
                return 1

            else

                # Parse JSON to get the input file ID
                local job_input
                job_input=$(echo "$job_describe" | jq -r '.runInput.ukb23374_vcf_list["$dnanexus_link"]' | sed 's/project-.*://' 2>/dev/null || echo "")

                # If that didn't work or the result isn't as expected, fail
                if [[ -z "$job_input" ]] || [[ ! "$job_input" =~ ^file- ]]; then
                    log_message "ERROR: Could not parse the input of job $job_id"
                    return 1

                else

                    # Get job status
                    local job_status
                    job_status=$(echo "$job_describe" | jq -r .state | sed 's/project-.*://' 2>/dev/null || echo "unknown")

                    # Get creation time
                    local job_created_ms
                    job_created_ms=$(echo "$job_describe" | jq -r .created | sed 's/project-.*://' 2>/dev/null || echo "0")

                    # Parse creation time
                    local job_created
                    job_created=$(date -u -d "@$(( job_created_ms / 1000))" 2>/dev/null || echo "unknown")

                fi


            fi

            # Check if input is indentical
            if [[ "$job_input" == "$VCF_LIST_HASH" ]]; then
                log_message "INFO: Found ongoing applet run with same input $VCF_LIST_HASH:"
                log_message "INFO:   Job ID: $job_id"
                log_message "INFO:   Job input: $job_input"
                log_message "INFO:   Created: $job_created"
                log_message "INFO:   Status: $job_status"

                case "$job_status" in
                  "terminated"|"failed"|"debug_hold")
                    echo "INFO: Previous job terminated or failed - Proceeding"
                    ;;
                  *)
                    log_message "ERROR: The previous job is running or finished without errors - Aborting to avoid duplicate processing"
                    return 1
                    ;;
                esac

            else
                continue
            fi
        fi
    done <<< "$ongoing_jobs"
    
    log_message "INFO: No ongoing runs using the same input"
    return 0
}
