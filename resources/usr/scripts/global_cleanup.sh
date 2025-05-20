###################
# Cleanup functions
###################

global_cleanup() {
    local exit_code=$?
    
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: Starting cleanup ..." >&2
    
    # Kill any remaining parallel jobs
    pkill -P $$ 2>/dev/null || true
    
    # Clean up any remaining temporary directories
    for dir in "${GLOBAL_TEMP_DIRS[@]}"; do
        if [[ -d "$dir" ]]; then
            rm -rf "$dir" 2>/dev/null || true
        fi
    done
    
    # Clean up any remaining temporary files
    for file in "${GLOBAL_TEMP_FILES[@]}"; do
        if [[ -f "$file" ]]; then
            rm -f "$file" 2>/dev/null || true
        fi
    done
    
    # Log final status
    if [[ $exit_code -eq 0 ]]; then
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: Cleanup completed successfully" >&2
    else
        echo "[$(date '+%Y-%m-%d %H:%M:%S')] WARNING: Cleanup completed after error (exit code: $exit_code)" >&2
    fi
    
    exit $exit_code
}

