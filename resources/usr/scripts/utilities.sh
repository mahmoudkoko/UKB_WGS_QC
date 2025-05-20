###################
# Utility functions
###################


# Function to print STDERR message with timestamp
log_message() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" >&2
}


# Function to safely add directories to cleanup list
track_temp_dir() {
    local dir="$1"
    GLOBAL_TEMP_DIRS+=("$dir")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: Tracking temp directory: $dir" >&2
}

# Function to safely add files to cleanup list
track_temp_file() {
    local file="$1"
    GLOBAL_TEMP_FILES+=("$file")
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] INFO: Tracking temp file: $file" >&2
}

