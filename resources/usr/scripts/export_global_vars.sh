#########################
# Global cleanup function
#########################

# Function to export global env variables

export_global_vars() {

# Going forward, CAPs will indicate variables not specific to a function.
export DX_APPLET_ID=$(cat /home/dnanexus/dnanexus-job.json | jq -r .executable | sed 's/project-.*://' )
export DX_EXECUTABLE_NAME=$(cat /home/dnanexus/dnanexus-job.json | jq -r .executableName | sed 's/project-.*://' )
export DX_EXECUTION_DIR=$(cat /home/dnanexus/dnanexus-job.json | jq -r .folder | sed 's/project-.*://' )

# This will export the number of parallel VCF QC jobs to run
export DX_CPUS=$(cat /home/dnanexus/dnanexus-executable.json | jq -r '.runSpec.systemRequirements."*".instanceType' | sed 's/mem.*x//')

if [[ DX_CPUS -gt 4 ]]; then 
	export N_JOBS=$(( DX_CPUS / 4 ))
else
	export N_JOBS=1
fi 


# This will export the input vcf hash, name and directory.
export VCF_LIST_HASH="$(cat job_input.json  | jq -r '.ukb23374_vcf_list.["$dnanexus_link"]' | sed 's/project-.*://' )"
export VCF_LIST_NAME=$(dx describe --json $VCF_LIST_HASH | jq -r .name | sed 's/project-.*://' )
export VCF_LIST_DIR=$(dx describe --json $VCF_LIST_HASH | jq -r .folder | sed 's/project-.*://' )

# This will export the output directory (where the QC data should be stored)
export OUTPUT_DIR="$(cat job_input.json  | jq -r .qc_output_dir | sed 's/project-.*://' )"


# BCFtools, bgzip, and tabix binaries are pre-packaged in this applet under /usr/bin which is already in $PATH. 
# This will export the path to some required libraries to run bcftools or its plugins
export LD_LIBRARY_PATH=/usr/lib
export BCFTOOLS_PLUGINS=/usr/lib

# Global cleanup tracking
declare -a GLOBAL_TEMP_DIRS=()
declare -a GLOBAL_TEMP_FILES=()

}


# Function to log execution context on STDOUT (log file)

log_execution_context() {
    echo "=== Execution Context ==="
    echo "Job ID: ${DX_JOB_ID:-'Not available'}"
    echo "Project ID: ${DX_PROJECT_CONTEXT_ID:-'Not available'}"
    echo "Applet ID: ${DX_APPLET_ID:-'Not available'}"
    echo "Applet Name: ${DX_EXECUTABLE_NAME:-'Not available'}"
    echo "Destination: ${DX_EXECUTION_DIR:-'Not available'}"
    echo "=== Inputs ==="
    echo "Directory: ${VCF_LIST_DIR:-'Not available'}"
    echo "File name: ${VCF_LIST_NAME:-'Not available'}"
    echo "File hash: ${VCF_LIST_HASH:-'Not available'}"
    echo "=== Output directories ==="
    echo "QC files (bgen, sites, stats) : ${OUTPUT_DIR:-'Not available'}"
    echo "Log files (summary.txt.gz): ${VCF_LIST_DIR:-'Not available'}"
}
