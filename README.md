# UKBiobank WGS QC


## Testing the code


Login using your password `dx login` or using an admin token `dx login --token an_admin_token `.


Create a set of test input files to test the applet and upload them to `WGS/chr21/batches`.


```bash
# Create a target directory
dx mkdir -p "WGS/chr21/batches"

# Select five files from chr21
my_first_test_file=$( (seq 1 5 |\
	xargs -I{} echo "Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]/chr21/ukb23374_c21_b100{}_v1.vcf.gz" )|\
	dx upload - -o "/WGS/chr21/batches/vcf_list_test_1_of_3.input.txt" --brief)

```


Clone the QC repository locally, then compress the resources folder (contains the binaries, libraries and scripts) and upload them to the project folder on `dnanexus` platform


```bash

git clone https://github.com/mahmoudkoko/UKB_WGS_QC.git

cp -r UKB_WGS_QC/resources/usr wgs_qc_code


# Now add the first test file
echo -e "{\"qc_output_dir\": \"WGS/chr21/\", \"ukb23374_vcf_list\": {\"\$dnanexus_link\": \"${my_first_test_file}\"}}" > ./wgs_qc_code/job_input.json

zip -r wgs_qc_code.zip wgs_qc_code


```


Upload this code to the project directory


```bash

dx mkdir -p "Applets/test/"

my_applet_code=$(dx upload wgs_qc_code.zip --path "Applets/test/" --brief)

```


Start a cloud workstation to test the applet; choose high priority to prevent interruptions or delays


```bash
dx run app-cloud_workstation \
--destination "project-GzKk3XjJZz4ZgXzP2v1029qB:WGS/chr21/" \
--priority high \
--input fids="$my_applet_code" \
--input max_session_length="1h" \
--instance-type mem1_ssd1_v2_x4 \
--ssh

```

Once logged in, you will see that the file we provided was downloaded. Do the following to emulate the run enviroment of the applet.

```bash
unzip wgs_qc_code.zip

# copy the binaries to /usr/bin/
sudo cp wgs_qc_code/bin/* /usr/bin/

# copy the libraries
sudo cp wgs_qc_code/lib/* /usr/lib/

# copy the scripts
sudo cp -r wgs_qc_code/scripts /usr/scripts

# replace the job input json
cat wgs_qc_code/job_input.json > /home/dnanexus/job_input.json

# install parallel
#sudo apt install parallel
```

Starts by sourcing the scripts and preparing the environemnts:


```bash
source /usr/scripts/export_global_vars.sh
source /usr/scripts/check_ongoing_applet_runs.sh
source /usr/scripts/process_graphtyper_vcf_block.sh
source /usr/scripts/global_cleanup.sh
source /usr/scripts/utilities.sh

# set the env
export_global_vars

export -f process_graphtyper_vcf_block
export -f log_message
export -f track_temp_dir
export -f track_temp_file


# Set up global cleanup trap
trap global_cleanup EXIT INT TERM
```


Now see if the applet can regonize the correct variables (should show your project and input parameters)
```bash

# Execution context
log_execution_context

```

If this prints correctly, proceed to downloading the input file.


```bash

# Initiate dirs and files
log_dir="./LOG_DIR"
jobsummary_dx_id=""

joblog="${log_dir}/${VCF_LIST_NAME%.input.txt}.log"
jobsummary="${log_dir}/${VCF_LIST_NAME%.input.txt}.summary.txt.gz"

track_temp_file "$joblog"
track_temp_file "$jobsummary"    

# Create log directory
mkdir -p "$log_dir"

    
track_temp_dir "$log_dir"
    

# Process the input file
# Process the input file

log_message "INFO: Getting the VCF list ..."

# Check for ongoing applet runs with same input
if ! check_ongoing_applet_runs ; then
    log_message "ERROR: Duplicate applet run detected - exiting"

elif [[ ! "${VCF_LIST_NAME}" =~ ^vcf_list_ ]] || [[ ! ${VCF_LIST_NAME} =~ .input.txt$ ]] ; then

    # If the name isn't right, issue an error and exit
    log_message "ERROR: $VCF_LIST_NAME does not match expected pattern (vcf_list_n_of_N.input.txt) - exiting"

    # Otherwise download the file
elif ! dx download "${DX_PROJECT_CONTEXT_ID}:${VCF_LIST_HASH}" -o "${log_dir}/${VCF_LIST_NAME}"; then

    # If the download fails, issue an error and exit
    log_message "ERROR: Failed to download VCF list file - exiting"

else

    log_message "INFO: List contains $(cat ${log_dir}/${VCF_LIST_NAME} | wc -l) elements"
fi

    
```

If the code identifies and downloads the input list correctly, proceed to testing the qc function.

```bash
test_file=$(head -n1 ${log_dir}/${VCF_LIST_NAME})

process_graphtyper_vcf_block $test_file
```




## Building the applet



Create a directory to keep your applets organized, then build the applet in this directory and capture the applet ID


```bash

my_applet_id=$(dx build -f . -d "Applets/test/" --brief | jq -r .id)

```



```bash

dx run ${my_applet_id} \
--instance-type mem1_ssd1_v2_x72 \
--destination "project-GzKk3XjJZz4ZgXzP2v1029qB:WGS/chr21/batches/" \
--input ukb23374_vcf_list=${my_first_test_file} \
--input qc_output_dir="WGS/chr21/" \
--tag "wgs_qc_batch" \
--name "WGS QC: chr21 - test 1 of 3" \
--priority low \
--debug-on All \
--allow-ssh \
--brief

```



```bash
# Create a target directory
dx mkdir -p "WGS/chr21/batches"


# Five files
my_first_test_file=$(seq 1 5 | xargs -I{} echo '"Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]/chr21/ukb23374_c21_b100{}_v1.vcf.gz"' |\
dx upload - -o "/WGS/chr21/batches/vcf_list_test_1_of_3.input.txt" --brief)

# 45 files
my_second_test_file=$(seq 6 50 | xargs -I{} echo '"Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]/chr21/ukb23374_c21_b100{}_v1.vcf.gz"' |\
dx upload - -o "/WGS/chr21/batches/vcf_input_list_2_of_3.input.txt" --brief)

# 450 files
my_third_test_file=$(seq 51 500 | xargs -I{} echo '"Bulk/GATK and GraphTyper WGS/GraphTyper population level WGS variants, pVCF format [500k release]/chr21/ukb23374_c21_b100{}_v1.vcf.gz"' |\
dx upload - -o "/WGS/chr21/batches/vcf_list_test_3_of_3.input.txt" --brief)


```




First test with a few files (~5 files). Invoke ssh access and debug mode if e.g., you would like to check the enviroment variables to sort out errors or understand some aspects of the code.

Note the instance is rather large to allow faster scheduling; the cost is minimal for a few files.


```bash

dx run ${my_applet_id} \
--instance-type mem1_ssd1_v2_x72 \
--destination "project-GzKk3XjJZz4ZgXzP2v1029qB:WGS/chr21/batches/" \
--input ukb23374_vcf_list=${my_first_test_file} \
--input qc_output_dir="WGS/chr21/" \
--tag "wgs_qc_batch" \
--name "WGS QC: chr21 - test 1 of 3" \
--priority low \
--debug-on All \
--allow-ssh \
--brief

```
