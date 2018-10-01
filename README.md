# Local Protein Sequence Design
This repository provides PyRosetta based  methods for designing a local part of a protein.

## Installation
I recommand to create a python virtual environment before installing the package.
```
virtualenv -p /usr/bin/python3 --system-site-packages venv
```
Activate the virtual environment and install the package by
```
pip install -e .
```

## Dependencies
Several external programs are required for filtering the designs. The list of dependencies can be found [here](https://github.com/xingjiepan/local_protein_sequence_design/blob/master/site_settings/site_settings.template.json). You need to create a `site_settings.json` file inside the `site_settings` directory and specify paths to all the dependencies.

## Running An Example
This is an example for designing sequences that stabilize the a backbone generated by the [loop helix loop reshaping](https://github.com/Kortemme-Lab/loop_helix_loop_reshaping) method.

### First iteration of design
Let's start by doing a fast and dirty round of design. The input backbones and insertion point definition files are in the `test_inputs/two_lhl_units_2lv8` directory. Do determine the designable and repackable residue, you need also to know the structure before the LHL reshaping. This information is provided by the files `test_inputs/2lv8_inputs/2lv8_cleaned.pdb` and `test_inputs/2lv8_inputs/2lv8_insertion_points.json`. Run the design script
```
./run_jobs.py sequence_design_for_LHL_reshaping_example job_scripts/sequence_design_for_LHL_reshaping_example.py
```
This command will create the directory `data/sequence_design_for_LHL_reshaping_example` and dump the designed structure into it.

### Filter the designs
After making the design, we want to evaluate the design by a set of metrics. To do so, run
```
./run_jobs.py sequence_design_for_LHL_reshaping_example job_scripts/filter.py
```
This command will create a `filter_info.json` file for each of the design. To pull the filter scores of all designs into a tsv table, run
```
./job_scripts/aggregate_data_into_a_table.py data/sequence_design_for_LHL_reshaping_example
```
This command will create a tsv file `data/sequence_design_for_LHL_reshaping_example/summary_table.tsv`. Now you may want to select a subset of designs for next iteration of design. In the demo, we only made one design. So we will use this design for the subsequent iteration. Copy the summary table by
```
cp data/sequence_design_for_LHL_reshaping_example/summary_table.tsv data/sequence_design_for_LHL_reshaping_example/selected_designs_for_next_iteration.tsv
```

### Second iteration of design
Do the second round of design using the result of the first round of design as the input. Run
```
./run_jobs.py sequence_design_for_LHL_reshaping_iteration_2_example job_scripts/sequence_design_for_LHL_reshaping_iteration_2_example.py
```
This command will create the directory `data/sequence_design_for_LHL_reshaping_examples` and dump the designed structure into it. Note that the second iteration of design will be finished in much longer time the the first iteration, because we are using `ex1 ex2` rotamers in this round.

### Select good designs
In this example, we only make 1 design, so there is no point select good designs. But in a real problem, you want to do the selection. First filter the designs using the command described above. Then you may want to select the good designs using the script `job_scripts/select_designs.py`.

## Production Run
For a real problem, copy the example scripts to the `job_scripts/user` directory and edit them for inputs and options. Run a job script by
```
./run_jobs.py my_output_data_set job_scripts/user/my_script.py
```
A real design problem is usually too computationaly expensive to ran on a single computer. You'd better use a cluster. To submit a job to the UCSF QB3 cluster. Just do
```
./run_jobs.py my_output_data_set job_scripts/user/my_script.py -d SGE -n 100
```
This command will use 100 CPUs to run the job.
