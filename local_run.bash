export OMP_NUM_THREADS=4  # or whatever number of threads 

OUTPUT_FILE='output_files/argon09/out.argon09_t_075'
XYZ_FILE='argon09_t_075.xyz'
lmp -in input_files/in.argon09_t_075 -log $OUTPUT_FILE
python Properties_run.py --filename $OUTPUT_FILE --output_dir output_files/argon09
python RDF.py --filename $XYZ_FILE --output_dir output_files/argon09 --box_size  6.665  6.665  6.665 

OUTPUT_FILE='output_files/argon1_0/out.argon1_t_075'
XYZ_FILE='argon1_t_075.xyz'
lmp -in input_files/in.argon1_t_075 -log output_files/argon1_0/out.argon1_t_075
python Properties_run.py --filename $OUTPUT_FILE --output_dir output_files/argon1_0
python RDF.py --filename $XYZ_FILE --output_dir output_files/argon1_0 --box_size  6.433  6.433  6.433 

OUTPUT_FILE='output_files/argon1_1/out.argon1_1_t_075'
XYZ_FILE='argon1_1_t_075.xyz'
lmp -in input_files/in.argon1_1_t_075 -log output_files/argon1_1/out.argon1_1_t_075
python Properties_run.py --filename $OUTPUT_FILE --output_dir output_files/argon1_1
python RDF.py --filename $XYZ_FILE --output_dir output_files/argon1_1 --box_size  6.228  6.228  6.228