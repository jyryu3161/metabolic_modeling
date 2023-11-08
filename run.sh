# python3 tf_metabolism.py -o output -c1 ./input_data/sample1.csv -c2 ./input_data/sample2.csv


python3 tf_metabolism.py -o output_aging_old -c1 ./input_data/young_group_rnaseq.csv -c2 ./input_data/old_age.csv
python3 tf_metabolism.py -o output_aging_mid -c1 ./input_data/young_group_rnaseq.csv -c2 ./input_data/mid_group_rnaseq.csv
python3 tf_metabolism.py -o output_NASH -c1 ./input_data/NASH_control.csv -c2 ./input_data/NASH_test.csv
