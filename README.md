#** Identification of transcription factors responsible for metabolic reprogramming **
#Project
This project is to develop a framework that systematically predicts transcription factors responsible for metabolic reprogramming.

##Installation
###Major dependencies
- [gurobipy](http://www.gurobi.com/)


###Procedure
**Note**: This source code was developed in Linux, and has been tested in Ubuntu 20.04.5 LTS

1. Clone the repository

2. Create and activate conda environment

        conda env create -f environment.yml
        conda activate gems

3. Activate gurobi license

###Implementation

```
python3 tf_metabolism.py -o output -c1 ./input_data/sample1.csv -c2 ./input_data/sample2.csv
```


###Input arguments and corresponding files

- `-o` : Output directory

- `-c1` : RNA-Seq data for condition 1
    - See example (./input_data/sample1.csv)

- `-c2` : RNA-Seq data for condition 2
    - See example (./input_data/sample2.csv)
    
##Publication

