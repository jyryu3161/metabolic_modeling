#** Identification of transcription factors responsible for metabolic reprogramming **
#Project
This project is to develop a framework that systematically predicts transcription factors responsible for metabolic reprogramming.

##Installation
###Major dependencies
- [gurobipy](http://www.gurobi.com/)


###Procedure
**Note**: This source code was developed in Linux, and has been tested in Ubuntu 14.04.5 LTS (i7-4770 CPU @ 3.40GHz)

1. Clone the repository

2. Create and activate virtual environment

        virtualenv venv
        source venv/bin/activate

3. Install packages at the root of the repository

        pip install pip --upgrade
        pip install -r requirements.txt

4. Install [gurobipy](http://www.gurobi.com/)

    In our case, we installed gurobipy in the root of a server, and created its symbolic link in `venv`:

        ln -s /usr/local/lib/python2.7/dist-packages/gurobipy/ $HOME/recon-manager/venv/lib/python2.7/site-packages/

###Implementation

For independent samples
```
python ./tf_metabolism.py -o ./output_sample -c1 ./input_data/sample1.csv -c2 ./input_data/sample2.csv
```

For related samples
```
python ./tf_metabolism.py -o ./output_sample -c1 ./input_data/sample1.csv -c2 ./input_data/sample2.csv -r
```

###Input arguments and corresponding files

- `-o` : Output directory

- `-c1` : RNA-Seq data for condition 1
    - See example (./input_data/sample1.csv)

- `-c2` : RNA-Seq data for condition 2
    - See example (./input_data/sample2.csv)
    
- `-r` : Flag for related samples
    - Use -r if samples are related
##Publication

