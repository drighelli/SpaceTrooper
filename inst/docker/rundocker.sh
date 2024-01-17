docker run \
    'mounting volume inside a docker i.e. for input/output data'
    'in order to write here the user in the docker need the writing authorization on your local folder'
    '-v <your_local_dir>:<link_in_docker> \'
    -it \
    --user rstudio \
    drighelli/spatialmoleculeanalysis:1.0 bash
