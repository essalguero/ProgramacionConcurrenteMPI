#!/bin/bash

export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH

#USER=Eugenio.Salguero
USER=eugenio
EXEC_FILES=convolutionMPI
DATA_FILES=

make clean
make

HOSTS=$(cat copy_host_file)

for HOST in $HOSTS
do
	#ssh $HOST rm -r ~/$USER/bin
	#ssh $HOST mkdir -p ~/$USER/bin
	ssh $HOST rm -r ~/$USER/bin/*
	scp $EXEC_FILES $HOST:~/$USER/bin
	scp host_file *.sh $HOST:~/$USER/bin
	ssh $HOST chmod 777 ~/$USER/bin/*
	if [ ! -z "$DATA_FILES" ]; then
		ssh $HOST mkdir -p ~/$USER/files
		scp $DATA_FILES $HOST:~/$USER/files
	fi
done

#mpirun --app app_file 
mpirun convolutionMPI $@
