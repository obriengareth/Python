#module load anaconda/3.7
#/apps/python/anaconda3.7/bin/mpiexec -machinefile nodes.txt -n 4 -N 1 python hello_world_mpi.py
#$(which mpiexec) -machinefile nodes.txt -n 4 -N 1 python hello_world_mpi.py
###############################################################################
###############################################################################

import sys
from mpi4py import MPI

###############################################################################
###############################################################################

# MPI basic information
comm = MPI.COMM_WORLD
myrank = comm.Get_rank()
mysize = comm.Get_size()
name = MPI.Get_processor_name() 

print("Hello World from Rank %d of %d, node %s" % (myrank, mysize, name))

###############################################################################
###############################################################################

#size = MPI.COMM_WORLD.Get_size()
#rank = MPI.COMM_WORLD.Get_rank()
#srnk = MPI.COMM_WORLD.Split(color=int(name[-2:]), key=rank).Get_rank() # unimplemented Split_type
#srnk = MPI.COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED).Get_rank()
#srnk = MPI.COMM_WORLD.Split_type(MPI.COMM_TYPE_SHARED).Get_rank()
#print("Rank %d of %d, node %s" % (myrank, mysize, name))

###############################################################################
###############################################################################
