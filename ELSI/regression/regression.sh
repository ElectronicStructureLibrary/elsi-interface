#!/bin/sh

# Construct and Write
cd construct_and_write
mpirun -n 4 ../../bin/construct_and_write_elpa 100 10 > /dev/null
fail=`h5diff elsi_eigenvalue_problem.hdf5 data/elsi_eigenvalue_problem.hdf5`
if [ "x$fail" != "x" ]; then
   echo "construct_and_write_elpa failed"
else
   echo "construct_and_write_elpa successfull"
   rm elsi_eigenvalue_problem.hdf5
fi
cd ..

# Set and Write
cd set_and_write
mpirun -n 4 ../../bin/set_and_write_elpa 100 10 > /dev/null
fail=`h5diff elsi_eigenvalue_problem.hdf5 data/elsi_eigenvalue_problem.hdf5`
if [ "x$fail" != "x" ]; then
   echo "set_and_write_elpa failed"
else
   echo "set_and_write_elpa successfull"
   rm elsi_eigenvalue_problem.hdf5
fi
cd ..

# Read and Write
cd read_and_write
mpirun -n 4 ../../bin/read_and_write_elpa > /dev/null
fail=`h5diff elsi_eigenvalue_problem_out.hdf5 data/elsi_eigenvalue_problem_out.hdf5`
if [ "x$fail" != "x" ]; then
   echo "read_and_write_elpa failed"
else
   echo "read_and_write_elpa successfull"
   rm elsi_eigenvalue_problem_out.hdf5
fi
cd ..

# Read and Solve
cd read_and_solve
mpirun -n 4 ../../bin/read_and_solve_elpa 100 > buffer.dat
tail -n 100 buffer.dat > eigenvalues.dat
rm buffer.dat
fail=`diff eigenvalues.dat data/eigenvalues.dat`
if [ "x$fail" != "x" ]; then
   echo "read_and_solve_elpa failed"
else
   echo "read_and_solve_elpa successfull"
   rm eigenvalues.dat
fi
cd ..
