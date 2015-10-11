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
mpirun -n 4 ../../bin/construct_and_write_omm 100 10 > /dev/null
fail=`h5diff elsi_eigenvalue_problem.hdf5 data/elsi_eigenvalue_problem.hdf5`
if [ "x$fail" != "x" ]; then
   echo "construct_and_write_omm failed"
else
   echo "construct_and_write_omm successfull"
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
mpirun -n 4 ../../bin/set_and_write_omm 100 10 > /dev/null
fail=`h5diff elsi_eigenvalue_problem.hdf5 data/elsi_eigenvalue_problem.hdf5`
if [ "x$fail" != "x" ]; then
   echo "set_and_write_omm failed"
else
   echo "set_and_write_omm successfull"
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
mpirun -n 4 ../../bin/read_and_write_omm > /dev/null
fail=`h5diff elsi_eigenvalue_problem_out.hdf5 data/elsi_eigenvalue_problem_out.hdf5`
if [ "x$fail" != "x" ]; then
   echo "read_and_write_omm failed"
else
   echo "read_and_write_omm successfull"
   rm elsi_eigenvalue_problem_out.hdf5
fi
cd ..

# Read and Solve
cd read_and_solve
mpirun -n 4 ../../bin/read_and_solve_elpa 10 > buffer.dat
grep "total energy" buffer.dat > e_tot.dat
rm buffer.dat
fail=`diff e_tot.dat data/e_tot.dat`
if [ "x$fail" != "x" ]; then
   echo "read_and_solve_elpa failed"
else
   echo "read_and_solve_elpa successfull"
   rm eigenvalues.dat
fi
mpirun -n 4 ../../bin/read_and_solve_omm 10 > buffer.dat
grep "total energy" buffer.dat > e_tot.dat
rm buffer.dat
fail=`diff e_tot.dat data/e_tot.dat`
if [ "x$fail" != "x" ]; then
   echo "read_and_solve_omm failed"
else
   echo "read_and_solve_omm successfull"
   rm eigenvalues.dat
fi
cd ..
