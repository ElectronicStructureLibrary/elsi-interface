#!/bin/sh

# Construct and Write
cd construct_and_write
mpirun -n 10 ../../bin/construct_and_write_elpa 100 16 > /dev/null
fail=`h5diff elsi_eigenvalue_problem.hdf5 data/elsi_eigenvalue_problem.hdf5`
if [ "x$fail" != "x" ]; then
   echo "construct_and_write_elpa failed"
else
   echo "construct_and_write_elpa successfull"
   rm elsi_eigenvalue_problem.hdf5
fi
mpirun -n 10 ../../bin/construct_and_write_omm 100 16 > /dev/null
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
mpirun -n 10 ../../bin/set_and_write_elpa 100 16 > /dev/null
fail=`h5diff elsi_eigenvalue_problem.hdf5 data/elsi_eigenvalue_problem.hdf5`
if [ "x$fail" != "x" ]; then
   echo "set_and_write_elpa failed"
else
   echo "set_and_write_elpa successfull"
   rm elsi_eigenvalue_problem.hdf5
fi
mpirun -n 10 ../../bin/set_and_write_omm 100 16 > /dev/null
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
mpirun -n 10 ../../bin/read_and_write_elpa > /dev/null
fail=`h5diff elsi_eigenvalue_problem_out.hdf5 data/elsi_eigenvalue_problem_out.hdf5`
if [ "x$fail" != "x" ]; then
   echo "read_and_write_elpa failed"
else
   echo "read_and_write_elpa successfull"
   rm elsi_eigenvalue_problem_out.hdf5
fi

mpirun -n 10 ../../bin/read_and_write_omm > /dev/null
fail=`h5diff elsi_eigenvalue_problem_out.hdf5 data/elsi_eigenvalue_problem_out.hdf5`
if [ "x$fail" != "x" ]; then
   echo "read_and_write_omm failed"
else
   echo "read_and_write_omm successfull"
   rm elsi_eigenvalue_problem_out.hdf5
fi

mpirun -n 10 ../../bin/read_and_write_pexsi > /dev/null
fail=`h5diff elsi_eigenvalue_problem_out.hdf5 data/elsi_eigenvalue_problem_out.hdf5`
if [ "x$fail" != "x" ]; then
   echo "read_and_write_pexsi failed"
else
   echo "read_and_write_pexsi successfull"
   rm elsi_eigenvalue_problem_out.hdf5
fi

cd ..

# Read and Solve
cd read_and_solve
mpirun -n 10 ../../bin/read_and_solve_elpa 360 > buffer.dat
grep "total energy" buffer.dat > e_tot.dat
rm buffer.dat
etot=`cat e_tot.dat | awk '{print $4}' | sed 's/[eE]+/\*10\^/'`
etot_ref=`cat data/e_tot-elpa.dat | awk '{print $4}' | sed 's/[eE]+/\*10\^/'`
diff=`echo "sqrt(($etot - $etot_ref)^2)" | bc -l` 
small=`echo "$diff < 1*10^-5" | bc -l`
if [ $small == "1" ]; then
   echo "read_and_solve_elpa successfull"
   rm e_tot.dat
else
   echo "read_and_solve_elpa failed, difference $diff"
fi

mpirun -n 10 ../../bin/read_and_solve_omm 360 > buffer.dat
grep "total energy" buffer.dat > e_tot.dat
rm buffer.dat
etot=`cat e_tot.dat | awk '{print $4}' | sed 's/[eE]+/\*10\^/'`
etot_ref=`cat data/e_tot-omm.dat | awk '{print $4}' | sed 's/[eE]+/\*10\^/'`
diff=`echo "sqrt(($etot - $etot_ref)^2)" | bc -l` 
small=`echo "$diff < 1*10^-4" | bc -l`
if [ $small == "1" ]; then
   echo "read_and_solve_omm successfull"
   rm e_tot.dat
else
   echo "read_and_solve_omm failed, difference $diff, $etot, $etot_ref"
fi

mpirun -n 10 ../../bin/read_and_solve_pexsi 360 > buffer.dat
grep "total energy" buffer.dat > e_tot.dat
rm buffer.dat
etot=`cat e_tot.dat | awk '{print $4}' | sed 's/[eE]+/\*10\^/'`
etot_ref=`cat data/e_tot-pexsi.dat | awk '{print $4}' | sed 's/[eE]+/\*10\^/'`
diff=`echo "sqrt(($etot - $etot_ref)^2)" | bc -l` 
small=`echo "$diff < 1*10^-4" | bc -l`
if [ $small == "1" ]; then
   echo "read_and_solve_pexsi successfull"
   rm e_tot.dat
else
   echo "read_and_solve_pexsi failed, difference $diff, $etot, $etot_ref"
fi

cd ..
