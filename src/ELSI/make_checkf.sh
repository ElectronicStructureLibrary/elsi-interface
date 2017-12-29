#!/bin/bash
# VY: This is a simple bash script used for "make checkf"

rm rw_real_serial.log 2> /dev/null
rm rw_complex_serial.log 2> /dev/null
rm rw_real_parallel.log 2> /dev/null
rm rw_complex_parallel.log 2> /dev/null
rm ev_real_serial.log 2> /dev/null
rm ev_real_serial_timings.json 2> /dev/null
rm ev_complex_serial.log 2> /dev/null
rm ev_complex_serial_timings.json 2> /dev/null
rm ev_real_elpa.log 2> /dev/null
rm ev_real_elpa_timings.json 2> /dev/null
rm ev_complex_elpa.log 2> /dev/null
rm ev_complex_elpa_timings.json 2> /dev/null
rm ev_real_sparse_elpa.log 2> /dev/null
rm ev_real_sparse_elpa_timings.json 2> /dev/null
rm ev_complex_sparse_elpa.log 2> /dev/null
rm ev_complex_sparse_elpa_timings.json 2> /dev/null
rm dm_real_elpa.log 2> /dev/null
rm dm_real_elpa_timings.json 2> /dev/null
rm dm_complex_elpa.log 2> /dev/null
rm dm_complex_elpa_timings.json 2> /dev/null
rm dm_real_sparse_elpa.log 2> /dev/null
rm dm_real_sparse_elpa_timings.json 2> /dev/null
rm dm_complex_sparse_elpa.log 2> /dev/null
rm dm_complex_sparse_elpa_timings.json 2> /dev/null
rm dm_real_libomm.log libOMM.log 2> /dev/null
rm dm_real_libomm_timings.json 2> /dev/null
rm dm_complex_libomm.log libOMM.log 2> /dev/null
rm dm_complex_libomm_timings.json 2> /dev/null
rm dm_real_sparse_libomm.log 2> /dev/null
rm dm_real_sparse_libomm_timings.json 2> /dev/null
rm dm_complex_sparse_libomm.log 2> /dev/null
rm dm_complex_sparse_libomm_timings.json 2> /dev/null
rm dm_real_pexsi.log logPEXSI0 2> /dev/null
rm dm_real_pexsi_timings.json 2> /dev/null
rm dm_complex_pexsi.log logPEXSI0 2> /dev/null
rm dm_complex_pexsi_timings.json 2> /dev/null
rm dm_real_sparse_pexsi.log 2> /dev/null
rm dm_real_sparse_pexsi_timings.json 2> /dev/null
rm dm_complex_sparse_pexsi.log 2> /dev/null
rm dm_complex_sparse_pexsi_timings.json 2> /dev/null
set -e # Stop on error

RED_ALART="false"
echo
echo "Test program output may be found in $PWD"

echo
echo -n "Running the 'serial reading/writing real matrices' Fortran test"
${MPI_EXEC} -n 1 ./test_rw_real.x ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc > rw_real_serial.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./rw_real_serial.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/rw_real_serial.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'serial reading/writing complex matrices' Fortran test"
${MPI_EXEC} -n 1 ./test_rw_complex.x ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc > rw_complex_serial.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./rw_complex_serial.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/rw_complex_serial.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel reading/writing real matrices' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_rw_real.x ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc > rw_real_parallel.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./rw_real_parallel.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/rw_real_parallel.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel reading/writing complex matrices' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_rw_complex.x ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc > rw_complex_parallel.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./rw_complex_parallel.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/rw_complex_parallel.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'serial real eigensolver' Fortran test"
${MPI_EXEC} -n 1 ./test_ev_real.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > ev_real_serial.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv ev_real_timings.json ev_real_serial_timings.json 2>/dev/null

if (! grep -q "Passed" <./ev_real_serial.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_real_serial.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'serial complex eigensolver' Fortran test"
${MPI_EXEC} -n 1 ./test_ev_complex.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > ev_complex_serial.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv ev_complex_timings.json ev_complex_serial_timings.json 2>/dev/null

if (! grep -q "Passed" <./ev_complex_serial.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_complex_serial.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real eigensolver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_real.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > ev_real_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv ev_real_timings.json ev_real_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./ev_real_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_real_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex eigensolver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_complex.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > ev_complex_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv ev_complex_timings.json ev_complex_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./ev_complex_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_complex_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real sparse eigensolver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_real_sparse.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > ev_real_sparse_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv ev_real_sparse_timings.json ev_real_sparse_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./ev_real_sparse_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_real_sparse_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex sparse eigensolver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_complex_sparse.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > ev_complex_sparse_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv ev_complex_sparse_timings.json ev_complex_sparse_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./ev_complex_sparse_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_complex_sparse_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

if [ "$ENABLE_SIPS" = "yes" ]
then
   echo
   echo -n "Running the 'parallel real eigensolver + SIPs' Fortran test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_real.x 5 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > ev_real_sips.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done

   if (! grep -q "Passed" <./ev_real_sips.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/ev_real_sips.log for details."
   else
      tput setaf 10
      echo " PASSED!"
      tput sgr0
   fi
fi

echo
echo -n "Running the 'parallel real density matrix solver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv dm_real_timings.json dm_real_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_real_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_real_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex density matrix solver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv dm_complex_timings.json dm_complex_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_complex_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_complex_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real sparse density matrix solver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real_sparse.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_sparse_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv dm_real_sparse_timings.json dm_real_sparse_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_real_sparse_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_real_sparse_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex sparse density matrix solver + ELPA' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex_sparse.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_sparse_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done
mv dm_complex_sparse_timings.json dm_complex_sparse_elpa_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_complex_sparse_elpa.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_complex_sparse_elpa.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real density matrix solver + libOMM' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real.x 2 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_libomm.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done
mv dm_real_timings.json dm_real_libomm_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_real_libomm.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_real_libomm.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex density matrix solver + libOMM' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex.x 2 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_libomm.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done
mv dm_complex_timings.json dm_complex_libomm_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_complex_libomm.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_complex_libomm.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real sparse density matrix solver + libOMM' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real_sparse.x 2 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_sparse_libomm.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done
mv dm_real_sparse_timings.json dm_real_sparse_libomm_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_real_sparse_libomm.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_real_sparse_libomm.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex sparse density matrix solver + libOMM' Fortran test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex_sparse.x 2 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_sparse_libomm.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done
mv dm_complex_sparse_timings.json dm_complex_sparse_libomm_timings.json 2>/dev/null

if (! grep -q "Passed" <./dm_complex_sparse_libomm.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_complex_sparse_libomm.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

if [ "$DISABLE_PEXSI" != "yes" ]
then
   echo
   echo -n "Running the 'parallel real density matrix solver + PEXSI' Fortran test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real.x 3 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_pexsi.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done
   mv dm_real_timings.json dm_real_pexsi_timings.json 2>/dev/null

   if (! grep -q "Passed" <./dm_real_pexsi.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/dm_real_pexsi.log for details."
   else
      tput setaf 10
      echo " PASSED!"
      tput sgr0
   fi

   echo
   echo -n "Running the 'parallel complex density matrix solver + PEXSI' Fortran test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex.x 3 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_pexsi.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done
   mv dm_complex_timings.json dm_complex_pexsi_timings.json 2>/dev/null

   if (! grep -q "Passed" <./dm_complex_pexsi.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/dm_complex_pexsi.log for details."
   else
      tput setaf 10
      echo " PASSED!"
      tput sgr0
   fi

   echo
   echo -n "Running the 'parallel real sparse density matrix solver + PEXSI' Fortran test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real_sparse.x 3 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_sparse_pexsi.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done
   mv dm_real_sparse_timings.json dm_real_sparse_pexsi_timings.json 2>/dev/null

   if (! grep -q "Passed" <./dm_real_sparse_pexsi.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/dm_real_sparse_pexsi.log for details."
   else
      tput setaf 10
      echo " PASSED!"
      tput sgr0
   fi

   echo
   echo -n "Running the 'parallel complex sparse density matrix solver + PEXSI' Fortran test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex_sparse.x 3 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_sparse_pexsi.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done
   mv dm_complex_sparse_timings.json dm_complex_sparse_pexsi_timings.json 2>/dev/null

   if (! grep -q "Passed" <./dm_complex_sparse_pexsi.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/dm_complex_sparse_pexsi.log for details."
   else
      tput setaf 10
      echo " PASSED!"
      tput sgr0
   fi

else
   tput setaf 5
   echo
   echo "PEXSI has been disabled by user's choice."
   tput sgr0
fi

echo
if [ "$RED_ALART" = "true" ]
then
   tput setaf 5
   echo "MAKE CHECK FAILED, CHECK YOUR COMPILATION SETTINGS!"
   tput sgr0
   exit 1
else
   tput setaf 10
   echo "MAKE CHECK PASSED, YOU'RE GOOD TO GO!"
   tput sgr0
   exit 0
fi
