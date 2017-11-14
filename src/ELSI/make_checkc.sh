#!/bin/bash
# VY: This is a simple bash script used for "make checkc"

if [ "$C_INTERFACE" = "yes" ]
then

rm ev_real_serial_c.log 2> /dev/null
rm ev_complex_serial_c.log 2> /dev/null
rm ev_real_elpa_c.log 2> /dev/null
rm ev_complex_elpa_c.log 2> /dev/null
rm dm_real_elpa_c.log 2> /dev/null
rm dm_complex_elpa_c.log 2> /dev/null
rm dm_real_libomm_c.log libOMM.log 2> /dev/null
rm dm_complex_libomm_c.log libOMM.log 2> /dev/null
rm dm_real_pexsi_c.log logPEXSI0 2> /dev/null
rm dm_complex_pexsi_c.log logPEXSI0 2> /dev/null
set -e # Stop on error

RED_ALART="false"
echo
echo "Test program output may be found in $PWD"

echo
echo -n "Running the 'serial real eigensolver' C test"
${MPI_EXEC} -n 1 ./test_ev_real_c.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > ev_real_serial_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./ev_real_serial_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_real_serial_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'serial complex eigensolver' C test"
${MPI_EXEC} -n 1 ./test_ev_complex_c.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > ev_complex_serial_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./ev_complex_serial_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_complex_serial_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real eigensolver + ELPA' C test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_real_c.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > ev_real_elpa_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./ev_real_elpa_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_real_elpa_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex eigensolver + ELPA' C test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_ev_complex_c.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > ev_complex_elpa_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./ev_complex_elpa_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/ev_complex_elpa_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real density matrix solver + ELPA' C test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real_c.x 1 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 1 > dm_real_elpa_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./dm_real_elpa_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_real_elpa_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex density matrix solver + ELPA' C test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex_c.x 1 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 1 > dm_complex_elpa_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

if (! grep -q "Passed" <./dm_complex_elpa_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_complex_elpa_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel real density matrix solver + libOMM' C test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real_c.x 2 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 2 > dm_real_libomm_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done

if (! grep -q "Passed" <./dm_real_libomm_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_real_libomm_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

echo
echo -n "Running the 'parallel complex density matrix solver + libOMM' C test"
${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex_c.x 2 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 2 > dm_complex_libomm_c.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done

if (! grep -q "Passed" <./dm_complex_libomm_c.log); then
   tput setaf 5
   RED_ALART="true"
   echo " FAILED!"
   tput sgr0
   echo "See `pwd`/dm_complex_libomm_c.log for details."
else
   tput setaf 10
   echo " PASSED!"
   tput sgr0
fi

if [ "$DISABLE_PEXSI" != "yes" ]
then
   echo
   echo -n "Running the 'parallel real density matrix solver + PEXSI' C test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_real_c.x 3 ${ELSI_DIR}/test/H_real.csc ${ELSI_DIR}/test/S_real.csc 3 > dm_real_pexsi_c.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done

   if (! grep -q "Passed" <./dm_real_pexsi_c.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/dm_real_pexsi_c.log for details."
   else
      tput setaf 10
      echo " PASSED!"
      tput sgr0
   fi

   echo
   echo -n "Running the 'parallel complex density matrix solver + PEXSI' C test"
   ${MPI_EXEC} -n ${MPI_SIZE} ./test_dm_complex_c.x 3 ${ELSI_DIR}/test/H_complex.csc ${ELSI_DIR}/test/S_complex.csc 3 > dm_complex_pexsi_c.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done

   if (! grep -q "Passed" <./dm_complex_pexsi_c.log); then
      tput setaf 5
      RED_ALART="true"
      echo " FAILED!"
      tput sgr0
      echo "See `pwd`/dm_complex_pexsi_c.log for details."
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

else
   tput setaf 5
   echo "C interface has been disabled by user's choice."
   tput sgr0
fi
