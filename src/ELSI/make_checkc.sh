#!/bin/bash
# VY: This is a simple bash script used for "make checkc"

if [ "$C_INTERFACE" = "yes" ]
then

rm ev_real_serial_c.log 2> /dev/null
rm ev_real_elpa_c.log 2> /dev/null
rm dm_real_elpa_c.log 2> /dev/null
rm dm_real_libomm_c.log libOMM.log 2> /dev/null
rm dm_real_pexsi_c.log logPEXSI0 2> /dev/null
set -e # Stop on error

RED_ALART="false"
echo
echo "Test program output may be found in $PWD"

echo
echo -n "Running the serial 'elsi_ev_real' C test"
${MPI_EXEC} -n 1 ./test_ev_real_c.x ${TOMATO_SEED} 1 > ev_real_serial_c.log &
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
echo -n "Running the parallel 'elsi_ev_real + ELPA' C test"
${MPI_EXEC} -n 4 ./test_ev_real_c.x ${TOMATO_SEED} 1 > ev_real_elpa_c.log &
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
echo -n "Running the parallel 'elsi_dm_real + ELPA' C test"
${MPI_EXEC} -n 4 ./test_dm_real_c.x ${TOMATO_SEED} 1 > dm_real_elpa_c.log &
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
echo -n "Running the parallel 'elsi_dm_real + libOMM' C test"
${MPI_EXEC} -n 4 ./test_dm_real_c.x ${TOMATO_SEED} 2 > dm_real_libomm_c.log &
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
if [ "$DISABLE_CXX" != "yes" ]
then
   echo -n "Running the parallel 'elsi_dm_real + PEXSI' C test"
   ${MPI_EXEC} -n 4 ./test_dm_real_c.x ${TOMATO_SEED} 3 > dm_real_pexsi_c.log &
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
else
   tput setaf 5
   echo "PEXSI has been disabled because ELSI is compiled without C++ support."
   echo "Please recompile ELSI with C++ support to enable PEXSI, or continue"
   echo "using ELSI without PEXSI."
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
   echo "C interface has been disabled by user. Please recompile ELSI with"
   echo "C interface, or continue using ELSI in Fortran."
   tput sgr0
fi
