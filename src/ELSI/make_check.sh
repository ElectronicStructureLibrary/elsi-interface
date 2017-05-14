#!/bin/bash
# VY: This is a simple bash script used for "make check"
# More tests can be done by using test_ev_real.x and test_dm_real.x

rm ev_real_elpa.log 2> /dev/null
rm dm_real_elpa.log 2> /dev/null
rm dm_real_libomm.log libOMM.log 2> /dev/null
rm dm_real_pexsi.log logPEXSI0 2> /dev/null
set -e # Stop on error

RED_ALART="false"
echo
echo "Test program output may be found in $PWD"

echo
echo -n "Running the 'elsi_ev_real + ELPA' test"
${MPI_EXEC} -n 4 ./test_ev_real.x ${TOMATO_SEED} 1 > ev_real_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

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
echo -n "Running the 'elsi_dm_real + ELPA' test"
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 1 > dm_real_elpa.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
    sleep 1
    echo -n '.'
done

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
echo -n "Running the 'elsi_dm_real + libOMM' test"
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 2 > dm_real_libomm.log &
PID=$!
while kill -0 $PID 2>/dev/null; do
   sleep 1
   echo -n '.'
done

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
if [ "$DISABLE_CXX" != "yes" ]
then
   echo -n "Running the 'elsi_dm_real + PEXSI' test"
   ${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 3 > dm_real_pexsi.log &
   PID=$!
   while kill -0 $PID 2>/dev/null; do
      sleep 1
      echo -n '.'
   done

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
else
   tput setaf 5
   echo "PEXSI has been disabled because ELSI is compiled without C++ support."
   echo "Please recompile ELSI with C++ support to enable PEXSI, or continue"
   echo "using ELSI without PEXSI."
   tput sgr0
fi

echo
# WHY AM I YELLING?
if [ "$RED_ALART" = "true" ] 
then
   tput setaf 5
   echo "MAKE CHECK FAILED, CHECK YOUR COMPILATION SETTINGS!"
   tput sgr0
else
   tput setaf 10
   echo "MAKE CHECK PASSED, YOU'RE GOOD TO GO!"
   tput sgr0
fi

echo
