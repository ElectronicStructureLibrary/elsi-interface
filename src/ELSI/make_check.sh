#!/bin/bash -l
# VY: This is a simple bash script used for "make check"
# More tests can be done by using test_ev_real.x and test_dm_real.x

set -e # Stop on error

echo
${MPI_EXEC} -n 4 ./test_ev_real.x ${TOMATO_SEED} 1 > ev_real_elpa.log
if (! grep -q "Passed" <./ev_real_elpa.log); then
   tput setaf 5
   echo "FAILED:  elsi_ev_real + ELPA (mp)"
   tput sgr0
   echo "See `pwd`/ev_real_elpa.log for details."
else
   tput setaf 10
   echo "PASSED:  elsi_ev_real + ELPA (mp)"
   tput sgr0
   rm ev_real_elpa.log
fi

if [ "$ELPA2_KERNEL" != "GPU" ]
then
   echo
   ${MPI_EXEC} -n 1 ./test_ev_real.x ${TOMATO_SEED} 1 > ev_real_elpa_sp.log
   if (! grep -q "Passed" <./ev_real_elpa_sp.log); then
      tput setaf 5
      echo "FAILED:  elsi_ev_real + ELPA (sp)"
      tput sgr0
      echo "See `pwd`/ev_real_elpa_sp.log for details."
   else
      tput setaf 10
      echo "PASSED:  elsi_ev_real + ELPA (sp)"
      tput sgr0
      rm ev_real_elpa_sp.log
   fi
fi

echo
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 1 > dm_real_elpa.log
if (! grep -q "Passed" <./dm_real_elpa.log); then
   tput setaf 5
   echo "FAILED:  elsi_dm_real + ELPA"
   tput sgr0
   echo "See `pwd`/dm_real_elpa.log for details."
else
   tput setaf 10
   echo "PASSED:  elsi_dm_real + ELPA"
   tput sgr0
   rm dm_real_elpa.log
fi

echo
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 2 > dm_real_libomm.log
if (! grep -q "Passed" <./dm_real_libomm.log); then
   tput setaf 5
   echo "FAILED:  elsi_dm_real + libOMM"
   tput sgr0
   echo "See `pwd`/dm_real_libomm.log for details."
else
   tput setaf 10
   echo "PASSED:  elsi_dm_real + libOMM"
   tput sgr0
   rm dm_real_libomm.log libOMM.log
fi

echo
if [ "$DISABLE_CXX" != "yes" ]
then
   ${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 3 > dm_real_pexsi.log
   if (! grep -q "Passed" <./dm_real_pexsi.log); then
      tput setaf 5
      echo "FAILED:  elsi_dm_real + PEXSI"
      tput sgr0
      echo "See `pwd`/dm_real_pexsi.log for details."
   else
      tput setaf 10
      echo "PASSED:  elsi_dm_real + PEXSI"
      tput sgr0
      rm dm_real_pexsi.log logPEXSI0
   fi
else
   tput setaf 5
   echo "PEXSI has been disabled because ELSI is compiled without C++ support."
   echo "Please recompile ELSI with C++ support to enable PEXSI, or continue"
   echo "using ELSI without PEXSI."
   tput sgr0
fi

echo
