#!/bin/bash -l
# VY: This is a simple bash script used for "make check"
# More tests can be done by using test_ev_real.x and test_dm_real.x

set -e # Stop on error

echo
${MPI_EXEC} -n 4 ./test_ev_real.x ${TOMATO_SEED} 1 > ev_real_elpa.log
if (! grep -q "Passed" <./ev_real_elpa.log); then
   echo "FAILED:  elsi_ev_real + ELPA (mp)"
   echo "See `pwd`/ev_real_elpa.log for details."
else
   echo "PASSED:  elsi_ev_real + ELPA (mp)"
   rm ev_real_elpa.log
fi

echo
${MPI_EXEC} -n 1 ./test_ev_real.x ${TOMATO_SEED} 1 > ev_real_elpa_sp.log
if (! grep -q "Passed" <./ev_real_elpa_sp.log); then
   echo "FAILED:  elsi_ev_real + ELPA (sp)"
   echo "See `pwd`/ev_real_elpa_sp.log for details."
   echo "- If you are using the GPU-accelerated version of ELPA and the log file contains the error"
   echo "- 'ELPA GPU version needs blocksize = 128', you may safely ignore the failure of this test."
else
   echo "PASSED:  elsi_ev_real + ELPA (sp)"
   rm ev_real_elpa_sp.log
fi

echo
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 1 > dm_real_elpa.log
if (! grep -q "Passed" <./dm_real_elpa.log); then
   echo "FAILED:  elsi_dm_real + ELPA"
   echo "See `pwd`/dm_real_elpa.log for details."
else
   echo "PASSED:  elsi_dm_real + ELPA"
   rm dm_real_elpa.log
fi

echo
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 2 > dm_real_libomm.log
if (! grep -q "Passed" <./dm_real_libomm.log); then
   echo "FAILED:  elsi_dm_real + libOMM"
   echo "See `pwd`/dm_real_libomm.log for details."
else
   echo "PASSED:  elsi_dm_real + libOMM"
   rm dm_real_libomm.log libOMM.log
fi

echo
${MPI_EXEC} -n 4 ./test_dm_real.x ${TOMATO_SEED} 3 > dm_real_pexsi.log
if (! grep -q "Passed" <./dm_real_pexsi.log); then
   echo "FAILED:  elsi_dm_real + PEXSI"
   echo "See `pwd`/dm_real_pexsi.log for details."
   echo "- If you compiled ELSI without C++ support, PEXSI is disabled, and this failure is expected."
   echo "- In this case, please recompile ELSI with C++ support to enable PEXSI, or ignore this error"
   echo "- message to continue using ELSI without PEXSI."
else
   echo "PASSED:  elsi_dm_real + PEXSI"
   rm dm_real_pexsi.log logPEXSI0
fi

echo
