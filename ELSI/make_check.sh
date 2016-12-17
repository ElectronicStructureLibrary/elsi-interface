#!/bin/bash -l
set -e # Stop on error

mpirun -n 4 ./test_ev_real.x $1 1 > ev_real_elpa.log
if (! grep -q "Passed" <./ev_real_elpa.log); then
   echo "FAILED:  elsi_ev_real + ELPA (mp)"
   echo "See `pwd`/ev_real_elpa.log for details."
else
   echo "PASSED:  elsi_ev_real + ELPA (mp)"
   rm ev_real_elpa.log
fi

mpirun -n 1 ./test_ev_real.x $1 1 > ev_real_elpa_sp.log
if (! grep -q "Passed" <./ev_real_elpa_sp.log); then
   echo "FAILED:  elsi_ev_real + ELPA (sp)"
   echo "See `pwd`/ev_real_elpa_sp.log for details."
else
   echo "PASSED:  elsi_ev_real + ELPA (sp)"
   rm ev_real_elpa_sp.log
fi

mpirun -n 4 ./test_dm_real.x $1 1 > dm_real_elpa.log
if (! grep -q "Passed" <./dm_real_elpa.log); then
   echo "FAILED:  elsi_dm_real + ELPA"
   echo "See `pwd`/dm_real_elpa.log for details."
else
   echo "PASSED:  elsi_dm_real + ELPA"
   rm dm_real_elpa.log
fi

mpirun -n 4 ./test_dm_real.x $1 2 > dm_real_libomm.log
if (! grep -q "Passed" <./dm_real_libomm.log); then
   echo "FAILED:  elsi_dm_real + libOMM"
   echo "See `pwd`/dm_real_libomm.log for details."
else
   echo "PASSED:  elsi_dm_real + libOMM"
   rm dm_real_libomm.log libOMM.log
fi

mpirun -n 4 ./test_dm_real.x $1 3 > dm_real_pexsi.log
if (! grep -q "Passed" <./dm_real_pexsi.log); then
   echo "FAILED:  elsi_dm_real + PEXSI"
   echo "See `pwd`/dm_real_pexsi.log for details."
else
   echo "PASSED:  elsi_dm_real + PEXSI"
   rm dm_real_pexsi.log logPEXSI0
fi
