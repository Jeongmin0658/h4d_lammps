# --- Move the ouput files to "output" folder ---
#
mv mc_statistics.out restart_after_MC_*.data logfile_hneqMDMC.* spce_lvmd_bulk.log log.lammps restart_end_run.data output

rm -rf __pycache__
