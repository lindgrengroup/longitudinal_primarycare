
cd /well/lindgren-ukbb/projects/ukbb-11867/georgenicholson/github_repos/longitudinal_primarycare/biobank_completion
qsub -P lindgren.prjc -q short.qc -cwd bash_script_files/bash_01_learn_local_autocorrelation.sh
qsub -P lindgren.prjc -q short.qc -cwd -hold_jid bash_01_learn_local_autocorrelation.sh bash_script_files/bash_02_fit_HMM.sh
qsub -P lindgren.prjc -q short.qc -cwd -hold_jid bash_02_fit_HMM.sh bash_script_files/bash_03_optimize_var_param.sh
qsub -P lindgren.prjc -q short.qc -cwd -hold_jid bash_03_optimize_var_param.sh bash_script_files/bash_04_genetic_assoc.sh
qsub -P lindgren.prjc -q short.qc -cwd -hold_jid bash_04_genetic_assoc.sh bash_script_files/bash_04a_begin_to_combine_outputs.sh
