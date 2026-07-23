# RUNBOOK - how to get results on this server

## Directory map
~/plan4res/p4r-env    = plan4res engine (built, HiGHS fixes compiled in). Never edit.
~/plan4res_work/data  = workshop: datasets + results. Run p4r HERE.
~/mosshoos-plan4res   = this git repo (the public record). Commit results HERE.

## A. plan4res (the fixed framework, on HiGHS - no licence)
cd ~/plan4res_work/data
p4r SSV toyDataset            # water values -> results_simul/cuts.txt
p4r SIM toyDataset            # stochastic simulation
p4r POSTTREAT toyDataset      # tables + IMG/ graphs
# or combined: p4r SSVandSIM toyDataset
# outputs: toyDataset/results_simul/  (ActivePower/, IMG/, MarginalCosts/, Volume/...)

## B. MOSSHOOS (our model)
cd ~/mosshoos-plan4res
/usr/bin/python3.12 src/mosshoos_run_all.py --inflow data/OMVS_inflows.csv --n_iter 200 --n_mc 500
/usr/bin/python3.12 storage_sweep.py --inflow data/OMVS_inflows.csv
/usr/bin/python3.12 wavelet_ablation.py --n_iter 200 --n_eval 300
/usr/bin/python3.12 compare_water_values_v2.py \
    --dataset ~/plan4res_work/data/toyDataset \
    --p4r-cuts ~/plan4res_work/data/toyDataset/results_simul/cuts.txt
/usr/bin/python3.12 residual_demand_check.py
# outputs land in results/ (inside the repo -> commit them)

## C. Dashboards
/usr/bin/python3.12 dashboard_server.py --port 8050    # interactive (run model, download data)
# static, no server: docs/index.html  /  Pages: https://belvinos91.github.io/mosshoos-plan4res/

## D. Publish new results
git add results/ && git commit -m "describe the run" && git push origin main

## Note on python
Use /usr/bin/python3.12 explicitly (the default /usr/local python lacks _ctypes).
