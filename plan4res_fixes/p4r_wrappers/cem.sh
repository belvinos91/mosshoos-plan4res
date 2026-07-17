#/bin/bash

source ${INCLUDE}/sh_utils.sh

if ! check_results_dir "$mode1"; then create_results_dir "invest";  fi
remove_previous_simulation_results "invest"
if check_ssv_output "invest"; then 
	ssv_output="invest"
else
	if check_ssv_output "optim"; then 
		ssv_output="optim"
		copy_ssv_output_from_to "optim" "invest"
		echo -e "${print_orange}        - No Bellman values found in ${INSTANCE}/results_invest/$OUT/ .${no_color}" 
		echo -e "${print_orange}        - Using Bellman values from ${INSTANCE}/results_optim/$OUT/ .${no_color}" 
	else
		if check_ssv_output "simul"; then 
			ssv_output="simul"
			copy_ssv_output_from_to "simul" "invest"
			echo -e "${print_orange}        - No Bellman values found in ${INSTANCE}/results_invest/$OUT/ .${no_color}" 
			echo -e "${print_orange}        - Using Bellman values from ${INSTANCE}/results_simul/$OUT/ .${no_color}" 
		else
			return 1
		fi
	fi
fi
filter_cuts "invest"
create_results_dir "invest"
solver=$(get_solver "${CONFIG}/uc_solverconfig.txt")
solver=$(echo "$solver" | sed 's/MILPSolver//')
echo -e "\n${print_blue}     - using $solver to solve MILPs ${no_color}"


NTASKS_PER_NODE=1
N_PARAL_CEM=$number_threads
N_CPUS_PER_TASK=$number_threads
OMP_NUM_THREADS=$number_threads
if check_param "${CONFIG}/BSPar-Investment.txt" "intMaxThread"; then
	intMaxThread=$(get_param_value "intMaxThread" "${CONFIG}/BSPar-Investment.txt")
	if [[ $N_PARAL_CEM -gt $intMaxThread ]]; then
		replace_param "${CONFIG}/BSPar-Investment.txt" "intMaxThread" "$N_PARAL_CEM"
		echo -e "${print_blue}        - BSPar-Investment.txt config file : replaced value of intMaxThread: $intMaxThread by $N_PARAL_CEM.${no_color}" 
	fi
else
	increment_int_param_count "${CONFIG}/BSPar-Investment.txt"
	add_str_param "${CONFIG}/BSPar-Investment.txt" "intMaxThread" "$N_PARAL_CEM"
	echo -e "${print_blue}        - BSPar-Investment.txt config file : added param intMaxThread: $N_PARAL_CEM .${no_color}" 
fi

# run investment solver
P4R_CMD="srun --wckey=${WCKEY} --nodes=${SLURM_JOB_NUM_NODES} --ntasks-per-node=${NTASKS_PER_NODE} --cpus-per-task=${N_CPUS_PER_TASK} --mpi=pmix -l"
if [[ "$HOTSTART" == *"save"* ]]; then  
    # run in hotstart
    # select which save_state is the most recent
    if [[ "${INSTANCE_IN_P4R}/results_invest/$OUT/save_state0.nc4" -nt "${INSTANCE_IN_P4R}/results_invest/$OUT/save_state1.nc4" ]]; then
		saved="${INSTANCE_IN_P4R}/results_invest/$OUT/save_state0.nc4"
    else
		saved="${INSTANCE_IN_P4R}/results_invest/$OUT/save_state1.nc4"
    fi
    echo -e "${print_blue}        - run in HOTSTART mode (with state saving)${no_color}"
	# always use save_state0 because there are indices 0 or 1 in the nc4 and indices 1 are not read properly
	saved="${INSTANCE_IN_P4R}/results_invest/$OUT/save_state0.nc4"															 

	echo -e "\n${print_green}        - run CEM [$start_time] using investment_solver:${no_color}${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -b ${saved} -a ${INSTANCE_IN_P4R}/results_invest/$OUT/save_state -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ InvestmentBlock.nc4"
	if [ "$solver" = "HiGHS" ]; then
		# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
		pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
		time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv -b ${saved} -a ${INSTANCE_IN_P4R}/results_invest/$OUT/save_state -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4 | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt
		popd > /dev/null
	else
		# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
		pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
		time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -b ${saved} -a ${INSTANCE_IN_P4R}/results_invest/$OUT/save_state -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4 | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt		
		popd > /dev/null
	fi
elif [[ "$HOTSTART" == *"HOTSTART"* ]]; then
	echo -e "\n${print_blue}        - run in HOTSTART mode ${no_color}"
	if [ -f "${INSTANCE}/results_invest$OUT/Solution_OUT.csv" ]; then
		echo -e "\n${print_green}        - run CEM [$start_time] using investment_solver:${no_color}${P4R_ENV} investment_solver -n ${N_PARAL_CEM} -x ${INSTANCE_IN_P4R}/results_invest$OUT/Solution_OUT.csv  -l ${INSTANCE_IN_P4R}/results_invest$OUT/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/ InvestmentBlock.nc4"
		if [ "$solver" = "HiGHS" ]; then
			# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
			pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
			time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM} -x ${INSTANCE_IN_P4R}/results_invest/$OUT/Solution_OUT.csv  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4 | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt
			popd > /dev/null
		else
			# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
			pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
			time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM} -x ${INSTANCE_IN_P4R}/results_invest/$OUT/Solution_OUT.csv  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4 | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt		
			popd > /dev/null
		fi
	else
		echo -e "\n${print_orange}        - File ${INSTANCE}/results_invest$OUT/Solution_OUT.csv not found  "
		echo -e "\n${print_orange}        - run without hotstart  "
		echo -e "\n${print_green}        - run CEM [$start_time] using investment_solver:${no_color} ${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ InvestmentBlock.nc4"
		if [ "$solver" = "HiGHS" ]; then
			# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
			pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
			time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4 | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt
			popd > /dev/null
		else
			# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
			pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
			time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4 | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt		
			popd > /dev/null
		fi
	fi
else
	echo -e "\n${print_blue}        - run CEM without hotstart ${no_color}"
    if [ "$LOOPCEM" = "" ]; then	
		remove_previous_investment_results
	fi
    echo -e "\n${print_green}        - run CEM [$start_time] using investment_solver:${no_color}${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ InvestmentBlock.nc4"
    if [ "$solver" = "HiGHS" ]; then
		# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
		pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
		time ${P4R_ENV} investment_solver -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4  | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt
		popd > /dev/null
	else
		# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
		pushd "${INSTANCE_IN_P4R}/results_invest/$OUT" > /dev/null
		time ${P4R_ENV} investment_solver -a ${INSTANCE_IN_P4R}/results_invest/$OUT/save_state -n ${N_PARAL_CEM}  -l ${INSTANCE_IN_P4R}/results_invest/$OUT/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_invest/$OUT/ ${INSTANCE_IN_P4R}/nc4_invest/$OUT/InvestmentBlock.nc4  | tee ${INSTANCE}/results_invest/$OUT/invest_out.txt		
		popd > /dev/null
	fi
fi
wait
invest_output=$(cat "${INSTANCE}/results_invest/$OUT/invest_out.txt")
INVEST_OUTPUT="$invest_output"
echo -e "\n${print_green}        - CEM ended [$start_time] ${no_color}"
move_simul_results "invest"
if ! investment_status; then return 1;fi
