#/bin/bash

source ${INCLUDE}/sh_utils.sh

#remove_previous_simulation_results "simul"
create_results_dir "simul"

if check_ssv_output "simul"; then 
    ssv_output="simul"
    echo -e "${print_grees}    - Bellman values found in ${INSTANCE}/results_simul${OUT}/ .${no_color}" 
else
    if check_ssv_output "optim"; then 
		ssv_output="optim"
		copy_ssv_output_from_to "optim" "simul"
		echo -e "${print_orange}    - No Bellman values found in ${INSTANCE}/results_simul${OUT}/ .${no_color}" 
		echo -e "${print_orange}    - Using Bellman values from ${INSTANCE}/results_optim${OUT}/ .${no_color}" 
    elif check_ssv_output "simul"; then 
		ssv_output="simul"
		copy_ssv_output_from_to "simul" "invest"
		echo -e "${print_orange}        - No Bellman values found in ${INSTANCE}/results_invest/$OUT/ .${no_color}" 
		echo -e "${print_orange}        - Using Bellman values from ${INSTANCE}/results_simul/$OUT/ .${no_color}" 
	else
       echo -e "${print_red}    - No Bellman values found in ${INSTANCE}/results_simul${OUT}/ nor $INSTANCE}/results_optim${OUT}/ .${no_color}" 
       return 1
    fi
fi

solver=$(get_solver "${CONFIG}/uc_solverconfig.txt")
solver=$(echo "$solver" | sed 's/MILPSolver//')


NTASKS_PER_NODE=1
N_PARAL_SIM=$number_threads
N_CPUS_PER_TASK=$number_threads
OMP_NUM_THREADS=$number_threads

if check_param "${CONFIG}/BSPar-Investment.txt" "intMaxThread"; then
	intMaxThread=$(get_param_value "intMaxThread" "${CONFIG}/BSPar-Investment.txt")
	if [[ $N_PARAL_SIM -gt $intMaxThread ]]; then
		replace_param "${CONFIG}/BSPar-Investment.txt" "intMaxThread" "$N_PARAL_SIM"
		echo -e "${print_blue}    - BSPar-Investment.txt config file : replaced value of intMaxThread: $intMaxThread by $N_PARAL_SIM.${no_color}" 
	fi   
else
	increment_int_param_count "${CONFIG}/BSPar-Investment.txt"
	add_str_param "${CONFIG}/BSPar-Investment.txt" "intMaxThread" "$N_PARAL_SIM"
	echo -e "${print_blue}    - BSPar-Investment.txt config file : added param intMaxThread: $N_PARAL_SIM .${no_color}" 
fi

# run investment solver
P4R_CMD="srun --wckey=${WCKEY} --nodes=1 --ntasks-per-node=${NTASKS_PER_NODE} --cpus-per-task=${N_CPUS_PER_TASK} --mpi=pmix -l"

# run investment solver
if [[ ($groupsim -eq 0) && ($onegroupsim -eq 0) ]]; then
	echo -e "\n${print_blue}  Run simulation  [$start_time] : ${no_color} ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l ${INSTANCE_IN_P4R}/results_simul/${OUT}/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_simul/$OUT/ InvestmentBlock.nc4"
	filter_cuts "simul"
	if [ "$solver" = "HiGHS" ]; then
		# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
		pushd "${INSTANCE_IN_P4R}/results_simul/${OUT}" > /dev/null
		time ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l ${INSTANCE_IN_P4R}/results_simul/${OUT}/bellmanvalues.csv -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_simul/$OUT/ ${INSTANCE_IN_P4R}/nc4_simul/$OUT/InvestmentBlock.nc4
		popd > /dev/null

	else
		# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
		pushd "${INSTANCE_IN_P4R}/results_simul/${OUT}" > /dev/null
		time ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l ${INSTANCE_IN_P4R}/results_simul/${OUT}/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_simul/$OUT/ ${INSTANCE_IN_P4R}/nc4_simul/$OUT/InvestmentBlock.nc4
		popd > /dev/null
	fi

	wait
	move_simul_results "simul"
fi

if [[ ($groupsim -eq 1) && ($sizegroupsim -ge 1) && ($onegroupsim -eq 0) ]]; then
	echo -e "\n${print_blue}     - run simulation at [$start_time] for groups of scenarios of size $sizegroupsim ${no_color}"
	update_yaml_param "${CONFIG}/settingsCreateInputPlan4res.yml" 2 "ParametersCreate invest" no
	
	# run groups of simulations in parallel
	
	# compute number of groups
	BASE=$(( NBSCEN_SIM / sizegroupsim ))
	RESTE=$(( NBSCEN_SIM % sizegroupsim ))
	if [ $RESTE -ge 1 ]; then
		nb_group_sim=$(( BASE + 1 ))
	else	
		nb_group_sim=$BASE
	fi
	
	current_first_scen=0
	sizecurrentgroup=$sizegroupsim
	number_remaining_scen=$NBSCEN_SIM
	for ((grp=0; grp<$nb_group_sim; grp++)); do
		# format for current group of scenarios
		if [ $number_remaining_scen -lt $sizegroupsim ] ; then
			sizecurrentgroup=$number_remaining_scen
		else
			sizecurrentgroup=$sizegroupsim
		fi
		echo -e "\n${print_blue}     - run group of scenario $grp of size $sizecurrentgroup starting at scenario $current_first_scen ${no_color}"

		if [ $sizecurrentgroup -gt 0 ]; then
			
			create_format_settings_group ${CONFIG_IN_P4R}/settings_format_${mode1}.yml ${INSTANCE_IN_P4R}/results_simul/${OUT}/settings_format_${mode1}_$grp.yml $grp $current_first_scen $sizecurrentgroup
					
			# create netcdf
			if [[ ! -d "${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}" ]]; then
				mkdir "${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}"
			fi
			if [[ ! -d "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}" ]]; then
				mkdir "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}"
			fi
			if [ "$FORMAT" = "FORMATGROUP" ]; then
				echo -e "\n${print_blue}        - Create netcdf input files: ${no_color}${P4R_ENV} python -W ignore ${PYTHONSCRIPTS_IN_P4R}/format.py ${INSTANCE_IN_P4R}/results_simul/${OUT}/settings_format_${mode1}_$grp.yml ${CONFIG_IN_P4R}/settingsCreateInputPlan4res.yml ${DATASET}"

				cp ${INSTANCE}/results_${mode1}/${OUT}/bellmanvalues.csv ${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}/
				P4R_CMD="srun --wckey=${WCKEY} --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 -J Format --mpi=pmix -l"
				${P4R_ENV} python -W ignore ${PYTHONSCRIPTS_IN_P4R}/format.py "${INSTANCE_IN_P4R}/results_simul/${OUT}/settings_format_${mode1}_${grp}.yml" "${CONFIG_IN_P4R}/settingsCreateInputPlan4res.yml" "${DATASET}" ${number_threads}
				wait
			fi
			filter_cuts "simul" "$grp"

			# run simulation
			if [[ -d "${INSTANCE_IN_P4R}/nc4_simul_${grp}" ]]; then 
				P4R_CMD="srun --wckey=${WCKEY} --nodes=1 --ntasks-per-node=${NTASKS_PER_NODE} --cpus-per-task=${N_CPUS_PER_TASK} --mpi=pmix -l"			
				echo -e "\n${print_blue}  Run simulation for group $grp [$start_time] : ${no_color} ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l ${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}/ InvestmentBlock.nc4"
				if [ "$solver" = "HiGHS" ]; then
					# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
					pushd "${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}" > /dev/null
					time ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l "${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}/bellmanvalues.csv" -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}/" ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}/InvestmentBlock.nc4
					popd > /dev/null

				else
					# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
					pushd "${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}" > /dev/null
					time ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l "${INSTANCE_IN_P4R}/results_simul/${OUT}_${grp}/bellmanvalues.csv"  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}/" ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp}/InvestmentBlock.nc4
					popd > /dev/null
				fi
				wait
				
				# renumber scenarios
				move_simul_results_group $grp $current_first_scen $sizecurrentgroup
			else
				echo -e "\n${print_red}        - ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${grp} does not exist, run with option -F ${no_color}"			
			fi
			old_current_first_scen=$current_first_scen
			current_first_scen=$(( old_current_first_scen + sizecurrentgroup ))
			number_remaining_scen=$(( number_remaining_scen - sizecurrentgroup ))
			
		fi
	done
	move_simul_results "simul"
fi

if [[  ($onegroupsim -eq 1) ]]; then
	echo -e "\n${print_blue}     - run simulation at [$start_time] for the group $numgroupsim of scenarios of size $sizegroupsim ${no_color}"
	update_yaml_param "${CONFIG}/settingsCreateInputPlan4res.yml" 2 "ParametersCreate invest" no
	
	# run 1 groups of simulations 
	
	# compute number of groups
	BASE=$(( NBSCEN_SIM / sizegroupsim ))
	RESTE=$(( NBSCEN_SIM % sizegroupsim ))
	if [ $RESTE -ge 1 ]; then
		nb_group_sim=$(( BASE + 1 ))
	else	
		nb_group_sim=$BASE
	fi
	
	current_first_scen=0
	sizecurrentgroup=$sizegroupsim
	number_remaining_scen=$NBSCEN_SIM
	for ((grp=0; grp<$numgroupsim; grp++)); do
		# format for current group of scenarios
		if [ $number_remaining_scen -lt $sizegroupsim ] ; then
			sizecurrentgroup=$number_remaining_scen
		else
			sizecurrentgroup=$sizegroupsim
		fi
		old_current_first_scen=$current_first_scen
		current_first_scen=$(( old_current_first_scen + sizecurrentgroup ))
		number_remaining_scen=$(( number_remaining_scen - sizecurrentgroup ))
	done
	if [ $sizecurrentgroup -gt 0 ]; then
		
		create_format_settings_group ${CONFIG_IN_P4R}/settings_format_${mode1}.yml ${INSTANCE_IN_P4R}/results_simul/${OUT}/settings_format_${mode1}_$numgroupsim.yml $numgroupsim $current_first_scen $sizecurrentgroup
				
		# create netcdf
		if [[ ! -d "${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}" ]]; then
			mkdir "${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}"
		fi
		if [[ ! -d "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}" ]]; then
			mkdir "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}"
		fi
		
		if [ "$FORMAT" = "FORMATGROUP" ]; then
			echo -e "\n${print_blue}        - Create netcdf input files for group ${numgroupsim} : ${no_color}${P4R_ENV} python -W ignore ${PYTHONSCRIPTS_IN_P4R}/format.py ${INSTANCE_IN_P4R}/results_simul/${OUT}/settings_format_${mode1}_${numgroupsim}.yml ${CONFIG_IN_P4R}/settingsCreateInputPlan4res.yml ${DATASET} ${number_threads}"

			P4R_CMD="srun --wckey=${WCKEY} --nodes=1 --ntasks=1 --ntasks-per-node=1 --cpus-per-task=1 -J Format --mpi=pmix -l"
			cp ${INSTANCE}/results_${mode1}/${OUT}/bellmanvalues.csv ${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}/
			${P4R_ENV} python -W ignore ${PYTHONSCRIPTS_IN_P4R}/format.py "${INSTANCE_IN_P4R}/results_simul/${OUT}/settings_format_${mode1}_${numgroupsim}.yml" "${CONFIG_IN_P4R}/settingsCreateInputPlan4res.yml" "${DATASET}" ${number_threads}
			wait
		fi
		filter_cuts "simul" "$numgroupsim"

		# run simulation 
		if [[ -d "${INSTANCE_IN_P4R}/nc4_simul_${numgroupsim}" ]]; then 
			P4R_CMD="srun --wckey=${WCKEY} --nodes=1 --ntasks-per-node=${NTASKS_PER_NODE} --cpus-per-task=${N_CPUS_PER_TASK} --mpi=pmix -l"
			echo -e "\n${print_blue}  Run simulation for group $numgroupsim [$start_time] : ${no_color} ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l ${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}/bellmanvalues.csv  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}/ InvestmentBlock.nc4"
			if [ "$solver" = "HiGHS" ]; then
				# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
				pushd "${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}" > /dev/null
				time ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l "${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}/bellmanvalues.csv" -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}/" ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}/InvestmentBlock.nc4
				popd > /dev/null

			else
				# investment_solver writes its *_OUT.csv to CWD (no -d in this build)
				pushd "${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}" > /dev/null
				time ${P4R_ENV} investment_solver -n ${N_PARAL_SIM} -s  -l "${INSTANCE_IN_P4R}/results_simul/${OUT}_${numgroupsim}/bellmanvalues.csv"  -o -S ${CONFIG_IN_P4R}/BSPar-Investment.txt -c ${CONFIG_IN_P4R}/ -p "${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}/" ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim}/InvestmentBlock.nc4
				popd > /dev/null
			fi
			wait
			
			# renumber scenarios
			move_simul_results_group $numgroupsim $current_first_scen $sizecurrentgroup		
			move_simul_results "simul"
		else
			echo -e "\n${print_red}        - ${INSTANCE_IN_P4R}/nc4_simul/$OUT_${numgroupsim} does not exist, run with option -F ${no_color}"			
		fi
	fi
	
	
fi


wait

echo -e "\n${print_green}- simulations completed [$start_time] ${no_color}"

