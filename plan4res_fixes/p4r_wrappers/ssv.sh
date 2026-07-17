#/bin/bash

source ${INCLUDE}/sh_utils.sh

create_results_dir "${mode1}"	
if [ "$HOTSTART" = "" ]; then
	remove_previous_ssv_results "${mode1}"
fi

NBTASKS=$NBSCEN_OPT
if [ "$NBTASKS" -gt "$NB_MAX_PARALLEL_SIMUL" ] ; then
	NBTASKS=$NB_MAX_PARALLEL_SIMUL
fi
NTASKS_PER_NODE=1  # used only to compute the number of CPUs to request
N_CPUS_PER_TASK=$((CPUS_PER_NODE / NTASKS_PER_NODE ))
OMP_NUM_THREADS=$number_threads
echo "TOTAL_NB_CPUS=$TOTAL_NB_CPUS; NBTASKS=$NBTASKS; N_CPUS_PER_TASK=$N_CPUS_PER_TASK"
 
if check_param "${CONFIG}/sddp_solver.txt" "strDirOUT"; then
	replace_param "${CONFIG}/sddp_solver.txt" "strDirOUT" "${INSTANCE}/results_${mode1}/${OUT}/"
	echo -e "${print_blue}        - sddp_solver.txt config file : replaced param strDirOUT: ${INSTANCE}/results_${mode1}/${OUT}/.${no_color}" 
else
	increment_str_param_count "${CONFIG}/sddp_solver.txt"
	add_str_param "${CONFIG}/sddp_solver.txt" "strDirOUT" "${INSTANCE}/results_${mode1}/${OUT}/"
	echo -e "${print_blue}        - sddp_solver.txt config file : added param strDirOUT: ${INSTANCE}/results_${mode1}/${OUT}/.${no_color}" 
fi

if check_param "${CONFIG}/sddp_solver.txt" "intNbSimulCheckForConv"; then
	intNbSimulCheckForConv=$(get_param_value "intNbSimulCheckForConv" "${CONFIG}/sddp_solver.txt")
  
	if [[ $intNbSimulCheckForConv -gt $NBSCEN_OPT ]]; then	
		# forbids to run more simulations than scenarios
		# this should be possible but it behaves strangely, to be checked with Rafael
												   
		replace_param "${CONFIG}/sddp_solver.txt" "intNbSimulCheckForConv" "$NBSCEN_OPT"
		echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNbSimulCheckForConv by $NBSCEN_OPT.${no_color}" 
	else
		echo -e "${print_blue}        - sddp_solver.txt config file : kept value of intNbSimulCheckForConv: $intNbSimulCheckForConv ${no_color}" 
	fi
else
	increment_int_param_count "${CONFIG}/sddp_solver.txt"
	param_value=$((NBSCEN_OPT / 5))
	add_int_param "${CONFIG}/sddp_solver.txt" "intNbSimulCheckForConv" "$param_value"
	echo -e "${print_blue}        - sddp_solver.txt config file : added param intNbSimulCheckForConv: $NBSCEN_OPT.${no_color}" 
fi

if [[ ! -z "${EpsilonSSV}" ]]; then
	if check_param "${CONFIG}/sddp_solver.txt" "dblAccuracy"; then
		replace_param "${CONFIG}/sddp_solver.txt" "dblAccuracy" "$EpsilonSSV"
		echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of dblAccuracy by $EpsilonSSV.${no_color}" 
	else
		increment_dbl_param_count "${CONFIG}/sddp_solver.txt"
		add_dbl_param "${CONFIG}/sddp_solver.txt" "dblAccuracy" "$EpsilonSSV"
		echo -e "${print_blue}        - sddp_solver.txt config file : added param dblAccuracy: $EpsilonSSV.${no_color}" 
	fi
else
	EpsilonSSV=$(get_param_value "dblAccuracy" "${CONFIG}/sddp_solver.txt")
fi

# run sddp solver
echo -e "\n${print_blue}    - Run SSV [$start_time] with sddp_solver to compute Bellman values for storages: ${no_color}"	

local oldHOTSTART="$HOTSTART"
if [ ! -f "${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt" ]; then
	if [ "${HOTSTART}" = "HOTSTART" ]; then
		echo -e "\n${print_orange}    - HOTSTART requested but not possible as ${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt does not exist ${no_color}"
		oldHOTSTART="HOTSTART"
		HOTSTART="" 
	fi	
fi
											  
P4R_CMD="srun  --wckey=${WCKEY} --nodes=${SLURM_JOB_NUM_NODES} --ntasks-per-node=${NTASKS_PER_NODE} --distribution=cyclic --cpus-per-task=${N_CPUS_PER_TASK} -l --mpi=pmix "

if [[ ! $STEPS = 0 ]]; then	
	echo -e "\n${print_blue}    - Run SSV in 2 steps with max $NumberSSVIterationsFirstStep iterations, $NumberSSVForwardFirstStep Forward scenarios per iterations and convergence check each $CheckConvEachXIterInFirstStep iterations in first step  ${no_color}"	

	# first step 
	############	
	# update sddp_solver to check convergence only each $ITERCONV iterations
	if check_param "${CONFIG}/sddp_solver.txt" "intNStepConv"; then
		replace_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "$ITERCONV"
		echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNStepConv by $ITERCONV.${no_color}" 
	else
		increment_int_param_count "${CONFIG}/sddp_solver.txt"
		add_int_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "$ITERCONV"
		echo -e "${print_blue}        - sddp_solver.txt config file : added param intNStepConv: $ITERCONV.${no_color}" 
	fi
	
	# update sddp_solver to run max $STEPS iterations
	if check_param "${CONFIG}/sddp_solver.txt" "intMaxIter"; then
		intMaxIter=$(get_param_value "intMaxIter" "${CONFIG}/sddp_solver.txt")
		replace_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "$STEPS"
		echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intMaxIter by $STEPS.${no_color}" 
	else
		increment_int_param_count "${CONFIG}/sddp_solver.txt"
		add_int_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "$STEPS"
		echo -e "${print_blue}        - sddp_solver.txt config file : added param intMaxIter: $STEPS.${no_color}" 
	fi
 
 	if check_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward"; then
		replace_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "$SSVFORWARD"
		echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNbSimulForward by $SSVFORWARD.${no_color}" 
	else
		increment_int_param_count "${CONFIG}/sddp_solver.txt"
		add_int_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "$SSVFORWARD"
		echo -e "${print_blue}        - sddp_solver.txt config file : added param intNbSimulForward: $SSVFORWARD.${no_color}" 
	fi
	
	# run first step of sddp
	# hotstart run is possible in 2 steps mode
	if [ ! "$HOTSTART" = "" ]; then
		# run in hotstart
		echo -e "${print_blue}\n    - Running in HOTSTART mode using /results_${mode1}/${OUT}/cuts.txt${no_color}"
		time ${P4R_ENV} sddp_solver -l ${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/ ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4 | tee ${INSTANCE}/results_${mode1}/${OUT}/ssv_out.txt
		wait
																											  
		SSV_OUTPUT=$(grep -oP 'ACCURACY\s+\K[0-9.]+' "${INSTANCE}/results_${mode1}/${OUT}/ssv_out.txt"| tail -n 1) 

	else	
		time ${P4R_ENV} sddp_solver -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/ ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4 | tee ${INSTANCE}/results_${mode1}/${OUT}/ssv_out.txt
		wait
    SSV_OUTPUT=$(grep -oP 'ACCURACY\s+\K[0-9.]+' "${INSTANCE}/results_${mode1}/${OUT}/ssv_out.txt"| tail -n 1) 
	fi
    SSV_OUTPUT=$(printf "%06.4f" "$SSV_OUTPUT")
																		 
	echo -e "SSV_OUTPUT=$SSV_OUTPUT" 
  
	if [[ ${SSV_OUTPUT} < $EpsilonSSV ]]; then
		echo -e "${print_green}        - SSV reached optimality, no need for second step.${no_color}" 
	else
		# second step 
		############	
		if [[ ! -z "${CheckConvEachXIter}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "$CheckConvEachXIter"
			echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNStepConv by $CheckConvEachXIter.${no_color}" 
		else
			replace_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "1"
			echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNStepConv by 1.${no_color}" 
		fi
		
		if [[ ! -z "${NumberSSVIterations}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "$NumberSSVIterations"
			echo -e "${print_blue}        - updated sddp_solver.txt config file : replaced value of intMaxIter by $NumberSSVIterations.${no_color}" 
		elif [[ ! -z "${intMaxIter}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "$intMaxIter"
			echo -e "${print_blue}        - successfully updated sddp_solver.txt config file : replaced value of intMaxIter by $intMaxIter.${no_color}" 
		else
			replace_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "500"
			echo -e "${print_blue}        - successfully updated sddp_solver.txt config file : replaced value of intMaxIter by 500.${no_color}" 
		fi
      
		if [[ ! -z "${NumberSSVForward}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "$NumberSSVForward"
			echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNbSimulForward by $NumberSSVForward.${no_color}" 
		else
			replace_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "1"
			echo -e "${print_blue}        - sddp_solver.txt config file : replaced value of intNbSimulForward by 1.${no_color}" 
		fi
   
		# run second step of sddp in hotstart mode
		if [ -f "${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt" ]; then
			echo -e "    time ${P4R_ENV} sddp_solver -l ${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT} ${INSTANCE_IN_P4R}/nc4_optim/${OUT} SDDPBlock.nc4"      
			time ${P4R_ENV} sddp_solver -l ${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/ ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4
		else
			echo -e "    time ${P4R_ENV} sddp_solver -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT} ${INSTANCE_IN_P4R}/nc4_optim/${OUT} SDDPBlock.nc4"      
			time ${P4R_ENV} sddp_solver -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/ ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4
		fi
  fi
else	
	# update number max of iterations	
	if check_param "${CONFIG}/sddp_solver.txt" "intMaxIter"; then
		intMaxIter=$(get_param_value "intMaxIter" "${CONFIG}/sddp_solver.txt")
		if [[ ! -z "${NumberSSVIterations}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "$NumberSSVIterations"
			echo -e "${print_blue}        - updated sddp_solver.txt config file : replaced value of intMaxIter by $NumberSSVIterations.${no_color}" 
		fi
	elif [[ ! -z "${NumberSSVIterations}" ]]; then
		increment_int_param_count "${CONFIG}/sddp_solver.txt"		
		add_int_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "$NumberSSVIterations"
		echo -e "${print_blue}        - updated sddp_solver.txt config file : added intMaxIter: $NumberSSVIterations.${no_color}" 
	else
		increment_int_param_count "${CONFIG}/sddp_solver.txt"		
		add_int_param "${CONFIG}/sddp_solver.txt" "intMaxIter" "500"
		echo -e "${print_blue}        - updated sddp_solver.txt config file : added intMaxIter: 500.${no_color}" 
	fi
            
	if check_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward"; then
		intNbSimulForward=$(get_param_value "intNbSimulForward" "${CONFIG}/sddp_solver.txt")
		if [[ ! -z "${NumberSSVForward}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "$NumberSSVForward"
			echo -e "${print_blue}        - updated sddp_solver.txt config file : replaced value of intNbSimulForward by $NumberSSVForward.${no_color}" 
		fi
	elif [[ ! -z "${NumberSSVForward}" ]]; then
		increment_int_param_count "${CONFIG}/sddp_solver.txt"		
		add_int_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "$NumberSSVForward"
		echo -e "${print_blue}        - updated sddp_solver.txt config file : added intNbSimulForward: $NumberSSVForward.${no_color}" 
	else
		increment_int_param_count "${CONFIG}/sddp_solver.txt"		
		add_int_param "${CONFIG}/sddp_solver.txt" "intNbSimulForward" "1"
		echo -e "${print_blue}        - updated sddp_solver.txt config file : added intNbSimulForward: 1.${no_color}" 
  fi
	# update test convergence each X iteration	
	if check_param "${CONFIG}/sddp_solver.txt" "intNStepConv"; then
		intNStepConv=$(get_param_value "intNStepConv" "${CONFIG}/sddp_solver.txt")
		if [[ ! -z "${CheckConvEachXIter}" ]]; then
			replace_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "$CheckConvEachXIter"
			echo -e "${print_blue}        - updated sddp_solver.txt config file : replaced value of intNStepConv by $CheckConvEachXIter.${no_color}" 
		fi
	elif [[ ! -z "${CheckConvEachXIter}" ]]; then
		increment_int_param_count "${CONFIG}/sddp_solver.txt"		
		add_int_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "$CheckConvEachXIter"
		echo -e "${print_blue}        - updated sddp_solver.txt config file : added intNStepConv: $CheckConvEachXIter.${no_color}" 
	else
		increment_int_param_count "${CONFIG}/sddp_solver.txt"
		add_int_param "${CONFIG}/sddp_solver.txt" "intNStepConv" "1"
		echo -e "${print_blue}        - updated sddp_solver.txt config file : added intNStepConv: 1.${no_color}" 
	fi	

	if [ "$HOTSTART" = "HOTSTART" ]; then
		# run in hotstart
		echo -e "${print_blue}\n    - Running in HOTSTART mode using /results_${mode1}/${OUT}/cuts.txt ${no_color}"
		time ${P4R_ENV} sddp_solver -n ${OMP_NUM_THREADS} -l ${INSTANCE_IN_P4R}/results_${mode1}/${OUT}/cuts.txt -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/ ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4
	else   
		time ${P4R_ENV} sddp_solver -n ${OMP_NUM_THREADS} -S ${CONFIG_IN_P4R}/sddp_solver.txt -c ${CONFIG_IN_P4R}/ -p ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/ ${INSTANCE_IN_P4R}/nc4_optim/${OUT}/SDDPBlock.nc4  

	fi
fi
HOTSTART="$oldHOTSTART"
