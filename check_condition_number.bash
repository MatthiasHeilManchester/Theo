#! /bin/bash



# Setup directories
main_dir=NEW_RUNS_CHECK_CONDITION_NUMBER
if [ -e $main_dir ]; then
    echo "Directory " $main_dir "already exists -- please move it out of the way."
    exit
fi
mkdir $main_dir


# Build and copy
code=two_d_poisson_with_face_element
src=two_d_poisson_with_face_element.cc
make $code
cp $code $src $main_dir
cd $main_dir

# Problem type
problem_type_list="real_problem fe_only fake_rc"
for problem_type in `echo $problem_type_list`; do
    
    
    # Default: real problem
    other_flags=" --check_condition_number 0 "
    if [ "$problem_type" == "fe_only" ]; then
        other_flags=" --check_condition_number 1 --dont_use_singularity "
    fi
    
    if [ "$problem_type" == "fake_rc" ]; then
        other_flags=" --check_condition_number 2 "
    fi
    
    
    
    # Dirichlet or Neumann conditions?
    dirichlet_neumann_list="dirichlet neumann"
    for dirichlet_neumann in `echo $dirichlet_neumann_list`; do
        
        dirichlet_neumann_flag=" "
        pin_method_list="not_dirichlet"
        dirichlet_neumann_dir_label="neumann"
        if [ "$dirichlet_neumann" == "dirichlet" ]; then
            pin_method_list="dirichlet_via_lagrange_multipliers"
            if [ "$problem_type" == "fe_only" ]; then
                pin_method_list="dirichlet_via_pin dirichlet_via_lagrange_multipliers"
            fi
        fi
        
        # How do we pin (if applicable!)
        for pin_method in `echo $pin_method_list`; do
            
            if [ "$pin_method" == "dirichlet_via_lagrange_multipliers" ]; then 
                dirichlet_neumann_dir_label="dirichlet_via_lagrange_multipliers"   
                dirichlet_neumann_flag=" --use_dirichlet_bcs --enforce_dirichlet_bcs_by_lagrange_multipliers"
            fi
            if [ "$pin_method" == "dirichlet_via_pin" ]; then 
                dirichlet_neumann_dir_label="dirichlet_via_pin"   
                dirichlet_neumann_flag=" --use_dirichlet_bcs "
            fi
            
            
            # Scaling factor for overall size of domain
            scaling_factor_list="1.0" # hierher 0.1 0.01"
            for scaling_factor in `echo $scaling_factor_list`; do
                
                
                case_without_area=CASE_`echo $problem_type`_scaling_factor`echo $scaling_factor`"_"`echo $dirichlet_neumann_dir_label`
                echo "ZONE T=\"$case_without_area\"" >> condition_number.dat
                
                # Loop over the element areas
                area_list="0.1 0.05 0.025 0.0125 0.006125 0.0030625 0.00153125 0.0005"
                for area in `echo $area_list`; do
                    
                    
                    run_dir=CASE_`echo $problem_type`_area`echo $area`_scaling_factor`echo $scaling_factor`"_"`echo $dirichlet_neumann_dir_label`
                    mkdir $run_dir
                    cp $code $src $run_dir
                    cd $run_dir
                    mkdir RESLT
                    output_file=OUTPUT
                    
                    echo "Doing: "$run_dir
                    echo "with flags: "
                    echo " " 
                    flags="$other_flags  --uniform_element_area $area --scaling_factor_for_domain $scaling_factor $dirichlet_neumann_flag"
                    echo "  "$flags
                    echo " " 
                    
                    ./$code $flags  > $output_file
                    
                    
                    # Post-process
                    #cd RESLT
                    #oomph-convert -p2 integrand_elements*.dat
                    #cd ..
                    
                    grep 'Matrix size and condition number'  $output_file | awk '{print $6 " " $7}' >> ../condition_number.dat
                    
                    # Out of the case directory
                    cd ..
                    
                    
                done # area
                
            done # scaling
            
        done # problem type
    
    done # pin method (if applicable)

done # dirichlet neumann
    






echo " "
echo "To post-process:" 
echo " " 
echo "   tecplot NEW_RUNS_CHECK_CONDITION_NUMBER/condition_number.lay " 
echo " " 


exit


#----------------
mesher_postfix_list="_gmsh" # "_tetgen  _gmsh" # _tetgen"
code=""
for mesher_postfix in `echo $mesher_postfix_list`; do
    code=assess_tet_mesh_volume$mesher_postfix
    make $code
    cp $code $main_dir
    cd $main_dir
                
    run_dir=CASE`echo $mesher_postfix`
    mkdir $run_dir
    cp $code $run_dir
    cd $run_dir
    mkdir RESLT
    output_file=OUTPUT
    
    nref=10 # 4 # 5
    nunref=10 # 4 # 5

    target_volume=0.1
    if [ $mesher_postfix == "_gmsh" ]; then
        target_volume=1.0
    fi

    ./$code --suppress_inner_boundaries --desired_target_volume $target_volume --check_uniform_refinement $nref --check_uniform_unrefinement $nunref --specify_gmsh_target_volume_via_file --refinement_volume_decrease_factor 1.2  > $output_file


# --unrefinement_volume_increase_factor 3.0

#  --refinement_volume_decrease_factor 1.2  
    cd ..
    cd ..
done
