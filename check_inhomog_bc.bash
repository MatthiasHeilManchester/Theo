#! /bin/bash



# Setup directories
main_dir=NEW_RUNS_CHECK_INHOMOGENEOUS_BCS
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


# Loop over the element areas
area_list="0.1 0.05 0.025 0.0125 0.006125"
for area in `echo $area_list`; do
    

    # Loop over scaling factor for domain
    constant_list="0.0 1.0" # 0.1 0.01"
    for constant in `echo $constant_list`; do
        
        
        run_dir=CASE_area`echo $area`_constant`echo $constant`
        mkdir $run_dir
        cp $code $src $run_dir
        cd $run_dir
        mkdir RESLT
        output_file=OUTPUT
        
        ./$code --uniform_element_area $area --constant_in_singular_solution $constant  --use_dirichlet_bcs --enforce_dirichlet_bcs_by_lagrange_multipliers > $output_file

#       

        # Post-process
        cd RESLT
        oomph-convert extended_soln*.dat; makePvd extended_soln extended_soln.pvd
        
        echo " "
        echo "To post-process:" 
        echo " " 
        echo "  paraview extended_soln.pvd "
        echo " " 
        cd ..
        
        # Out of the case directory
        cd ..

    done # constant

done # area

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
