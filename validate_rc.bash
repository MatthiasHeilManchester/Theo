#! /bin/bash



# Setup directories
main_dir=NEW_RUNS_VALIDATE_RC
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
area_list="0.1 0.05 0.025 0.0125 0.006125 0.0030625 0.00153125"
for area in `echo $area_list`; do
    
        
    run_dir=CASE_area`echo $area`
    mkdir $run_dir
    cp $code $src $run_dir
    cd $run_dir
    mkdir RESLT
    output_file=OUTPUT
    
    ./$code --check_rc_equation --uniform_element_area $area > $output_file
       
    # Post-process
    cd RESLT
    oomph-convert -p2 integrand_elements*.dat
    cd ..

    # Out of the case directory
    cd ..
    
done # area

    grep 'Residual if we assign the exact non-singular part of the solution for area' CASE*/$output_file | awk '{print $(NF-7) " " $(NF-2) " " $NF}' > rc_vs_refinement_trace_for_exact_non_singular_part.dat

    grep 'Residual if we assign the exact non-singular part of the solution to nodal values' CASE*/$output_file | awk '{print $(NF-7) " " $(NF-2) " " $NF}' > rc_vs_refinement_trace_for_nodes_are_set_to_non_singular_part.dat

    ~/bin/make_combined_pvd.bash CASE*/RESLT/integrand_elements0.vtp
    mv combined.pvd exact_integrand.pvd
    ~/bin/make_combined_pvd.bash CASE*/RESLT/integrand_elements1.vtp
    mv combined.pvd fe_integrand.pvd
    
    echo " "
    echo "To post-process:" 
    echo " " 
    echo "  Load exact_integrand.pvd and fe_integrand.pvd into paraview"
    echo "  and warp by last-but-one ****check order!*** variable (integrand!) visualise with"
    echo "  glyphs."
    echo " "
    echo "  tecplot rc_vs_refinement_trace_for_exact_non_singular_part.dat rc_vs_refinement_trace_for_nodes_are_set_to_non_singular_part.dat"
    echo " "
    echo "  shows r_c with exact and FE-computed derivatives as fct of mesh"
    echo "  refinement" 
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
