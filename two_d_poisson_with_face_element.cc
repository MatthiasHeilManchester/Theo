//LIC// ====================================================================
//LIC// This file forms part of oomph-lib, the object-oriented,
//LIC// multi-physics finite-element library, available
//LIC// at http://www.oomph-lib.org.
//LIC//
//LIC//    Version 1.0; svn revision $LastChangedRevision: 1097 $
//LIC//
//LIC// $LastChangedDate: 2015-12-17 11:53:17 +0000 (Thu, 17 Dec 2015) $
//LIC//
//LIC// Copyright (C) 2006-2016 Matthias Heil and Andrew Hazel
//LIC//
//LIC// This library is free software; you can redistribute it and/or
//LIC// modify it under the terms of the GNU Lesser General Public
//LIC// License as published by the Free Software Foundation; either
//LIC// version 2.1 of the License, or (at your option) any later version.
//LIC//
//LIC// This library is distributed in the hope that it will be useful,
//LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
//LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//LIC// Lesser General Public License for more details.
//LIC//
//LIC// You should have received a copy of the GNU Lesser General Public
//LIC// License along with this library; if not, write to the Free Software
//LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
//LIC// 02110-1301  USA.
//LIC//
//LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
//LIC//
//LIC//====================================================================
//Driver for Poisson in backward step domain -- meshed with triangle

//Generic includes
#include "generic.h"
#include "poisson.h"

// The mesh
#include "meshes/triangle_mesh.h"


#include "new_poisson_sing_face_element.h"

#include "power_method.h"

using namespace std;

using namespace oomph;

//==start_of_namespace==============================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{

 /// Number that identifies which
 /// case we're doing when checking the condition number
 /// 0: real problem; 1: FE-only; 2: fake r_c equation
 unsigned Problem_type_for_check_condition_number=0;

 // Default uniform element area
 double Uniform_element_area=0.1;

 // If true: Flux bc on boundaries 0-3; pin on boundaries 4-6
 bool Do_flux_problem=true;

 // Default rescaling factor for domain
 double Scaling_factor_for_domain = 1.0;
 
 // Dimensionless domain values
 // ---------------------------

 /// Dimless width of inflow channel
 double H_up = 1.0; 

 /// Dimless length of channel upstream of step
 double L_up = 2.0; 

 /// \short Dimless length of channel downstream of step
 double L_down = 2.0; 

 /// Dimless width of outflow channel
 double H_down = 2.0; 

 /// Radius of internal boundary (surrounding the singularity)
 double Radius_of_internal_boundary=0.5;

 // hierher kill
 // /// Test for assigning C directly
 // double Imposed_amplitude = 0.0;

// Bit of a hack but it facilitates reuse...
#include "unstructured_backward_step_mesh.h"


 /// Non-singular part of the solution for testing
 double u_non_singular(const Vector<double>& x);

 /// Non-singular part of the solution for testing (as vector)
 void u_non_singular_as_vector(const Vector<double>& x,
                               Vector<double>& u)
 {
  u.resize(1);
  u[0]=u_non_singular(x);
 }

 /// \short "Singular" function and gradient
 void singular_fct_and_gradient(const Vector<double>& x,
                                double& u, Vector<double>& du_dx)
 {
  // Radius & polar angle
  double r=sqrt((x[0]-L_up)*(x[0]-L_up)+x[1]*x[1]);

  // Little hack to make sure atan2 doesn't overshoot
  double y=x[1];
  double tol_y=-1.0e-12;
  if ((y<0.0)&&(y>tol_y)) y=0.0;
  double phi=atan2(y,(x[0]-L_up));

  u =                       pow(r, 2.0/3.0)*sin((2.0/3.0)*(phi+MathematicalConstants::Pi/2.0));
  double dudr =   (2.0/3.0)*pow(r,-1.0/3.0)*sin((2.0/3.0)*(phi+MathematicalConstants::Pi/2.0));
  double dudphi = (2.0/3.0)*pow(r, 2.0/3.0)*cos((2.0/3.0)*(phi+MathematicalConstants::Pi/2.0));

  du_dx[0]=dudr*cos(phi)-1.0/r*dudphi*sin(phi);
  du_dx[1]=dudr*sin(phi)+1.0/r*dudphi*cos(phi);
 }


 /// \short "Singular" function
 double singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(2);
  singular_fct_and_gradient(x,u,du_dx);
  return u;
 }

 /// \short Gradient of "Singular" function
 Vector<double> gradient_of_singular_fct(const Vector<double>& x)
 {
  double u=0.0;
  Vector<double> du_dx(2);
  singular_fct_and_gradient(x,u,du_dx);
  return du_dx;
 }

 /// Exact solution
 double u_exact(const Vector<double>& x)
 {
  double u=0.0;
  if (!CommandLineArgs::command_line_flag_has_been_set
       ("--suppress_sing_in_exact_soln"))
   {
    u+=singular_fct(x);
   }

  // Add non-singular part of the solution
  u+=u_non_singular(x);
  
  return u;
 }


 /// Non-singular part of the solution for testing
 void u_and_gradient_non_singular(const Vector<double>& x,
                                  double& u,
                                  Vector<double>& dudx)
 {
  u=0.1*(sin(MathematicalConstants::Pi * x[0]/2.0)*
         sinh(MathematicalConstants::Pi * x[1]/2.0));
  dudx[0]=MathematicalConstants::Pi/2.0*
   0.1*(cos(MathematicalConstants::Pi * x[0]/2.0)*
        sinh(MathematicalConstants::Pi * x[1]/2.0));
  dudx[1]=MathematicalConstants::Pi/2.0*
   0.1*(sin(MathematicalConstants::Pi * x[0]/2.0)*
        cosh(MathematicalConstants::Pi * x[1]/2.0));
 }


 /// Non-singular part of the solution for testing
 double u_non_singular(const Vector<double>& x)
 { 
  double u=0.0;
  Vector<double> dudx(2);
  u_and_gradient_non_singular(x,u,dudx);
  return u;
 }


 /// Exact solution as vector
 void u_exact_as_vector(const Vector<double>& x, Vector<double>& u)
 {
  u[0]=u_exact(x);
 }




 /// Flux
 void prescribed_flux(const Vector<double>& x, 
                      double& flux)
 {

  // ici fix this!

  double tol=1.0e-6;
  Vector<double> normal(2);

  // rightmost boundary
  if(fabs(x[0]-2.0-2.0*Scaling_factor_for_domain)<
     Scaling_factor_for_domain*tol)
  {
    normal[0]=1.0;
    normal[1]=0.0;
  }

  // Vertical boundary below the step
  if(fabs(x[0]-2.0)<Scaling_factor_for_domain*tol)
  {
    normal[0]=-1.0;
    normal[1]=0.0;
  }


  // 
  if(fabs(x[0]-2.0*(1.0-Scaling_factor_for_domain))<Scaling_factor_for_domain*tol)
  {
    normal[0]=-1.0;
    normal[1]=0.0;
  }

  if(fabs(x[1]-1.0*Scaling_factor_for_domain)<Scaling_factor_for_domain*tol)
  {
    normal[0]=0.0;
    normal[1]=1.0;
  }

  if(fabs(x[1])<Scaling_factor_for_domain*tol)
  {
    normal[0]=0.0;
    normal[1]=-1.0;
  }

  if(fabs(x[1]+1.0*Scaling_factor_for_domain)<Scaling_factor_for_domain*tol)
  {
    normal[0]=0.0;
    normal[1]=-1.0;
  }


  Vector<double> gradient_non_sing_part(2);

  double dummy_u;

  u_and_gradient_non_singular(x,dummy_u,gradient_non_sing_part);

  Vector<double> gradient_sing_part(2);

  singular_fct_and_gradient(x,dummy_u,gradient_sing_part);

  flux=normal[0]*(gradient_sing_part[0] + gradient_non_sing_part[0]);
  flux+=normal[1]*(gradient_sing_part[1] + gradient_non_sing_part[1]);
 }



 /// Function to specify boundary conditions; currently
 /// the exact solution but can be changed
 double u_BC(const Vector<double>& x)
 {

  // if(fabs(x[0]-2.0)<0.001 && x[1]<=0.0)
  // {
  //   // if(x[1]>-0.5)
  //   // {
  //   //   return 0.0;
  //   // }
  //   // else
  //   // {
  //   //   double X = -(2.0)*x[1] - 1.0;
  //   //   return 1.0-(2.0*X*X*X - 3.0*X*X + 1.0);
  //   // }
  //   return 0.0;
  // }

  // if(fabs(x[1])<0.001&& x[0] <= 2.0)
  // {
  //   // if(x[0]>1.5)
  //   // {
  //   //   return 0.0;
  //   // }
  //   // else
  //   // {
  //   //   double X = 1.0-(2.0/3.0)*x[0];
  //   //   return 1.0-(2.0*X*X*X - 3.0*X*X + 1.0);      
  //   // }
  //   return 0.0;
  // }

  // if(fabs(x[0]-4.0)<0.001)
  // {
  //   return 1.0;
  // }

  // if(fabs(x[1]-1.0)<0.001)
  // {
  //   return 1.0;
  // }

  // if(fabs(x[0])<0.001)
  // {
  //   return -2.0*x[1]*x[1]*x[1]+3.0*x[1]*x[1];
  // }

  // if(fabs(x[1]+1.0)<0.001)
  // {
  //   return -(x[0]-2.0)*(x[0]-2.0)*(x[0]-2.0)/4.0 + 3.0*(x[0]-2.0)*(x[0]-2.0)/4.0;
  // }
  // else
  //   {return 0.0;}

  // return (x[0]-2.0)*x[1]*sin(0.2*(x[0]+x[1]));

  return u_exact(x);

 }

} // end_of_namespace



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////



//==start_of_problem_class============================================
/// Poisson in backward-facing step. Dirichlet boundary conditions along
/// all boundaries.
//====================================================================
template<class ELEMENT>
class StepProblem : public Problem
{

public:


  /// Constructor
  StepProblem();

  /// Destructor 
 ~StepProblem()
  {
   // hierher: at some point delete things properly and do memory leak
   // check
   delete Bulk_mesh_pt->spatial_error_estimator_pt();
   delete Bulk_mesh_pt;
  }
 
 /// Update the after solve (empty)
 void actions_after_newton_solve(){}
 
 /// \short Update the problem specs before solve (empty)
 void actions_before_newton_solve() {}
 
 // Perform actions after mesh adaptation
 void actions_after_adapt()
  {
   // Recreate face elements
   create_face_elements();
   
   // Complete problem setup
   complete_problem_setup();

   // Rebuild global mesh
   rebuild_global_mesh();
  }
 
 /// Perform actions after mesh adaptation (empty)
 void actions_before_adapt()
  {
   // Kill face elements
   delete_face_elements();

   // Rebuild global mesh
   rebuild_global_mesh();
  }
 
 /// Access function for the specific mesh
 RefineableTriangleMesh<ELEMENT>* mesh_pt()
  {
   return dynamic_cast<RefineableTriangleMesh<ELEMENT>*>(Problem::mesh_pt());
  }
 
 /// Doc the solution
 void doc_solution(DocInfo& doc_info);

 /// Validation
 void impose_amplitude_runs(DocInfo& doc_info);

 /// hierher more validation
 void check_residual_for_exact_non_singular_fe_soln(DocInfo& doc_info);

 /// Check condition number
 void check_condition_number();

private:

 /// Do what it says
 void complete_problem_setup();
 
 /// Helper function to apply boundary conditions
 void apply_boundary_conditions();

 /// hierher Delete face elements and flush meshes
 void delete_face_elements()
  {

   // Loop over the flux elements
   unsigned n_element = Flux_boundary_condition_mesh_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill
     delete Flux_boundary_condition_mesh_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Flux_boundary_condition_mesh_pt->flush_element_and_node_storage();


   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
    {
     return;
    }


   // hierher
   // // Loop over the flux jump elements
   // unsigned n_element = Face_mesh_for_flux_jump_pt->nelement();
   // for(unsigned e=0;e<n_element;e++)
   //  {
   //   delete Face_mesh_for_flux_jump_pt->element_pt(e);
   //  }
   
   // // hierher: actually kill nodes too because they've been duplicated

   // // Wipe the mesh
   // Face_mesh_for_flux_jump_pt->flush_element_and_node_storage();



   // Loop over the bc elements
   n_element = Face_mesh_for_bc_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Kill
     delete Face_mesh_for_bc_pt->element_pt(e);
    }
   
   // Wipe the mesh
   Face_mesh_for_bc_pt->flush_element_and_node_storage();



   
   // Loop over the integral face elements
   n_element = Face_mesh_for_singularity_integral_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     delete Face_mesh_for_singularity_integral_pt->element_pt(e);
    }
   Face_mesh_for_singularity_integral_pt->flush_element_and_node_storage();

  }

 /// Create face elements
 void create_face_elements()
  {
   
   // Flux boundary conditions?
   if (Global_Physical_Variables::Do_flux_problem)
    {     
     // Flux boundariies
     for(unsigned i_bound = 0;i_bound<4;i_bound++)
      {
       unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
       for(unsigned e=0; e<n_element; e++)
        {
         //Create Pointer to bulk element adjacent to the boundary
         ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
          (Bulk_mesh_pt->boundary_element_pt(i_bound,e));
         
         //Get Face index of boundary in the bulk element
         int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound,e);
         
         //Create corresponding face element
         PoissonWithSingularityFluxElement<ELEMENT>* 
          flux_element_pt =
          new PoissonWithSingularityFluxElement<ELEMENT>(
           bulk_elem_pt,face_index);
         
         // Set the pointer to the prescribed flux function
         flux_element_pt->flux_fct_pt() = 
          &Global_Physical_Variables::prescribed_flux;
         
         if (!CommandLineArgs::command_line_flag_has_been_set
             ("--dont_use_singularity"))
          {
           // We pass the pointer of singular function element to the 
           // face element (Set function because it also declares 
           // the amplitude to be
           // external data for that element).
           flux_element_pt->set_poisson_sing_el_pt(
            dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
             Singular_fct_element_mesh_pt->element_pt(0)));
          }


         //Attach it to the mesh
         Flux_boundary_condition_mesh_pt->add_element_pt(flux_element_pt);
        }
      }
    }
   // Dirichlet boundary conditions
   else
    {
     if (CommandLineArgs::command_line_flag_has_been_set
         ("--enforce_dirichlet_bcs_by_lagrange_multipliers"))
      {
       // BC elements live on the whole outer boundary:
       {
        for(unsigned i_bound = 0;i_bound<7;i_bound++)
         {
          unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
          
          // Loop over the bulk elements adjacent to boundary b
          for(unsigned e=0;e<n_element;e++)
           {
            // Get pointer to the bulk element that is adjacent to boundary b
            ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
             Bulk_mesh_pt->boundary_element_pt(i_bound,e));
            //Find the index of the face of element e along boundary b 
            int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound,e);
            
            // Build the corresponding bc element
            PoissonWithSingularityBCFaceElement<ELEMENT>* bc_element_pt = new 
             PoissonWithSingularityBCFaceElement<ELEMENT>
             (bulk_elem_pt,face_index,BC_el_id);
            
            // Tell the element about the singular fct
            if (!CommandLineArgs::command_line_flag_has_been_set
                ("--dont_use_singularity"))
             {
              bc_element_pt->set_poisson_sing_el_pt(
               dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
                Singular_fct_element_mesh_pt->element_pt(0)));
             }
            
            //Add the bc element to the surface mesh
            Face_mesh_for_bc_pt->add_element_pt(bc_element_pt);
           }
         }
       }
      }
    }


   if (CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
    {
     return;
    }


   // hierher renable


   // // Map keeps a running count of duplicate nodes already created;
   // // existing_duplicate_node_pt[orig_node_pt]=new_node_pt.
   // std::map<Node*,Node*> existing_duplicate_node_pt;


   // Flux jump elements on boundary 7 of region 1:
   //----------------------------------------------
   // NOTE: Since these duplicate nodes, these elements must be
   //----------------------------------------------------------
   //       constructed first!
   //       ------------------
   // {
   //  // Where are we? Inside region 1 on boundary 7
   //  unsigned b=7;
   //  unsigned region_id=1;
   //  unsigned nel=Bulk_mesh_pt->nboundary_element_in_region(b,region_id);
   //  for (unsigned e=0;e<nel;e++)
   //   {
   //    FiniteElement* el_pt=
   //     Bulk_mesh_pt->boundary_element_in_region_pt(b,region_id,e);
      
   //    // What is the index of the face of the bulk element at the boundary
   //    int face_index = Bulk_mesh_pt->
   //     face_index_at_boundary_in_region(b,region_id,e);
      
   //    // Build the corresponding flux jump element
   //    PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* flux_jump_element_pt 
   //     = new PoissonWithSingularityFluxJumpFaceElement<ELEMENT>
   //     (el_pt,face_index,existing_duplicate_node_pt,Flux_jump_el_id);

   //    //Add the flux jump element to the mesh
   //    Face_mesh_for_flux_jump_pt->add_element_pt(flux_jump_element_pt);
   //   }
   // }
   
   // Now add all new (duplicated) nodes to mesh
   // for (std::map<Node*,Node*>::iterator it=
   //       existing_duplicate_node_pt.begin();
   //      it!=existing_duplicate_node_pt.end();it++)
   //  {
   //   Face_mesh_for_flux_jump_pt->add_node_pt((*it).second);
   //  }






   
   // Create the face elements needed to compute the amplitude of
   // the singular function,
   
   // Only outer boundaries
   for(unsigned i_bound = 0;i_bound<7;i_bound++)
    {
      unsigned n_element = Bulk_mesh_pt->nboundary_element(i_bound);
      for(unsigned e=0; e<n_element; e++)
      {
       //Create Pointer to bulk element adjacent to the boundary
       ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>
         (Bulk_mesh_pt->boundary_element_pt(i_bound,e));
       
       //Get Face index of boundary in the bulk element
       int face_index = Bulk_mesh_pt->face_index_at_boundary(i_bound,e);
       
       //Create corresponding face element
       PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>* 
        flux_element_pt =
        new PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>(
         bulk_elem_pt,face_index);

       //We pass the pointer of singular function element to the face element
       flux_element_pt->poisson_sing_el_pt()=
        dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
         Singular_fct_element_mesh_pt->element_pt(0));
       
       //Attach it to the mesh
       Face_mesh_for_singularity_integral_pt->add_element_pt(flux_element_pt);
      }
    }
   
   // Update the pointer to the face elements (currently needed so
   // this GeneralisedElement can assemble the contributions to the
   // r_C residual from the face elements!
   dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
    Singular_fct_element_mesh_pt->element_pt(0))->
    set_mesh_of_face_elements(Face_mesh_for_singularity_integral_pt);



  

  
   // // Now loop over bulk elements 
   // //----------------------------
   // // and swap over any of their nodes have been replaced
   // //----------------------------------------------------
   // unsigned n_el=Bulk_mesh_pt->nelement();
   // for (unsigned e=0;e<n_el;e++)
   //  {
   //   ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(
   //    Bulk_mesh_pt->element_pt(e));
   
   //   // Loop over all nodes and check if they're amongst the replaced
   //   // ones
   //   unsigned nnod=bulk_el_pt->nnode();
   //   for (unsigned j=0;j<nnod;j++)
   //    {
   //     Node* nod_pt=bulk_el_pt->node_pt(j);
     
   //     // Find this original node in the map; if we find it
   //     // it's already been duplicated
   //     std::map<Node*,Node*>::iterator it=
   //      existing_duplicate_node_pt.find(nod_pt);
   //     if (it!=existing_duplicate_node_pt.end())
   //      {
   //       // Use the existing duplicate node
   //       bulk_el_pt->node_pt(j)=(*it).second;
   //      }
   //    }   
   //  }


  }
 
 /// Pointer to the bulk mesh
 RefineableTriangleMesh<ELEMENT> *Bulk_mesh_pt;
 
 /// Face element mesh for jump in interior of domain
 // hierher Mesh* Face_mesh_for_flux_jump_pt;

 /// Face element mesh for BC (Lagrange multiplier!) 
 Mesh* Face_mesh_for_bc_pt;

 /// \short Face elements used to compute the amplitude of the singular
 /// function
 Mesh* Face_mesh_for_singularity_integral_pt;
 
 /// Mesh for (single) element containing singular fct
 Mesh* Singular_fct_element_mesh_pt;
 
 /// Mesh of face elements applying flux bc 
 Mesh* Flux_boundary_condition_mesh_pt;

 /// \short Enumeration for IDs of FaceElements (used to figure out
 /// who's added what additional nodal data...)
 enum{Flux_jump_el_id,BC_el_id};
 
}; // end_of_problem_class


//==start_of_constructor==================================================
/// Constructor for StepProblem problem
//========================================================================
template<class ELEMENT>
StepProblem<ELEMENT>::StepProblem()
{

  // Build the mesh
  Bulk_mesh_pt=Global_Physical_Variables::build_the_mesh<ELEMENT>
   (Global_Physical_Variables::Uniform_element_area);


  // Set error estimator for bulk mesh
  Z2ErrorEstimator* error_estimator_pt=new Z2ErrorEstimator;
  Bulk_mesh_pt->spatial_error_estimator_pt()=error_estimator_pt;
  
  // Set element size limits
  Bulk_mesh_pt->max_element_size()=0.1;
  Bulk_mesh_pt->min_element_size()=1e-30;
  Bulk_mesh_pt->max_permitted_error()=0.005;
  Bulk_mesh_pt->min_permitted_error()=0.0;

  // Rescale to shrink domain
  if(Global_Physical_Variables::Scaling_factor_for_domain!=1.0)
  {
    unsigned nnode = Bulk_mesh_pt->nnode();
    for(unsigned i=0; i<nnode; i++)
    {
     Node* nod_pt = Bulk_mesh_pt->node_pt(i);
     nod_pt->x(0) = Global_Physical_Variables::Scaling_factor_for_domain*
      (nod_pt->x(0)-2.0) +2.0;
     nod_pt->x(1) = Global_Physical_Variables::Scaling_factor_for_domain*
      (nod_pt->x(1));
    }
  }
    
  // Let's have a look at the boundary enumeration
  Bulk_mesh_pt->output_boundaries("boundaries.dat");

  // Add sub-mesh
  add_sub_mesh(Bulk_mesh_pt);

  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    
    // Create element that stores the singular fct and its amplitude
    //---------------------------------------------------------------
    ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
     new ScalableSingularityForPoissonElement<ELEMENT>;
    
    // Pass fct pointers:
    el_pt->unscaled_singular_fct_pt()
     =&Global_Physical_Variables::singular_fct;
    el_pt->gradient_of_unscaled_singular_fct_pt()=
     &Global_Physical_Variables::gradient_of_singular_fct;
    
    // Add to mesh
    Singular_fct_element_mesh_pt=new Mesh;
    Singular_fct_element_mesh_pt->add_element_pt(el_pt);
    add_sub_mesh(Singular_fct_element_mesh_pt);


    // Create face elements that compute contribution
    //-----------------------------------------------
    // to r_c
    //-------
    Face_mesh_for_singularity_integral_pt=new Mesh; 

    
    // Create face elements for flux jump
    //-----------------------------------
    // hierher Face_mesh_for_flux_jump_pt=new Mesh;
    

   }
  
  // Create face elements for imposition of BC
  Face_mesh_for_bc_pt=new Mesh;
  
  // Flux boundary condition 
  Flux_boundary_condition_mesh_pt=new Mesh;

  // Build the face elements
  create_face_elements();
  
  // Add to mesh
  add_sub_mesh(Face_mesh_for_bc_pt);
  add_sub_mesh(Flux_boundary_condition_mesh_pt);


  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    // hierher add_sub_mesh(Face_mesh_for_flux_jump_pt);

    
    // hierher currently not needed because the contributiions from these elements are accumulated in
    // Singular_fct_element_mesh_pt but will need to add this back in 
    // when we optimise the code
    // add_sub_mesh(Face_mesh_for_singularity_integral_pt); 
   }

  // Build global mesh
  build_global_mesh();
  
  // Complete problem setup
  complete_problem_setup();
    


  // // hierher kill
  // {
  //  oomph_info << "imposing amplitude; removve thiis!\n";
   
  //  ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
  //   dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>
  //   (Singular_fct_element_mesh_pt->element_pt(0));
   
   
  //  // Change r_C so that C is assigned directly
  //  double imposed_amplitude=1.0;
  //  el_pt->impose_singular_fct_amplitude(imposed_amplitude);
  // }


  // Setup equation numbering scheme
  oomph_info <<"Number of equations: " 
             << this->assign_eqn_numbers() 
             << std::endl;
  
} // end_of_constructor


//==start_of_complete======================================================
/// Set boundary condition, and complete the build of
/// all elements
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::complete_problem_setup()
{
 
 if (!CommandLineArgs::command_line_flag_has_been_set
       ("--dont_use_singularity"))
  {
   
   // Loop over the elements to set up element-specific
   // things that cannot be handled by constructor
   unsigned n_el=Bulk_mesh_pt->nelement();
   for (unsigned e=0;e<n_el;e++)
    {
     ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(
      Bulk_mesh_pt->element_pt(e));
     
     // Tell the bulk element about the singular fct
     bulk_el_pt->poisson_sing_el_pt()=
      dynamic_cast<TemplateFreeScalableSingularityForPoissonElement*>(
       Singular_fct_element_mesh_pt->element_pt(0));
     bulk_el_pt->exact_non_singular_fct_pt()=
      &Global_Physical_Variables::u_and_gradient_non_singular;
    }
   

   // hierher move this to create_face_elements as for the others
   // Flux jump elements
   // unsigned n_element = Face_mesh_for_flux_jump_pt->nelement();
   // for(unsigned e=0;e<n_element;e++)
   //  {
   //   // Upcast from GeneralisedElement to the present element
   //   PoissonWithSingularityFluxJumpFaceElement<ELEMENT>* el_pt = dynamic_cast<
   //    PoissonWithSingularityFluxJumpFaceElement<ELEMENT>*>(
   //     Face_mesh_for_flux_jump_pt->element_pt(e));
     
   //   // Tell the element about the singular fct
   //   el_pt->set_poisson_sing_el_pt(
   //    dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
   //     Singular_fct_element_mesh_pt->element_pt(0)));
   //  }
  


   //BC elements
   unsigned n_element =  Face_mesh_for_bc_pt->nelement();
   for(unsigned e=0;e<n_element;e++)
    {
     // Upcast from GeneralisedElement to the present element
     PoissonWithSingularityBCFaceElement<ELEMENT>* el_pt = dynamic_cast<
      PoissonWithSingularityBCFaceElement<ELEMENT>*>(
       Face_mesh_for_bc_pt->element_pt(e));
     
     // Tell the element about the singular fct
     el_pt->set_poisson_sing_el_pt(
      dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>(
       Singular_fct_element_mesh_pt->element_pt(0)));
    }


  }
 
 // Apply bcs
 apply_boundary_conditions();
 
}

//==start_of_apply_bc=====================================================
/// Helper function to apply boundary conditions
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::apply_boundary_conditions()
{
 

 //---------------------------------------------------------------
 // DEFAULT SETUP FOR STANDARD POISSON PROBLEM WITHOUT LAGRANGE
 // MULTIPLIERS (CHANGED BELOW IF THE SINGULARITY IS INCLUDED)
 //---------------------------------------------------------------


 // Set the boundary conditions for this problem: All nodes are
 // free by default -- just pin the ones that have Dirichlet conditions
 // here.
 unsigned num_bound = Bulk_mesh_pt->nboundary();
 for(unsigned ibound=0;ibound<num_bound;ibound++)
  {   
   bool pin_it=false;
   
   // Flux problem
   if (Global_Physical_Variables::Do_flux_problem)
    {
     if ((ibound==4)||(ibound==5)||(ibound==6))
      {
       pin_it=true;
      }
    }
   // Otherwise pin everywhere but leave internal boundary alone
   else
    {
     if (ibound!=7) 
      {
       pin_it=true;
      }
    }
   // Pin it?
   if (pin_it)
    { 
     unsigned num_nod=Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Bulk_mesh_pt->boundary_node_pt(ibound,inod)->pin(0);
      }
    }
  } // end loop over boundaries
 
 
 // Now set boundary values
 for (unsigned ibound=0;ibound<num_bound;ibound++)
  {
   // Leave internal boundary alone
   if (ibound!=7)
    {
     unsigned num_nod= Bulk_mesh_pt->nboundary_node(ibound);
     for (unsigned inod=0;inod<num_nod;inod++)
      {
       Node* nod_pt=Bulk_mesh_pt->boundary_node_pt(ibound,inod);
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       double u=Global_Physical_Variables::u_BC(x);
       nod_pt->set_value(0,u);
      }
    }
  }



 // Free vaues on boundary if Lagrange multiplier is used
 // to enforce bc for sum of FE soln and singuar fct
 // hierher
 // if (!CommandLineArgs::command_line_flag_has_been_set
 //       ("--dont_use_singularity"))
  {

   // Now unpin nodal values where the bc conditions are enforceed
   // by Lagrange multiplier to ensure that the sum of fe and singular
   // solution is correct
   unsigned nel=Face_mesh_for_bc_pt->nelement();
   for (unsigned e=0;e<nel;e++)
    {
     // Get element
     PoissonWithSingularityBCFaceElement<ELEMENT>* el_pt=
      dynamic_cast<PoissonWithSingularityBCFaceElement<ELEMENT>*>(
       Face_mesh_for_bc_pt->element_pt(e));
     
     // Specify desired nodal values for compound solution
     unsigned nnod=el_pt->nnode();
     Vector<double> nodal_boundary_value(nnod,0.0);

     // Unpin the fe part of the solution
     for (unsigned j=0;j<nnod;j++)
      {
       el_pt->unpin_u_fe_at_specified_local_node(j);
       
       // Now here's another subtle one! If this node stores three values
       // they come from: original Poisson dof; one Lagrange multiplier
       // to ensure continuity of solution across "jump"; one Lagrange
       // multiplier from enforcing Dirichlet boundary conditions. 
       // The latter two enforce the same thing, so we can only apply
       // the constraint once --> Pin the Lagrange multiplier that tries to 
       // enforce the Dirichlet BC
       unsigned nval=el_pt->node_pt(j)->nvalue();
       if (nval==3)
        {
         oomph_info << "hierher this shouldn't currently be used \n"
                    << " and ought to be double checked when \n"
                    << "flux jump comes back in...\n";
         abort();
         el_pt->pin_lagrange_multiplier_at_specified_local_node(j);
        }
       
       Node* nod_pt=el_pt->node_pt(j);
       Vector<double> x(2);
       x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
       double u=Global_Physical_Variables::u_BC(x);
       nodal_boundary_value[j]=u;
      }
     
     // Tell the element about the desired nodal boundary values
     el_pt->set_nodal_boundary_values(nodal_boundary_value);
    }

  }



} // end set bc


//==start_of_doc_solution=================================================
/// Doc the solution
//========================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::doc_solution(DocInfo& doc_info)
{

  ofstream some_file;
  char filename[100];

  // Number of plot points
  unsigned npts=10;

  // Output solution
  sprintf(filename,"%s/soln%i.dat",doc_info.directory().c_str(),
  doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();


  // Output solution just using vertices so we can see the mesh
  sprintf(filename,"%s/coarse_soln%i.dat",doc_info.directory().c_str(),
  doc_info.number());
  some_file.open(filename);
  npts=2;
  Bulk_mesh_pt->output(some_file,npts);
  some_file.close();
  

  // Plot "extended solution" showing contributions; also work out
  // average element size
  double av_el_size=0.0;
  sprintf(filename,"%s/extended_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned nel=Bulk_mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    npts=10;
    ELEMENT* el_pt=
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(e));
    el_pt->output_with_various_contributions(some_file,npts);
    av_el_size+=el_pt->size();
  }
  some_file.close();
  av_el_size/=double(nel);

  // Get error
  double error,norm; 
  sprintf(filename,"%s/error%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  Bulk_mesh_pt->compute_error(some_file,
                              Global_Physical_Variables::u_exact_as_vector,
                              error,norm);
  some_file.close();

  // Doc error norm:
  oomph_info << "\n av el size, av h, Ndof, # bulk els, Norm of error    : "   
             << av_el_size << " " 
             << sqrt(av_el_size) << " " 
             << ndof() << " " 
             << Bulk_mesh_pt->nelement() << " " 
             << sqrt(error) << std::endl;
  oomph_info << "Norm of solution : " << sqrt(norm) << std::endl << std::endl;
  oomph_info << std::endl;
  
  // Exact solution
  sprintf(filename,"%s/exact_soln%i.dat",doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  unsigned nplot=5;
  Bulk_mesh_pt->output_fct(some_file,nplot,
                           Global_Physical_Variables::u_exact_as_vector);
  some_file.close();


  // Non-singular part of the solution
  sprintf(filename,"%s/non_singular_part_of_soln%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  nplot=5;
  Bulk_mesh_pt->output_fct(some_file,nplot,
                           Global_Physical_Variables::u_non_singular_as_vector);
  some_file.close();



  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    // hierher put into one file and flush 
    // trace file for va lue of C
    sprintf(filename,"%s/suivi_C.dat",doc_info.directory().c_str());
    some_file.open(filename,std::ios::out | std::ios::app);
    some_file << sqrt(av_el_size) <<" " 
              << dynamic_cast<TemplateFreeScalableSingularityForPoissonElement*>(
               Singular_fct_element_mesh_pt->element_pt(0)
               )->amplitude_of_singular_fct()<<endl;
    some_file.close();
   }


  // DoubleVector r;
  // CRDoubleMatrix jac;

  // sprintf(filename,"%s/most_recent_jacobian%i.dat",doc_info.directory().c_str(),doc_info.number());
  // some_file.open(filename);
  // oomph_info << "SETTING UP JAC FOR OUTPUT OF MOST RECENT JAC\n";
  // get_jacobian(r,jac);
  // oomph_info << "DONE SETTING UP JAC FOR OUTPUT OF MOST RECENT JAC\n";
  // jac.sparse_indexed_output(some_file);



  // sprintf(filename,"%s/suivi_derivee.dat",doc_info.directory().c_str());
  // some_file.open(filename,std::ios::out | std::ios::app);
  // // We're looking for a bulk element that has a node on the
  // // the corner of the backward step; all candidate elements live on
  // // boundary 4
  // unsigned b = 4;
  // unsigned n_element = Bulk_mesh_pt->nboundary_element(b);
   
  // // Loop over the bulk elements adjacent to boundary b
  // for(unsigned e=0;e<n_element;e++)
  // {
  //  // Get pointer to the bulk element that is adjacent to boundary b
  //  ELEMENT* bulk_elem_pt = dynamic_cast<ELEMENT*>(
  //   Bulk_mesh_pt->boundary_element_pt(b,e));
  //  unsigned sortie = 0;
     
  //  unsigned nnod=bulk_elem_pt->nnode();
  //  for (unsigned j=0;j<nnod;j++)
  //  {
  //   Node* nod_pt=bulk_elem_pt->node_pt(j);
      
  //   // Are we at the corner?
  //   double distance=sqrt((nod_pt->x(0)-Global_Physical_Variables::L_up)*
  //                         (nod_pt->x(0)-Global_Physical_Variables::L_up)+
  //                         (nod_pt->x(1)-0.0)*(nod_pt->x(1)-0.0));
  //   double tol=1.0e-10;
  //   if (distance<tol)
  //   {
  //    Vector<double> outward_vector(2);
  //    outward_vector[0] = 1.0/sqrt(2.0);
  //    outward_vector[1] = 1.0/sqrt(2.0);
  //    Vector<double> Flux(1);
  //    bulk_elem_pt->get_flux(outward_vector,Flux);
  //    some_file << sqrt(av_el_size) << " " << Flux[0] << endl;
  //    sortie = 1;
  //    break;
  //   }
  //  }
  //  if (sortie==1) break;
  // }
  // some_file.close();


  if (!CommandLineArgs::command_line_flag_has_been_set
      ("--dont_use_singularity"))
   {
    oomph_info 
     << "Amplitude of singular function: "
     << dynamic_cast<TemplateFreeScalableSingularityForPoissonElement*>(
      Singular_fct_element_mesh_pt->element_pt(0))->
     amplitude_of_singular_fct() << std::endl;


    // Output face elements used to compute amplitude of singularity
    sprintf(filename,"%s/integrand_elements%i.dat",
            doc_info.directory().c_str(),
            doc_info.number());
    some_file.open(filename);
    nel=Face_mesh_for_singularity_integral_pt->nelement();
    for (unsigned e=0;e<nel;e++)
     {
      PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>* el_pt=
       dynamic_cast<PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>*>
       (Face_mesh_for_singularity_integral_pt->element_pt(e));
      el_pt->get_contribution_integral(some_file); //output(some_file,npts);
     }
    some_file.close();
   }
  
  
  // Output flux elements
  sprintf(filename,"%s/flux_elements%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file.open(filename);
  nel=Flux_boundary_condition_mesh_pt->nelement();
  for (unsigned e=0;e<nel;e++)
   {
    npts=5;
    Flux_boundary_condition_mesh_pt->finite_element_pt(e)->output(some_file,npts);
   }
  some_file.close();


} // end_of_doc_solution




//=======================================================================
// hierher
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::check_condition_number()
{
 
 // Fake r_c equation?
 if (Global_Physical_Variables::Problem_type_for_check_condition_number==2)
  {
   ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
    dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>
    (Singular_fct_element_mesh_pt->element_pt(0));
   
   // Change r_C so that C is assigned directly
   double imposed_amplitude=1.0;
   el_pt->impose_singular_fct_amplitude(imposed_amplitude);
  }


 // Storage for the Jacobian
 CRDoubleMatrix matrix;
 
 // Storage for the (to be unused) residual vector
 DoubleVector residual;
 
 // Storage for the eigenvector associated with the largest eigenvalue
 DoubleVector eigenvector_max;
 
 // Storage for the eigenvector associated with the smallest eigenvalue
 DoubleVector eigenvector_min;
 
 // Get the Jacobian matrix and residual vector
 get_jacobian(residual,matrix);
 
 // Tell the user
 oomph_info << "\n=================================================="
	    << "\nBeginning eigenvalue/singular value calculation..."
	    << "\n=================================================="
	    << std::endl;
 
 // Storage for the largest eigenvalue/singular value
 double lambda_max=power_method(&matrix,eigenvector_max);
 
 // Storage for the smallest eigenvalue/singular value
 double lambda_min=inverse_power_method(&matrix,eigenvector_min);
 
 // Set the output precision
 oomph_info.stream_pt()->precision(10);
 
 // If we're just computing eigenvalues
 if (!PowerIterationHelperNamespace::Compute_singular_values)
  {
   // Output the eigenvalues
   oomph_info << "\nLargest eigenvalue: " << lambda_max
              << "\nSmallest eigenvalue: " << lambda_min
              << std::endl;
  }
 // If we're computing singular values instead
 else
  {
   // Output the singular values
   oomph_info << "\nLargest singular value: " << lambda_max
              << "\nSmallest singular value: " << lambda_min
              << "\nMatrix size and condition number: " 
              << residual.nrow() << " " << lambda_max/lambda_min
              << std::endl;
  }
 
 // Create a method of outputting the Jacobian
 // std::ofstream matrix_file("jacobian.dat");
 
 // Dump the matrix so we can check that we're doing the right thing...
 //matrix.sparse_indexed_output(matrix_file);
 
 // Close the file
 // matrix_file.close();
 
 // Output the eigenvector associated with the largest eigenvalue
 eigenvector_max.output("eigenvector_max.dat");
 
 // Also output the eigenvector associated with the smallest eigenvalue
 eigenvector_min.output("eigenvector_min.dat"); 
 

}


//=======================================================================
// hierher
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::check_residual_for_exact_non_singular_fe_soln
(DocInfo& doc_info)
{
 
 oomph_info << std::endl << std::endl;

 // Assign exact non-singular solution to nodal values
 unsigned nnod=Bulk_mesh_pt->nnode();
 for (unsigned j=0;j<nnod;j++)
  {
   Node* nod_pt=Bulk_mesh_pt->node_pt(j);
   Vector<double> x(2);
   x[0]=nod_pt->x(0);
   x[1]=nod_pt->x(1);
   double u=Global_Physical_Variables::u_non_singular(x); // hierher break deliberately...
   nod_pt->set_value(0,u);
  }

 // STAGE 1: Evaluate the integral with exact non-singular part
 //------------------------------------------------------------

 // Check the computation of the integral that determines C.
 // If we set u_fe and its gradient to the non-singular part of
 // the solution this should be close to zero!
 unsigned nel=Face_mesh_for_singularity_integral_pt->nelement();
 for (unsigned e=0;e<nel;e++)
  {
   PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>* el_pt=
    dynamic_cast<PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>*>
    (Face_mesh_for_singularity_integral_pt->element_pt(e));
   el_pt->exact_non_singular_fct_pt()=
    &Global_Physical_Variables::u_and_gradient_non_singular;
  }
 
 
 ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
  dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>
  (Singular_fct_element_mesh_pt->element_pt(0));
 unsigned n_dof=el_pt->ndof();
 Vector<double> residuals(n_dof,0.0);  
 el_pt->fill_in_contribution_to_residuals(residuals);
 
 oomph_info << "Residual if we assign the exact non-singular part of the "
            << "solution for area "
            << Global_Physical_Variables::Uniform_element_area 
            << " and number of elements: " 
            << Bulk_mesh_pt->nelement() << " : " 
            << fabs(residuals[0]) << endl;
 
 doc_solution(doc_info);
 doc_info.number()++;

 // Reset
 for (unsigned e=0;e<nel;e++)
  {
   PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>* el_pt=
    dynamic_cast<PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>*>
    (Face_mesh_for_singularity_integral_pt->element_pt(e));
   el_pt->exact_non_singular_fct_pt()=0;
  }


 // STAGE 2: Compute integral from (previously assigned) exact nodal values
 //------------------------------------------------------------------------
 residuals.clear();
 residuals.resize(n_dof);
 el_pt->fill_in_contribution_to_residuals(residuals);
 
 oomph_info << "Residual if we assign the exact non-singular part of the "
            << "solution to nodal values for area "
            << Global_Physical_Variables::Uniform_element_area 
            << " and number of elements: " 
            << Bulk_mesh_pt->nelement() << " : " 
            << fabs(residuals[0]) << endl;
 
 doc_solution(doc_info);
 doc_info.number()++;
 
}



//=======================================================================
// Impose amplitude of singular fct via fake eqn
//=======================================================================
template<class ELEMENT>
void StepProblem<ELEMENT>::impose_amplitude_runs(DocInfo& doc_info)
{

  ofstream some_file;
  ofstream some_file_error;
  ofstream some_file_integral_of_error;
  char filename[100];

  sprintf(filename,"%s/rc_vs_amplitude.dat",
          doc_info.directory().c_str());
  some_file.open(filename);

  sprintf(filename,"%s/integral_of_Z2_error.dat",
          doc_info.directory().c_str());
  some_file_integral_of_error.open(filename);


  // Loop over imposed amplitudes
  double a_min=-1.0;
  unsigned nstep_aux=4;
  double d_ampl=(1.0-a_min)/double(nstep_aux);
  unsigned nstep=2*nstep_aux;
  double imposed_amplitude=a_min;
  for(unsigned incr = 0; incr < nstep; incr++)
  {
   double integral_of_Z2error = 0.0;
    
   sprintf(filename,"%s/elementwise_Z2error%i.dat",
          doc_info.directory().c_str(),
           doc_info.number());
   some_file_error.open(filename);
   

   ScalableSingularityForPoissonElement<ELEMENT>* el_pt=
    dynamic_cast<ScalableSingularityForPoissonElement<ELEMENT>*>
    (Singular_fct_element_mesh_pt->element_pt(0));

   // Change r_C so that C is assigned directly
   el_pt->impose_singular_fct_amplitude(imposed_amplitude);
   some_file << imposed_amplitude << " ";
   some_file_integral_of_error << imposed_amplitude << " ";
   
   // Solve the eqns
   this->newton_solve();

   oomph_info << "have solved for imposed_amplitude = "
              << imposed_amplitude << " output soln in "
              << doc_info.number() << std::endl;

   doc_solution(doc_info);
   doc_info.number()++;

   // Check residual (computed by actual integral!) in singular element  
   el_pt->dont_impose_singular_fct_amplitude();
  
  unsigned n_dof=el_pt->ndof();
  Vector<double> residuals(n_dof,0.0);  
  el_pt->fill_in_contribution_to_residuals(residuals);

  some_file <<  residuals[0] << endl;

  unsigned n = Bulk_mesh_pt->nelement();
  Vector<double> elementwise_Z2error(n);
  Mesh* mesh_copy_pt = Bulk_mesh_pt;
  
  // Use actual value without normalisation!
  Z2ErrorEstimator* z2_pt=dynamic_cast<Z2ErrorEstimator*>(
   Bulk_mesh_pt->spatial_error_estimator_pt());
  double backup=z2_pt->reference_flux_norm();
  z2_pt->reference_flux_norm()=1.0;
  
  //Bulk_mesh_pt->spatial_error_estimator_pt()->
   z2_pt->get_element_errors(mesh_copy_pt,elementwise_Z2error);

  // Reset
   z2_pt->reference_flux_norm()=backup;
    
  // Plot constant pressure th'out element
  unsigned nplot = 3;
  for(unsigned i=0; i<n;i++)
  {
   FiniteElement* el_pt=
     dynamic_cast<ELEMENT*>(Bulk_mesh_pt->element_pt(i));
   unsigned num_plot_points=el_pt->nplot_points(nplot);
   Vector<double> s(2);
   some_file_error << el_pt->tecplot_zone_string(nplot);
   for(unsigned j=0; j<num_plot_points;j++)
   {
    el_pt->get_s_plot(j,nplot,s);
    Vector<double> x(2);
    x[0]=el_pt->interpolated_x(s,0);
    x[1]=el_pt->interpolated_x(s,1);
    some_file_error << x[0] << " " << x[1] << " " 
                    << elementwise_Z2error[i] << endl;
   }
   el_pt->write_tecplot_zone_footer(some_file_error,nplot);
   integral_of_Z2error += elementwise_Z2error[i]*el_pt->size();
  }

  some_file_integral_of_error << integral_of_Z2error << endl;
  some_file_error.close();
  imposed_amplitude += d_ampl;
  }
  
  some_file_integral_of_error.close();
  some_file.close();
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//==start_of_main======================================================
/// Driver for backward step with impedance outflow bc
//=====================================================================
int main(int argc, char **argv)
{

 // Store command line arguments
 CommandLineArgs::setup(argc,argv);
  
 // Don't subtract off singularity
 CommandLineArgs::specify_command_line_flag(
  "--dont_use_singularity");

 // Suppress singular term in exact solution
 CommandLineArgs::specify_command_line_flag(
  "--suppress_sing_in_exact_soln");

 // Do runs where we impose the amplitude via the "fake" r_c equation
 CommandLineArgs::specify_command_line_flag(
  "--check_overall_behaviour_when_setting_amplitude_via_fake_rc_equation");

 // Check condition number and specify number that identifies which
 // case we're doing: 0: real problem; 1: FE-only; 2: fake r_c equation
 CommandLineArgs::specify_command_line_flag(
  "--check_condition_number",
  &Global_Physical_Variables::Problem_type_for_check_condition_number);

 // Check the r_c equation
 CommandLineArgs::specify_command_line_flag
  ("--check_rc_equation");

 // Scaling factor for domain (defaults to 1.0)
 CommandLineArgs::specify_command_line_flag(
  "--scaling_factor_for_domain",
  &Global_Physical_Variables::Scaling_factor_for_domain);
 
 // Use Dirichlet boundary conditions (allegedly makes condition
 // number (even) worse)
 CommandLineArgs::specify_command_line_flag(
  "--enforce_dirichlet_bcs_by_lagrange_multipliers");

 // Use Dirichlet boundary conditions (allegedly makes condition
 // number (even) worse)
 CommandLineArgs::specify_command_line_flag(
  "--use_dirichlet_bcs");

 // Uniform element size in mesh
 CommandLineArgs::specify_command_line_flag(
  "--uniform_element_area",
  &Global_Physical_Variables::Uniform_element_area);

 // // Max. number of adaptations in code
 // unsigned max_adapt = 0;
 // CommandLineArgs::specify_command_line_flag(
 //  "--max_adapt",&max_adapt);

 // Parse command line
 CommandLineArgs::parse_and_assign(); 
 
 // Doc what has actually been specified on the command line
 CommandLineArgs::doc_specified_flags();

 
 // hierher tidy
 if (CommandLineArgs::command_line_flag_has_been_set
     ("--use_dirichlet_bcs"))
  {
   Global_Physical_Variables::Do_flux_problem=false;
  }
 
 // Set up doc info
 // ---------------
 
 // Label for output
 DocInfo doc_info;
 
 // Set output directory
 doc_info.set_directory("RESLT");
 
 // Step number
 doc_info.number()=0;
 
  // Build the problem with 
 StepProblem<ProjectablePoissonElement<MyTPoissonElement<2,3> > > problem;




 // Check condition number
 //=======================
  if (CommandLineArgs::command_line_flag_has_been_set
      ("--check_condition_number"))
   {
    
    // Sorry this is a mess!
    if ((Global_Physical_Variables::Problem_type_for_check_condition_number==1)
     &&(!CommandLineArgs::command_line_flag_has_been_set
        ("--dont_use_singularity")))
      {
       oomph_info << "Please specify --dont_use_singularity\n";
       abort();
      }

    // Check condition number
    problem.check_condition_number();
    
    // Solve the bloody thing
    problem.newton_solve();
    problem.doc_solution(doc_info);
    exit(0);
   }



 // Impose amplitude via fake eqn?
 //================================
  if (CommandLineArgs::command_line_flag_has_been_set
       ("--check_overall_behaviour_when_setting_amplitude_via_fake_rc_equation"))
   {
    problem.impose_amplitude_runs(doc_info);
    exit(0);
   }
  

 // Check residual r_c
 //===================
  if (CommandLineArgs::command_line_flag_has_been_set
       ("--check_rc_equation"))
   {
    problem.check_residual_for_exact_non_singular_fe_soln(doc_info);
    exit(0);
   }
  

  cout << "not done yet\n";
  abort();

  

//##################



 // bool doc_resid=false;

 // Global_Physical_Variables::shrink_the_mesh=false;
 // Global_Physical_Variables::scale=0.5;

 // // Solve, refine uniformly and keep going
 //  for (unsigned i=0;i<max_adapt;i++)
 //   {

 //    if (doc_resid)
 //     {
 //      if(i==max_adapt-1)
 //      {
 //       problem.impose_amplitude_runs(doc_info);
 //       problem.check_residual_for_exact_non_singular_fe_soln();        
 //      }

 //     }
 //    else
 //     {
 //      // Solve the bloody thing
 //      problem.newton_solve();
 //      problem.doc_solution(doc_info);
 //      doc_info.number()++;
 //     }


 //    if (i!=(max_adapt-1)) problem.refine_uniformly();
 //   }


}
