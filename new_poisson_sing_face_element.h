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
// Header file for elements that are used to ...hierher
#ifndef OOMPH_POISSON_SING_FACE_ELEMENTS_HEADER
#define OOMPH_POISSON_SING_FACE_ELEMENTS_HEADER

namespace oomph
{

 //============================================================================
 // TemplateFreeScalableSingularityForPoissonElement defines the elements managing
 // the singular function : it is essentially a pointer to the singular function, 
 // its gradient and its amplitude
 //============================================================================
 class TemplateFreeScalableSingularityForPoissonElement :
  public virtual GeneralisedElement
  {
  public:

    typedef double(*UnscaledSingSolnFctPt) (const Vector<double>& x);


    typedef Vector<double>(*GradientOfUnscaledSingSolnFctPt) 
      (const Vector<double>& x);

    ///Constructor
    TemplateFreeScalableSingularityForPoissonElement()
    {
      //data to store amplitude
      add_internal_data(new Data(1));
    }

    ///Function to get pointer to unscaled version of singular function
    UnscaledSingSolnFctPt& unscaled_singular_fct_pt()
      {return Unscaled_singular_fct_pt;}

    ///Function to get pointer to unscaled version of gradient of singular function
    GradientOfUnscaledSingSolnFctPt& gradient_of_unscaled_singular_fct_pt() 
      {return Gradient_of_unscaled_singular_fct_pt;}

    ///Function to compute unscaled version of unscaled version
    double unscaled_singular_fct(const Vector<double>& x) const
    {
      if(Unscaled_singular_fct_pt==0){return 0.0;}
      return Unscaled_singular_fct_pt(x);
    }

    ///Compute unscaled version of gradient of singular function
    Vector<double> gradient_of_unscaled_singular_fct(const Vector<double>& x)
    const
    {
      Vector<double> grad;
      if(Gradient_of_unscaled_singular_fct_pt==0){return grad;}
      return Gradient_of_unscaled_singular_fct_pt(x);
    }

    ///Compute scaled version of singular function
    double singular_fct(const Vector<double>& x) const
    {
      return amplitude_of_singular_fct()*unscaled_singular_fct(x);
    }

    ///Compute scaled version of gradient of singular function
    Vector<double> gradient_of_singular_fct(const Vector<double>& x) const
    {
      Vector<double> grad(gradient_of_unscaled_singular_fct(x));
      unsigned n = grad.size();

      for(unsigned i=0;i<n;i++)
      {
        grad[i]*=amplitude_of_singular_fct();
      }
      return grad;
    }

    ///Access the amplitude of the singular function
    double amplitude_of_singular_fct() const
    {
      return data_that_stores_amplitude_of_singular_fct()
        ->value(index_of_value_that_stores_amplitude_of_singular_fct());
    }

    ///Set the amplitude of thz singular function
    void set_amplitude_of_singular_fct(const double& value)
    {
      data_that_stores_amplitude_of_singular_fct()
        ->set_value(index_of_value_that_stores_amplitude_of_singular_fct(),value);
    }

    ///pin amplitude of singular function
    void pin_amplitude_of_singular_fct()
    {
      data_that_stores_amplitude_of_singular_fct()
        ->pin(index_of_value_that_stores_amplitude_of_singular_fct());
    }

    ///Pointer to data that stores the amplitude of singular function
    Data* data_that_stores_amplitude_of_singular_fct() const
    {
      return internal_data_pt(0);
    }

    ///Gives the index of the amplitude value : default is 0
    unsigned index_of_value_that_stores_amplitude_of_singular_fct() const 
    {
      return 0;
    }

  private:

    ///Pointer to singular function
    UnscaledSingSolnFctPt Unscaled_singular_fct_pt;

    ///Pointer to gradient of singular funcion
    GradientOfUnscaledSingSolnFctPt Gradient_of_unscaled_singular_fct_pt;
  };



 
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//===========================================================================
/// PoissonWithSingularityBoundaryIntegralFaceElement is a class of face elements 
///used to compute the contribution to the residuals from the the singuar function
//===========================================================================
template <class ELEMENT>
class PoissonWithSingularityBoundaryIntegralFaceElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{
 
public:

 // hierher kill
 /* /// \short Function pointer to the prescribed-flux function fct(x,f(x)) --  */
 /* /// x is a Vector!  */
 /* typedef void (*PoissonPrescribedFluxFctPt) */
 /*  (const Vector<double>& x, double& flux); */

 /// \short Function pointer to the "exact" non-singular function fct(x,u,grad u)
 typedef void (*ExactNonSingularFctPt)
  (const Vector<double>& x, double& u, Vector<double>& grad_u);


 /// \short Constructor, takes the pointer to the "bulk" element and the 
 /// index of the face to which the element is attached.
 PoissonWithSingularityBoundaryIntegralFaceElement(FiniteElement* const &bulk_el_pt, 
                    const int& face_index);

 ///\short  Broken empty constructor
 PoissonWithSingularityBoundaryIntegralFaceElement()
  {
   throw OomphLibError(
    "Don't call empty constructor for PoissonWithSingularityBoundaryIntegralFaceElement",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Broken copy constructor
 PoissonWithSingularityBoundaryIntegralFaceElement(const PoissonWithSingularityBoundaryIntegralFaceElement& dummy) 
  { 
   BrokenCopy::broken_copy("PoissonWithSingularityBoundaryIntegralFaceElement");
  } 
 
 /// Broken assignment operator
 void operator=(const PoissonWithSingularityBoundaryIntegralFaceElement&) 
  {
   BrokenCopy::broken_assign("PoissonWithSingularityBoundaryIntegralFaceElement");
  }

 /// \short Specify the value of nodal zeta from the face geometry
 /// The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default (needed to break
 /// indeterminacy if bulk element is SolidElement)
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                          const unsigned &i) const 
  {return FaceElement::zeta_nodal(n,k,i);}

 /// Pointer to element that computes singular function related stuff
 TemplateFreeScalableSingularityForPoissonElement*& poisson_sing_el_pt() 
  { 
   return Poisson_sing_el_pt; 
  } 
 
 /// Add the element's contribution to its residual vector
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_poisson_flux(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 /// \short Add the element's contribution to its residual vector and its 
 /// Jacobian matrix
 inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                          DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_poisson_flux(residuals,jacobian,1);
  }

 /// Output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(std::ostream &outfile) {FiniteElement::output(outfile);}

 /// Output with various contributions
 void  output(std::ostream &outfile, 
              const unsigned &nplot)
  {
   //Vector of local coordinates
   Vector<double> s(Dim-1);
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);
     
     Vector<double> x(Dim);
     for(unsigned i=0;i<Dim;i++) 
      {
       x[i]=this->interpolated_x(s,i);
       outfile << x[i] << " ";
      }

     // Compute outer unit normal at the specified local coordinate
     Vector<double> unit_normal(Dim);
     outer_unit_normal(s,unit_normal);

     outfile << unit_normal[0] << " " << unit_normal[1] << " " << endl;

     // // Get gradient of FE solution from bulk element
     // ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());
     // Vector<double> flux(Dim);
     // Vector<double> s_bulk(Dim);
     // s_bulk=local_coordinate_in_bulk(s);
     // bulk_el_pt->get_flux(s_bulk,flux);
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
   
  }


 /// C-style output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(FILE* file_pt) {FiniteElement::output(file_pt);}

 /// \short C-style output function -- forward to broken version in 
 /// FiniteElement until somebody decides what exactly they want to plot 
 /// here...
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}



 /// \short Compute this element's contribution to the integral that determines C
 /// Also output into file
 double get_contribution_integral(std::ofstream &outfile);
 
 /// \short Compute this element's contribution to the integral that determines C
 double get_contribution_integral()
 {
  // Call with (non-open) dummy file
  ofstream some_file;
  return get_contribution_integral(some_file);
 }

 /// Pointer to exact non singular fct (and its gradient) only used
 /// to validate the computation of the integral. Ignored if null 
 /// which is the default
 // hierher make read only and use set/unset fct to enable/disable;
 // currently we'd have to reset this to null!
 ExactNonSingularFctPt& exact_non_singular_fct_pt()
  {
   return Exact_non_singular_fct_pt;
  }


protected:

 /// \short Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape(s,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian(s);
  }


 /// \short Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test_at_knot(const unsigned &ipt,
                                      Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape_at_knot(ipt,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian_at_knot(ipt);
  }


private:


 /// \short Add the element's contribution to its residual vector.
 /// flag=1(or 0): do (or don't) compute the contribution to the
 /// Jacobian as well. 
 void fill_in_generic_residual_contribution_poisson_flux(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag);
 
 ///The spatial dimension of the problem
 unsigned Dim;

 ///The index at which the unknown is stored at the nodes
 unsigned U_index_poisson;

 /// Pointer to element that stores the singular fcts etc.
 TemplateFreeScalableSingularityForPoissonElement* Poisson_sing_el_pt;

 /// Pointer to exact non singular fct (and its gradient) only used
 /// to validate the computation of the integral. Ignored if null 
 /// which is the default
 ExactNonSingularFctPt Exact_non_singular_fct_pt;
 
}; 

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//===========================================================================
/// Constructor, takes the pointer to the "bulk" element, the 
/// index of the fixed local coordinate and its value represented
/// by an integer (+/- 1), indicating that the face is located
/// at the max. or min. value of the "fixed" local coordinate
/// in the bulk element.
//===========================================================================
template<class ELEMENT>
PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
PoissonWithSingularityBoundaryIntegralFaceElement(FiniteElement* const &bulk_el_pt, 
                   const int &face_index) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  { 

   // Null out the fct pointer so integral is computed with
   // actual finite element representation of the non singular
   // part
   Exact_non_singular_fct_pt=0;


   // Let the bulk element build the FaceElement, i.e. setup the pointers 
   // to its nodes (by referring to the appropriate nodes in the bulk
   // element), etc.
   bulk_el_pt->build_face_element(face_index,this);

#ifdef PARANOID
   {
    //Check that the element is not a refineable 3d element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
    //If it's three-d
    if(elem_pt->dim()==3)
     {
      //Is it refineable
      RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>(elem_pt);
      if(ref_el_pt!=0)
       {
        if (this->has_hanging_nodes())
         {
          throw OomphLibError(
           "This flux element will not work correctly if nodes are hanging\n",
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         }
       }
     }
   }
#endif

   // Initialising the pointer to the singularity function
   this->Poisson_sing_el_pt = 0;
 
   // Extract the dimension of the problem from the dimension of 
   // the first node
   Dim = this->node_pt(0)->ndim();

   //Set up U_index_poisson. Initialise to zero, which probably won't change
   //in most cases, oh well, the price we pay for generality
   U_index_poisson = 0;

   //Cast to the appropriate PoissonEquation so that we can
   //find the index at which the variable is stored
   //We assume that the dimension of the full problem is the same
   //as the dimension of the node, if this is not the case you will have
   //to write custom elements, sorry
   switch(Dim)
    {
     //One dimensional problem
    case 1:
    {
     PoissonEquations<1>* eqn_pt = 
      dynamic_cast<PoissonEquations<1>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are one dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<1>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
     //Otherwise read out the value
     else
      {
       //Read the index from the (cast) bulk element
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;
    
    //Two dimensional problem
    case 2:
    {
     PoissonEquations<2>* eqn_pt = 
      dynamic_cast<PoissonEquations<2>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are two dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<2>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }
     else
      {
       //Read the index from the (cast) bulk element.
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;
    
    //Three dimensional problem
    case 3:
    {
     PoissonEquations<3>* eqn_pt = 
      dynamic_cast<PoissonEquations<3>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are three dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<3>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
       
      }
     else
      {
       //Read the index from the (cast) bulk element.
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;

    //Any other case is an error
    default:
     std::ostringstream error_stream; 
     error_stream <<  "Dimension of node is " << Dim 
                  << ". It should be 1,2, or 3!" << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     break;
    }
  }



//===========================================================================
/// Compute the element's residual vector and the (zero) Jacobian matrix.
//===========================================================================
template<class ELEMENT>
void PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
fill_in_generic_residual_contribution_poisson_flux(
 Vector<double> &residuals, DenseMatrix<double> &jacobian, 
 const unsigned& flag)
{

 // hierher populate when contribution is split up
 oomph_info << "This shouldn't be called at the moment\n";
 abort();
}

//===========================================================================
/// Calculate the contribution of the face element to the integral that
/// determines the amplitude hierher: this will be split into contributions
/// to the residual later on!
//===========================================================================
template<class ELEMENT>
double PoissonWithSingularityBoundaryIntegralFaceElement<ELEMENT>::
get_contribution_integral(std::ofstream &outfile)
{

 bool do_output=false;
 if (outfile.is_open())
  {
   do_output=true;
  }

 //Find out how many nodes there are
 const unsigned n_node = nnode();
  
 //Set up memory for the shape and test functions
 Shape psif(n_node), testf(n_node);
 
 //Set the value of Nintpt
 const unsigned n_intpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(Dim-1);
 
 // Saves result of integration
 double integral_result = 0.0;

 // Output stuff?
 if (do_output)
  {
   // hierher this won't work in 3D!
   outfile << "ZONE I=" << n_intpt << std::endl;
  }

 //Loop over the integration points
 //--------------------------------
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Find the shape and test functions and return the Jacobian
   //of the mapping
   double J = shape_and_test(s,psif,testf);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;

   // compute outer normal unit vector
   Vector<double> unit_normal(Dim);
   outer_unit_normal(s,unit_normal);

   // Get the gradient of u_fe and global coordinates
   Vector<double> flux(Dim);
   Vector<double> s_bulk(Dim);
   Vector<double> x(Dim); 

   for(unsigned i=0;i<Dim;i++)  
    { 
     x[i]=this->interpolated_x(s,i); 
    } 
   
   // Get gradient of singular function
   Vector<double> grad_u_sing(Dim);
   grad_u_sing = Poisson_sing_el_pt->gradient_of_unscaled_singular_fct(x);

   // Get the value of the singular function at our current location
   double sing_fct_value = Poisson_sing_el_pt->unscaled_singular_fct(x);

   // Get the local bulk coordinates    
   s_bulk=local_coordinate_in_bulk(s);
   ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(this->bulk_element_pt());

   // Get FE part of the solution
   double value_fe = bulk_el_pt->raw_interpolated_u_poisson(s_bulk);

   // Get the Gradient of the FE part of the solution
   bulk_el_pt->get_flux(s_bulk,flux);

   // Overwrite with exact version!
   double backed_up_value_fe=value_fe;
   Vector<double> backed_up_flux(flux);
   if (Exact_non_singular_fct_pt!=0)
    {
     Exact_non_singular_fct_pt(x,value_fe,flux);
    }

   // Output stuff?
   if (do_output)
    {
     for(unsigned i=0;i<Dim;i++)  
      { 
       outfile << x[i] << " ";
      }
     for(unsigned i=0;i<Dim;i++)  
      { 
       outfile << unit_normal[i] << " ";
      }
     outfile << value_fe << " ";
     for(unsigned i=0;i<Dim;i++)  
      { 
       outfile << flux[i] << " ";
      }
     outfile << sing_fct_value << " ";
     for(unsigned i=0;i<Dim;i++)  
      { 
       outfile << grad_u_sing[i] << " ";
      }
     double integrand=0.0;
     double fe_based_integrand=0.0;
     for(unsigned i=0;i<Dim;i++)
      {
       integrand += 
        unit_normal[i]*(flux[i]*sing_fct_value 
                        - grad_u_sing[i]*value_fe); 
       fe_based_integrand += 
        unit_normal[i]*(backed_up_flux[i]*sing_fct_value 
                        - grad_u_sing[i]*backed_up_value_fe); 
      }
     double diff=fe_based_integrand-integrand;
     outfile << integrand << " " 
             << fe_based_integrand << " "
             << diff << " ";
     outfile << std::endl;
    } 

    // Now we compute the contibution of the integral
    for(unsigned i=0;i<Dim;i++)
    {
     integral_result += unit_normal[i]*W*(flux[i]*sing_fct_value 
       - grad_u_sing[i]*value_fe); 
    } 
  }
 return integral_result;
}

//======================================================================
/// \short Class for elements that handle singularities
/// in Poisson/Laplace equations. Templated by bulk element within
/// which we impose regularity on the FE solution by insisting that
/// the slope of the solution at a specified local coordinate, and in
/// in a specified direction is zero. Nodal values of that element
/// become external data for the current element whose equation
/// (zero slope of the FE solution, as discussed) determines the 
/// amplitude of the singular function.
//======================================================================
template<class BULK_ELEMENT> 
class ScalableSingularityForPoissonElement : 
  public virtual TemplateFreeScalableSingularityForPoissonElement
 {
 
   public:

  /// Constructor
   ScalableSingularityForPoissonElement() :
  Bulk_element_pt(0), Face_element_mesh_pt(0), 
   Impose_singular_fct_amplitude(false)
   {
   }



  /// Set pointer to mesh containing the FaceElements (and flush
  /// the previuos ones first!)
  void set_mesh_of_face_elements(Mesh* const& face_mesh_pt)
   {
    Face_element_mesh_pt=face_mesh_pt;
    flush_external_data();
    unsigned nel=face_mesh_pt->nelement();
    oomph_info << "nmber of face elements used to compute C: "
               << nel << std::endl;
    for (unsigned e=0;e<nel;e++)
     {
      FiniteElement* el_pt=
       dynamic_cast<PoissonWithSingularityBoundaryIntegralFaceElement<BULK_ELEMENT>*>(
        face_mesh_pt->element_pt(e))->bulk_element_pt();
      unsigned nnod=el_pt->nnode();
      for (unsigned j=0;j<nnod;j++)
       {
        add_external_data(el_pt->node_pt(j));
       }
     }
   }



  /// Call this to bypass the correct computation of the
  /// residual and replace it by r_c = C-ampl
  void impose_singular_fct_amplitude(double const& ampl)
  {
   Impose_singular_fct_amplitude=true;
   Imposed_amplitude = ampl;
  } 


  /// Reset to compute r_c properly via integral
  void dont_impose_singular_fct_amplitude()
  {
   Impose_singular_fct_amplitude=false;
  } 
  

  /// Add the element's contribution to its residual vector
  inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_poisson_sing_fct(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }
  
   private:

  /// Add the element's contribution to its residual vector
  inline void fill_in_generic_residual_contribution_poisson_sing_fct(
   Vector<double> &residuals,
   DenseMatrix<double> &jacobian,
   const unsigned& flag)
  {

   if (ndof()==0) return;

#ifdef PARANOID
   // hierher paranoid check null pointers and zero sized vectors    
#endif
   

   // Get eqn number of residual that determines C
   int c_local_eqn=internal_local_eqn(0,0);
   if (c_local_eqn>=0)
    {

     // Bypass actual computation?
     if (Impose_singular_fct_amplitude)
      {
       residuals[c_local_eqn] = this->amplitude_of_singular_fct() - 
        Imposed_amplitude;
      }
     // Do it properly
     else
      {
       unsigned n_element = Face_element_mesh_pt->nelement();
       for(unsigned e = 0; e<n_element; e++)
        {
         residuals[c_local_eqn]+= 
          dynamic_cast<PoissonWithSingularityBoundaryIntegralFaceElement
          <BULK_ELEMENT>*>(Face_element_mesh_pt->finite_element_pt(e))
          ->get_contribution_integral();
        }
      }
    }


   /* oomph_info  << "Residual r_C : "  */
   /*             << residuals[c_local_eqn] */
   /*             << std::endl; */
  }

   private:
  
  
  /// Pointer to bulk element where FE solution is regularised
  BULK_ELEMENT* Bulk_element_pt;

  /// Pointer to mesh of face elements that contribute to the surface
  /// integral that determines the amplitude of the unkown function
  Mesh* Face_element_mesh_pt;
  
  /// Imposed amplitude (only used if Impose_singular_fct_amplitude=true)
  double Imposed_amplitude;  

  /// \short Boolean to bypass the correct computation of the
  /// residual and replace it by r_c = C-ampl
  bool Impose_singular_fct_amplitude;

 };




////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

// hierher really need to tidy this up! Should only need one class 
// for T and Q
//
//====================================================================
/// New class. Mainly overloads output-related functions to add
/// "singular function" (which is assumed to satisfy the Laplace
/// equation; therefore no change to the governing (bulk) equations) 
/// to the FE solution. 
//====================================================================
template<unsigned DIM, unsigned NNODE_1D>
class MyTPoissonElement : public virtual TPoissonElement<DIM,NNODE_1D>
{
 

public:

 typedef void (*ExactNonSingularFctPt)
  (const Vector<double>& x, double& u, Vector<double>& grad_u);

 ExactNonSingularFctPt& exact_non_singular_fct_pt()
  {
   return Exact_non_singular_fct_pt;
  }  

 /// Constructor
 MyTPoissonElement() : Poisson_sing_el_pt(0),Exact_non_singular_fct_pt(0)
  {
  }

 /// \short Return FE representation of function value u_poisson(s) 
 /// plus scaled singular fct (if provided) at local coordinate s
 inline double interpolated_u_poisson(const Vector<double> &s) const
  {
   double u_fe=TPoissonElement<DIM,NNODE_1D>::interpolated_u_poisson(s);

   // hierher renable
    if (Poisson_sing_el_pt!=0) 
     {      
      Vector<double> x(DIM); 
      for(unsigned i=0;i<DIM;i++)  
       { 
        x[i]=this->interpolated_x(s,i); 
       } 
      u_fe+=Poisson_sing_el_pt->singular_fct(x); 
     } 
   // hierher end renable
   
   return u_fe;
  } 

 /// \short Return FE representation of function value u_fe
 inline double raw_interpolated_u_poisson(const Vector<double> &s) const
  {
    double u_fe=TPoissonElement<DIM,NNODE_1D>::interpolated_u_poisson(s);
    return u_fe;
  } 

 /// Output with various contributions
 void  output_with_various_contributions(std::ostream &outfile, 
                                         const unsigned &nplot)
  {
   //Vector of local coordinates
   Vector<double> s(DIM);
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);
     
     Vector<double> x(DIM);
     for(unsigned i=0;i<DIM;i++) 
      {
       x[i]=this->interpolated_x(s,i);
       outfile << x[i] << " ";
      }
     double u_sing=0.0;
     // hierher renable
     if (Poisson_sing_el_pt!=0) 
       { 
        u_sing=Poisson_sing_el_pt->singular_fct(x); 
       }

     double u_non_sing = 0.0;

     Vector<double> flux(DIM);

     // Overwrite with exact version!
     if (Exact_non_singular_fct_pt!=0)
      {
       Exact_non_singular_fct_pt(x,u_non_sing,flux);
      }  
     // hierher end renable
     outfile << this->interpolated_u_poisson(s) << " "
             << TPoissonElement<DIM,NNODE_1D>::interpolated_u_poisson(s) << " "
             << u_sing << " " << u_non_sing << " "
             << std::endl;   
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   this->write_tecplot_zone_footer(outfile,nplot);
   
  }

 
 /// Pointer to element that stores singular fct 
  TemplateFreeScalableSingularityForPoissonElement*& poisson_sing_el_pt() 
  { 
    return Poisson_sing_el_pt; 
  } 

private:


  /// Pointer to element that stores singular fct 
  TemplateFreeScalableSingularityForPoissonElement* Poisson_sing_el_pt; 
 
  /// Pointer to exact non-singular fct (only for post-processing!)
  ExactNonSingularFctPt Exact_non_singular_fct_pt;

};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////




//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the MyTPoissonElement elements: The spatial 
/// dimension of the face elements is one lower than that of the
/// bulk element but they have the same number of points
/// along their 1D edges.
//=======================================================================
template<unsigned DIM, unsigned NNODE_1D>
class FaceGeometry<MyTPoissonElement<DIM,NNODE_1D> >: 
 public virtual TElement<DIM-1,NNODE_1D>
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : TElement<DIM-1,NNODE_1D>() {}

};

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=======================================================================
/// Face geometry for the 1D MyTPoissonElement elements: Point elements
//=======================================================================
template<unsigned NNODE_1D>
class FaceGeometry<MyTPoissonElement<1,NNODE_1D> >: 
 public virtual PointElement
{

  public:
 
 /// \short Constructor: Call the constructor for the
 /// appropriate lower-dimensional TElement
 FaceGeometry() : PointElement() {} 

};



///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////


//======================================================================
/// \short A class for elements that imposes Dirichlet boundary 
/// conditions on complete solution (such that u_fe + C u_sing = u_bc) using a
/// Lagrange multiplier. Thus the element introduce an additional
/// unknown at the nodes it's attached to. C and u_sing are specified
/// via a ScalableSingularityForPoissonElement.
//======================================================================
 template <class ELEMENT>
  class PoissonWithSingularityBCFaceElement : 
  public virtual FaceGeometry<ELEMENT>, 
  public virtual FaceElement 
 {
 
   public:

  /// \short Constructor, takes the pointer to the "bulk" element and the 
  /// index of the face to which the element is attached. Optional final
  /// arg is the identifier for the additional unknowns multiplier
  PoissonWithSingularityBCFaceElement(FiniteElement* const &bulk_el_pt, 
                                    const int& face_index,
                                    const unsigned &id=0); 
  
  ///\short  Broken empty constructor
  PoissonWithSingularityBCFaceElement()
   {
    throw OomphLibError(
     "Don't call empty constructor for PoissonWithSingularityBCFaceElement",
     OOMPH_CURRENT_FUNCTION,
     OOMPH_EXCEPTION_LOCATION);
   }
  
  /// Broken copy constructor
  PoissonWithSingularityBCFaceElement(
   const PoissonWithSingularityBCFaceElement& dummy) 
   { 
    BrokenCopy::broken_copy("PoissonWithSingularityBCFaceElement");
   } 
  
  /// Broken assignment operator
  void operator=(const PoissonWithSingularityBCFaceElement&) 
   {
    BrokenCopy::broken_assign("PoissonWithSingularityBCFaceElement");
   }
  
  /// \short Specify the value of nodal zeta from the face geometry
  /// The "global" intrinsic coordinate of the element when
  /// viewed as part of a geometric object should be given by
  /// the FaceElement representation, by default (needed to break
  /// indeterminacy if bulk element is SolidElement)
  double zeta_nodal(const unsigned &n, const unsigned &k,           
                    const unsigned &i) const 
  {return FaceElement::zeta_nodal(n,k,i);}     


   /// Pointer to element that handles singular fct
   ScalableSingularityForPoissonElement<ELEMENT>* poisson_sing_el_pt() const
   {
    return Poisson_sing_el_pt;
   }

   /// \short Set pointer to element that stores singular fct. Data that stores
   /// the amplitude of the singular fct and its index is retrieved from
   /// that element so the Data can be used as external Data in this
   /// element.
   void set_poisson_sing_el_pt(ScalableSingularityForPoissonElement<ELEMENT>* 
                               poisson_sing_el_pt) 
   {
    Poisson_sing_el_pt=poisson_sing_el_pt;
    C_external_data_index=add_external_data(
     poisson_sing_el_pt->data_that_stores_amplitude_of_singular_fct());
    C_external_data_value_index=
     poisson_sing_el_pt->index_of_value_that_stores_amplitude_of_singular_fct();
   }


  /// Add the element's contribution to its residual vector
  inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_poisson_sing(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }


  /// \short Add the element's contribution to its residual vector and its
  /// Jacobian matrix
  inline void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                               DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_poisson_sing(residuals,jacobian,1);
  }

  /// Output function
  void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   output(outfile,n_plot);
  }

  /// \short Output function
  void output(std::ostream &outfile, const unsigned &nplot)
  {
   //oomph_info << "hierher need to update output fct" << std::endl;
       //Vector of local coordinates
   Vector<double> s(Dim-1);
   
   // Tecplot header info
   outfile << this->tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=this->nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     // Get local coordinates of plot point
     this->get_s_plot(iplot,nplot,s);
     
     Vector<double> x(Dim);
     for(unsigned i=0;i<Dim;i++) 
      {
       x[i]=this->interpolated_x(s,i);
       outfile << x[i] << " ";
      }
     outfile << endl;
    }
   return;
  }

  /// \short Provide nodal values of desired boundary values.
  /// They're imposed by Lagrange multipliers.
  void set_nodal_boundary_values(const Vector<double>& nodal_boundary_value)
  {
#ifdef PARANOID
   if (nodal_boundary_value.size()!=nnode())
    {
     std::stringstream error;
     error << "nodel_boundary_value is a vector of size " 
           << nodal_boundary_value.size() 
           << " but should have the same size as the number of nodes, "
           << nnode();
     throw OomphLibError(error.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
    }
#endif
   Nodal_boundary_value=nodal_boundary_value;
  }

  /// Pin Lagrange multiplier at specified local node
  void pin_lagrange_multiplier_at_specified_local_node(const unsigned& j)
  {
   node_pt(j)->pin(Lambda_index[j]);
  }

  /// Unpin FE part of the solution at specified local node
  void unpin_u_fe_at_specified_local_node(const unsigned&j )
  {   
   node_pt(j)->unpin(U_index_poisson);
  }

  /// C-style output function -- forward to broken version in FiniteElement
  /// until somebody decides what exactly they want to plot here...
  void output(FILE* file_pt) {FiniteElement::output(file_pt);}

  /// \short C-style output function -- forward to broken version in 
  /// FiniteElement until somebody decides what exactly they want to plot 
  /// here...
  void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}


   protected:

  /// \short Function to compute the shape and test functions and to return 
  /// the Jacobian of mapping between local and global (Eulerian)
  /// coordinates
  inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
   const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape(s,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian(s);
  }


  /// \short Function to compute the shape and test functions and to return 
  /// the Jacobian of mapping between local and global (Eulerian)
  /// coordinates
  inline double shape_and_test_at_knot(const unsigned &ipt,
                                       Shape &psi, Shape &test)
   const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape_at_knot(ipt,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian_at_knot(ipt);
  }


   private:


  /// \short Add the element's contribution to its residual vector.
  /// flag=1(or 0): do (or don't) compute the contribution to the
  /// Jacobian as well. 
  void fill_in_generic_residual_contribution_poisson_sing(
   Vector<double> &residuals, DenseMatrix<double> &jacobian, 
   const unsigned& flag);
 
 
  ///The spatial dimension of the problem
  unsigned Dim;

  ///The index at which the Poisson unknown is stored at the nodes
  unsigned U_index_poisson;

  /// \short The index at which the Lagrange multiplier that enforces
  /// the Dirichlet BC is stored at the nodes
  Vector<unsigned> Lambda_index;

  /// Desired boundary values at nodes
  Vector<double> Nodal_boundary_value;

  /// \short Index of external Data that stores the value of the amplitude of
  /// the singular function
  unsigned C_external_data_index;
  
  /// \short Index of value (within external Data) that stores the
  /// value of the amplitude of the singular function
  unsigned C_external_data_value_index;
  
  /// \short Pointer to element that stores pointer to singular fct 
  /// (and its gradients etc.) as well as amplitude
  ScalableSingularityForPoissonElement<ELEMENT>* Poisson_sing_el_pt;

 }; 

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//===========================================================================
/// Constructor, takes the pointer to the "bulk" element, the 
/// index of the fixed local coordinate and its value represented
/// by an integer indicating which face we're on.
/// Optional final arg is the identifier for the new values created
/// by this face element
//===========================================================================
 template<class ELEMENT>
  PoissonWithSingularityBCFaceElement<ELEMENT>::
  PoissonWithSingularityBCFaceElement(FiniteElement* const &bulk_el_pt, 
                                    const int &face_index, 
                                    const unsigned& id) : 
 FaceGeometry<ELEMENT>(), FaceElement()
  { 

   // Initialise
   Poisson_sing_el_pt=0;

   // Let the bulk element build the FaceElement, i.e. setup the pointers 
   // to its nodes (by referring to the appropriate nodes in the bulk
   // element), etc.
   bulk_el_pt->build_face_element(face_index,this);
 
#ifdef PARANOID
   {
    //Check that the element is not a refineable 3d element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
    //If it's three-d
    if(elem_pt->dim()==3)
     {
      //Is it refineable
      RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>(elem_pt);
      if(ref_el_pt!=0)
       {
        if (this->has_hanging_nodes())
         {
          throw OomphLibError(
           "This face element will not work correctly if nodes are hanging\n",
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         }
       }
     }
   }
#endif   

   // Extract the dimension of the problem from the dimension of 
   // the first node
   Dim = this->node_pt(0)->ndim();

   //Set up U_index_poisson. Initialise to zero, which probably won't change
   //in most cases, oh well, the price we pay for generality
   U_index_poisson = 0;

   //Cast to the appropriate PoissonEquation so that we can
   //find the index at which the variable is stored
   //We assume that the dimension of the full problem is the same
   //as the dimension of the node, if this is not the case you will have
   //to write custom elements, sorry
   switch(Dim)
    {
     //One dimensional problem
    case 1:
    {
     PoissonEquations<1>* eqn_pt = 
      dynamic_cast<PoissonEquations<1>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are one dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<1>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
     //Otherwise read out the value
     else
      {
       //Read the index from the (cast) bulk element
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;
    
    //Two dimensional problem
    case 2:
    {
     PoissonEquations<2>* eqn_pt = 
      dynamic_cast<PoissonEquations<2>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are two dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<2>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
     else
      {
       //Read the index from the (cast) bulk element.
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;
    
    //Three dimensional problem
    case 3:
    {
     PoissonEquations<3>* eqn_pt = 
      dynamic_cast<PoissonEquations<3>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are three dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<3>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
       
      }
     else
      {
       //Read the index from the (cast) bulk element.
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;

    //Any other case is an error
    default:
     std::ostringstream error_stream; 
     error_stream <<  "Dimension of node is " << Dim 
                  << ". It should be 1,2, or 3!" << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     break;
    }



   // Where is the extra dof representing the Lagrange multiplier stored?
   // Initially store number of values stored right now
   unsigned nnod=nnode();
   Lambda_index.resize(nnod);
   for (unsigned j=0;j<nnod;j++)
    {
     Lambda_index[j]=node_pt(j)->nvalue();
    }

   // Make space for one Lagrange multiplier
   Vector<unsigned> n_additional_values(nnod,1);
   this->add_additional_values(n_additional_values,id);


   // Now check if we've added a new value. If so, that's
   // the Lagrange multiplier; it not, it was already stored
   // there so the actual index is one less
   for (unsigned j=0;j<nnod;j++)
    {
     if (Lambda_index[j]==node_pt(j)->nvalue())
      {
       Lambda_index[j]--;
      }
    }

  }


//===========================================================================
/// Compute the element's residual vector and the Jacobian matrix.
//===========================================================================
 template<class ELEMENT>
  void PoissonWithSingularityBCFaceElement<ELEMENT>::
  fill_in_generic_residual_contribution_poisson_sing(
   Vector<double> &residuals, DenseMatrix<double> &jacobian, 
   const unsigned& flag) 
  {

   /* oomph_info << "in here: Poisson_sing_el_pt = "  */
   /*            << Poisson_sing_el_pt << std::endl; */

   // hierher this is all for 1D; keep around until DIM-dimensional 
   // version works
   if (this->dim()==0)
    {
     
     // Compute various quantities
     Vector<double> x(1);
     x[0]=node_pt(0)->x(0);
     double u_sing=Poisson_sing_el_pt->singular_fct(x);
     double lambda=node_pt(0)->value(Lambda_index[0]);
     double u_fe=node_pt(0)->value(U_index_poisson);
     double u_bc=Nodal_boundary_value[0];
     
     // Various local equation numbers
     int local_eqn_lagr=nodal_local_eqn(0,Lambda_index[0]);
     int local_eqn_u_fe=nodal_local_eqn(0,U_index_poisson);
     
     // Get flux (=gradient) vector in the bulk element. We're
     // making the fe solution regular by setting this to zero!
     Vector<double> s_flux(1); // hierher local coordinate where flux is to be
     // evaluated should get passed in, together with normal
     s_flux[0]=-1.0;
     Vector<double> flux(1);
     ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());
     bulk_el_pt->get_flux(s_flux,flux);
     
     // Derivatives of flux (in bulk element) w.r.t. to nodal unknowns
     // in bulk element
     unsigned nnod_bulk=bulk_el_pt->nnode();
     Vector<Vector<double> > dflux_dnodal_u(1); // hierher loop properly
     dflux_dnodal_u[0].resize(nnod_bulk);
     bulk_el_pt->get_dflux_dnodal_u(s_flux,dflux_dnodal_u);
     
     
#ifdef PARANOID
     // Lagrange multiplier active but u_fe pinned won't work!
     if ( (local_eqn_lagr>=0) && (local_eqn_u_fe<0) )
      {
       throw OomphLibError(
        "Lagrange multiplier active but u_fe pinned won't work!",
        OOMPH_CURRENT_FUNCTION,
        OOMPH_EXCEPTION_LOCATION);
      }
#endif

     // Lagrange multiplier for BC residual: It's determined by enforcing
     // that u_fe + C u_sing = u_bc
     if (local_eqn_lagr>=0)
      {
       residuals[local_eqn_lagr]+=((u_fe+u_sing)-u_bc); 
       if (flag==1)
        {
         if (local_eqn_u_fe>=0) jacobian(local_eqn_lagr,local_eqn_u_fe)=1.0;
        }
      }
     
     // Contribution to bulk equation: Lagrange multiplier
     if (local_eqn_u_fe>=0)
      {
       residuals[local_eqn_u_fe]+=lambda;
       if (flag==1)
        {
         if (local_eqn_lagr>=0) jacobian(local_eqn_u_fe,local_eqn_lagr)+=1.0;
        }
      }
    }
   else
    {
     //Find out how many nodes there are
     const unsigned n_node = nnode();
     
     //Set up memory for the shape and test functions
     Shape psi(n_node), test(n_node);
     
     //Set the value of Nintpt
     const unsigned n_intpt = integral_pt()->nweight();
     
     //Set the Vector to hold local coordinates
     Vector<double> s(Dim-1);
     
     
     //Loop over the integration points
     //--------------------------------
     for(unsigned ipt=0;ipt<n_intpt;ipt++)
      {
       
       //Assign values of s
       for(unsigned i=0;i<(Dim-1);i++) 
        {
         s[i] = integral_pt()->knot(ipt,i);
        }
       
       //Get the integral weight
       double w = integral_pt()->weight(ipt);
       
       //Find the shape and test functions and return the Jacobian
       //of the mapping
       double J = shape_and_test(s,psi,test);
       
       //Premultiply the weights and the Jacobian
       double W = w*J;
       
       //Calculate stuff at integration point
       Vector<double> interpolated_x(Dim,0.0);
       double u_fe=0.0;
       double u_bc=0.0;
       double lambda=0.0;
       for(unsigned l=0;l<n_node;l++) 
        {         
         u_fe+=this->nodal_value(l,U_index_poisson)*psi[l];
         u_bc+=Nodal_boundary_value[l]*psi[l];
         lambda+=this->nodal_value(l,Lambda_index[l])*psi[l];
         for(unsigned i=0;i<Dim;i++)
          {
           interpolated_x[i] += this->nodal_position(l,i)*psi[l];
          }
        }
       
       // Stuff related to singular fct
       double u_sing=0.0;
       double u_sing_unscaled=0.0;
       if (Poisson_sing_el_pt!=0)
        {
         u_sing=Poisson_sing_el_pt->singular_fct(interpolated_x);  
         u_sing_unscaled=Poisson_sing_el_pt->
          unscaled_singular_fct(interpolated_x);
        }
       
       //Now add to the appropriate equations
       
       //Loop over the test functions
       for(unsigned l=0;l<n_node;l++)
        {
         int local_eqn_lagr=nodal_local_eqn(l,Lambda_index[l]);
         int local_eqn_u_fe=nodal_local_eqn(l,U_index_poisson);
         int local_eqn_c=-1;         
         if (Poisson_sing_el_pt!=0)
          {
           local_eqn_c=external_local_eqn(C_external_data_index,
                                          C_external_data_value_index);
          }

#ifdef PARANOID
         // Lagrange multiplier active but u_fe pinned won't work!
         if ( (local_eqn_lagr>=0) && (local_eqn_u_fe<0) )
          {
           throw OomphLibError(
            "Lagrange multiplier active but u_fe pinned won't work!",
            OOMPH_CURRENT_FUNCTION,
            OOMPH_EXCEPTION_LOCATION);
          }
#endif
         // Lagrange multiplier for BC residual: It's determined by enforcing
         // that u_fe + C u_sing = u_bc
         if(local_eqn_lagr >= 0)
          {
           residuals[local_eqn_lagr] += ((u_fe+u_sing)-u_bc)*test[l]*W;
           
           // Jacobian?
           if (flag==1)
            {
             for(unsigned l2=0;l2<n_node;l2++)
              {
               int local_unknown_u_fe=nodal_local_eqn(l2,U_index_poisson);
               if (local_unknown_u_fe>=0)
                {
                 jacobian(local_eqn_lagr,local_unknown_u_fe)+=psi[l2]*test[l]*W;
                }
              }
             
             // Deriv. w.r.t. amplitude is simply the unscaled singular fct.
             if (local_eqn_c>=0)
              {
               jacobian(local_eqn_lagr,local_eqn_c)+=u_sing_unscaled*test[l]*W;
              }
            }
          }
         
         // Contribution of Lagrange multiplier to bulk eqn:
         if (local_eqn_u_fe>=0)
          {
           residuals[local_eqn_u_fe]+=lambda*test[l]*W;
           if (flag==1)
            {
             for(unsigned l2=0;l2<n_node;l2++)
              {
               int local_unknown_lambda=nodal_local_eqn(l2,Lambda_index[l2]);
               if (local_unknown_lambda>=0)
                {
                 jacobian(local_eqn_u_fe,local_unknown_lambda)+=
                  psi[l2]*test[l]*W;
                }
              }
            }
          }

        }
      }
    }
  }


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////



//======================================================================
/// \short A class for elements that allow the imposition of an 
/// applied flux on the boundaries of Poisson elements.
/// The element geometry is obtained from the  FaceGeometry<ELEMENT> 
/// policy class.
//======================================================================
template <class ELEMENT>
class PoissonWithSingularityFluxElement : public virtual FaceGeometry<ELEMENT>, 
public virtual FaceElement
{
 
public:

 /// \short Function pointer to the prescribed-flux function fct(x,f(x)) -- 
 /// x is a Vector! 
 typedef void (*PoissonPrescribedFluxFctPt)
  (const Vector<double>& x, double& flux);

 /// \short Constructor, takes the pointer to the "bulk" element and the 
 /// index of the face to which the element is attached.
 PoissonWithSingularityFluxElement(FiniteElement* const &bulk_el_pt, 
                    const int& face_index); 

 ///\short  Broken empty constructor
 PoissonWithSingularityFluxElement()
  {
   throw OomphLibError(
    "Don't call empty constructor for PoissonWithSingularityFluxElement",
    OOMPH_CURRENT_FUNCTION,
    OOMPH_EXCEPTION_LOCATION);
  }

 /// Broken copy constructor
 PoissonWithSingularityFluxElement(const PoissonWithSingularityFluxElement& dummy) 
  { 
   BrokenCopy::broken_copy("PoissonWithSingularityFluxElement");
  } 
 
 /// Broken assignment operator
 void operator=(const PoissonWithSingularityFluxElement&) 
  {
   BrokenCopy::broken_assign("PoissonWithSingularityFluxElement");
  }

 /// \short Specify the value of nodal zeta from the face geometry
 /// The "global" intrinsic coordinate of the element when
 /// viewed as part of a geometric object should be given by
 /// the FaceElement representation, by default (needed to break
 /// indeterminacy if bulk element is SolidElement)
 double zeta_nodal(const unsigned &n, const unsigned &k,           
                          const unsigned &i) const 
  {return FaceElement::zeta_nodal(n,k,i);}     

 /// Access function for the prescribed-flux function pointer
 PoissonPrescribedFluxFctPt& flux_fct_pt() {return Flux_fct_pt;}


 /// Add the element's contribution to its residual vector
 inline void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_poisson_flux(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }



 // hierher forced this to be done by finite differencing (for now)
 // because the amplitude of the singular function does provide
 // a contribution to the Jacobian
 /* /// \short Add the element's contribution to its residual vector and its  */
 /* /// Jacobian matrix */
 /* inline void fill_in_contribution_to_jacobian(Vector<double> &residuals, */
 /*                                          DenseMatrix<double> &jacobian) */
 /*  { */
 /*   //Call the generic routine with the flag set to 1 */
 /*   fill_in_generic_residual_contribution_poisson_flux(residuals,jacobian,1); */
 /*  } */

 /// Output function
 void output(std::ostream &outfile)
 {
  const unsigned n_plot=5;
  output(outfile,n_plot);
 }

 /// \short Output function
 void output(std::ostream &outfile, const unsigned &nplot)
  {
   
   // Dimension of element 
   unsigned el_dim=dim();

   //Vector of local coordinates
   Vector<double> s(el_dim);
   
   // Tecplot header info
   outfile << tecplot_zone_string(nplot);
   
   // Loop over plot points
   unsigned num_plot_points=nplot_points(nplot);
   for (unsigned iplot=0;iplot<num_plot_points;iplot++)
    {
     
     // Get local coordinates of plot point
     get_s_plot(iplot,nplot,s);
     
     Vector<double> x(el_dim+1);
     for(unsigned i=0;i<el_dim+1;i++) 
      {
       x[i]=interpolated_x(s,i);
       outfile << x[i] << " ";
      }

     // Compute outer unit normal at the specified local coordinate
     Vector<double> unit_normal(Dim);
     outer_unit_normal(s,unit_normal);
     
     // Get FE flux
     ELEMENT* bulk_el_pt=dynamic_cast<ELEMENT*>(bulk_element_pt());
     Vector<double> s_bulk(Dim);
     s_bulk=local_coordinate_in_bulk(s);
     Vector<double> fe_flux(Dim);
     bulk_el_pt->get_flux(s_bulk,fe_flux);
     
     
     // Get gradient of singular fct (incl. amplitude)
     Vector<double> grad_sing(Dim,0.0);
     if (Poisson_sing_el_pt!=0)
      {
       grad_sing=Poisson_sing_el_pt->gradient_of_singular_fct(x);
      }

     // Get actual flux 
     double actual_flux=0.0;
     for (unsigned i=0;i<Dim;i++)
      {
       actual_flux+=unit_normal[i]*(fe_flux[i]+grad_sing[i]);
      }
     
     double imposed_flux=0.0;
     get_flux(x,imposed_flux);
     outfile << imposed_flux << " " 
             << actual_flux << std::endl;   
    }
   
   // Write tecplot footer (e.g. FE connectivity lists)
   write_tecplot_zone_footer(outfile,nplot);

  }


 /// C-style output function -- forward to broken version in FiniteElement
 /// until somebody decides what exactly they want to plot here...
 void output(FILE* file_pt) {FiniteElement::output(file_pt);}

 /// \short C-style output function -- forward to broken version in 
 /// FiniteElement until somebody decides what exactly they want to plot 
 /// here...
 void output(FILE* file_pt, const unsigned &n_plot)
  {FiniteElement::output(file_pt,n_plot);}


 /// \short Set pointer to element that stores singular fct. Data that stores
 /// the amplitude of the singular fct and its index is retrieved from
 /// that element so the Data can be used as external Data in this
 /// element.
 void set_poisson_sing_el_pt(ScalableSingularityForPoissonElement<ELEMENT>* 
                             poisson_sing_el_pt) 
 {
  Poisson_sing_el_pt=poisson_sing_el_pt;
  C_external_data_index=add_external_data(
   poisson_sing_el_pt->data_that_stores_amplitude_of_singular_fct());
  C_external_data_value_index=
   poisson_sing_el_pt->index_of_value_that_stores_amplitude_of_singular_fct();
   }
 

protected:

 /// \short Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test(const Vector<double> &s, Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape(s,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian(s);
  }


 /// \short Function to compute the shape and test functions and to return 
 /// the Jacobian of mapping between local and global (Eulerian)
 /// coordinates
 inline double shape_and_test_at_knot(const unsigned &ipt,
                                      Shape &psi, Shape &test)
  const
  {
   //Find number of nodes
   unsigned n_node = nnode();

   //Get the shape functions
   shape_at_knot(ipt,psi);

   //Set the test functions to be the same as the shape functions
   for(unsigned i=0;i<n_node;i++) {test[i] = psi[i];}

   //Return the value of the jacobian
   return J_eulerian_at_knot(ipt);
  }


 /// Function to calculate the prescribed flux at a given spatial
 /// position
 void get_flux(const Vector<double>& x, double& flux)
  {
   //If the function pointer is zero return zero
   if(Flux_fct_pt == 0)
    {
     flux=0.0;
    }
   //Otherwise call the function
   else
    {
     (*Flux_fct_pt)(x,flux);
    }
  }
 
 /// Pointer to element that handles singular fct
 ScalableSingularityForPoissonElement<ELEMENT>* poisson_sing_el_pt() const
  {
   return Poisson_sing_el_pt;
  }
 
 
  private:


 /// \short Add the element's contribution to its residual vector.
 /// flag=1(or 0): do (or don't) compute the contribution to the
 /// Jacobian as well. 
 void fill_in_generic_residual_contribution_poisson_flux(
  Vector<double> &residuals, DenseMatrix<double> &jacobian, 
  const unsigned& flag);
 
 
 /// Function pointer to the (global) prescribed-flux function
 PoissonPrescribedFluxFctPt Flux_fct_pt;

 ///The spatial dimension of the problem
 unsigned Dim;

 ///The index at which the unknown is stored at the nodes
 unsigned U_index_poisson;
 
 /// \short Index of external Data that stores the value of the amplitude of
 /// the singular function
 unsigned C_external_data_index;
 
 /// \short Index of value (within external Data) that stores the
 /// value of the amplitude of the singular function
 unsigned C_external_data_value_index;
 
 /// \short Pointer to element that stores pointer to singular fct 
 /// (and its gradients etc.) as well as amplitude
 ScalableSingularityForPoissonElement<ELEMENT>* Poisson_sing_el_pt;
 
}; 

//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////



//===========================================================================
/// Constructor, takes the pointer to the "bulk" element, the 
/// index of the fixed local coordinate and its value represented
/// by an integer (+/- 1), indicating that the face is located
/// at the max. or min. value of the "fixed" local coordinate
/// in the bulk element.
//===========================================================================
template<class ELEMENT>
PoissonWithSingularityFluxElement<ELEMENT>::
PoissonWithSingularityFluxElement(FiniteElement* const &bulk_el_pt, 
                   const int &face_index) : 
  FaceGeometry<ELEMENT>(), FaceElement()
  {  
   // Initialise
   Poisson_sing_el_pt=0;

   // Let the bulk element build the FaceElement, i.e. setup the pointers 
   // to its nodes (by referring to the appropriate nodes in the bulk
   // element), etc.
   bulk_el_pt->build_face_element(face_index,this);
 
#ifdef PARANOID
   {
    //Check that the element is not a refineable 3d element
    ELEMENT* elem_pt = dynamic_cast<ELEMENT*>(bulk_el_pt);
    //If it's three-d
    if(elem_pt->dim()==3)
     {
      //Is it refineable
      RefineableElement* ref_el_pt=dynamic_cast<RefineableElement*>(elem_pt);
      if(ref_el_pt!=0)
       {
        if (this->has_hanging_nodes())
         {
          throw OomphLibError(
           "This flux element will not work correctly if nodes are hanging\n",
           OOMPH_CURRENT_FUNCTION,
           OOMPH_EXCEPTION_LOCATION);
         }
       }
     }
   }
#endif   

   // Initialise the prescribed-flux function pointer to zero
   Flux_fct_pt = 0;

   // Extract the dimension of the problem from the dimension of 
   // the first node
   Dim = this->node_pt(0)->ndim();

   //Set up U_index_poisson. Initialise to zero, which probably won't change
   //in most cases, oh well, the price we pay for generality
   U_index_poisson = 0;

   //Cast to the appropriate PoissonEquation so that we can
   //find the index at which the variable is stored
   //We assume that the dimension of the full problem is the same
   //as the dimension of the node, if this is not the case you will have
   //to write custom elements, sorry
   switch(Dim)
    {
     //One dimensional problem
    case 1:
    {
     PoissonEquations<1>* eqn_pt = 
      dynamic_cast<PoissonEquations<1>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are one dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<1>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                           OOMPH_EXCEPTION_LOCATION);
      }
     //Otherwise read out the value
     else
      {
       //Read the index from the (cast) bulk element
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;
    
    //Two dimensional problem
    case 2:
    {
     PoissonEquations<2>* eqn_pt = 
      dynamic_cast<PoissonEquations<2>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are two dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<2>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
      }
     else
      {
       //Read the index from the (cast) bulk element.
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;
    
    //Three dimensional problem
    case 3:
    {
     PoissonEquations<3>* eqn_pt = 
      dynamic_cast<PoissonEquations<3>*>(bulk_el_pt);
     //If the cast has failed die
     if(eqn_pt==0)
      {
       std::string error_string =
        "Bulk element must inherit from PoissonEquations.";
       error_string += 
        "Nodes are three dimensional, but cannot cast the bulk element to\n";
       error_string += "PoissonEquations<3>\n.";
       error_string += 
        "If you desire this functionality, you must implement it yourself\n";
       
       throw OomphLibError(error_string,
                           OOMPH_CURRENT_FUNCTION,
                        OOMPH_EXCEPTION_LOCATION);
       
      }
     else
      {
       //Read the index from the (cast) bulk element.
       U_index_poisson = eqn_pt->u_index_poisson();
      }
    }
    break;

    //Any other case is an error
    default:
     std::ostringstream error_stream; 
     error_stream <<  "Dimension of node is " << Dim 
                  << ". It should be 1,2, or 3!" << std::endl;
     
     throw OomphLibError(error_stream.str(),
                         OOMPH_CURRENT_FUNCTION,
                         OOMPH_EXCEPTION_LOCATION);
     break;
    }
  }


//===========================================================================
/// Compute the element's residual vector and the (zero) Jacobian matrix.
//===========================================================================
template<class ELEMENT>
void PoissonWithSingularityFluxElement<ELEMENT>::
fill_in_generic_residual_contribution_poisson_flux(
 Vector<double> &residuals, DenseMatrix<double> &jacobian, 
 const unsigned& flag)
{



 /* oomph_info  */
 /*  << "In PoissonWithSingularityFluxElement:: ...residual... Poisson_sing_el_pt = "  */
 /*  << Poisson_sing_el_pt << " ndof  = " << ndof() << std::endl; */

 if (flag==1) 
  {
   oomph_info << "Never get here -- include derivs w.r.t. C\n";
   abort();
  }

 //Find out how many nodes there are
 const unsigned n_node = nnode();
  
 //Set up memory for the shape and test functions
 Shape psif(n_node), testf(n_node);
 
 //Set the value of Nintpt
 const unsigned n_intpt = integral_pt()->nweight();
 
 //Set the Vector to hold local coordinates
 Vector<double> s(Dim-1);
 
 //Integers to hold the local equation and unknown numbers
 int local_eqn=0;

 // Locally cache the index at which the variable is stored
 const unsigned u_index_poisson = U_index_poisson;

 //Loop over the integration points
 //--------------------------------
 for(unsigned ipt=0;ipt<n_intpt;ipt++)
  {

   //Assign values of s
   for(unsigned i=0;i<(Dim-1);i++) {s[i] = integral_pt()->knot(ipt,i);}
   
   //Get the integral weight
   double w = integral_pt()->weight(ipt);
   
   //Find the shape and test functions and return the Jacobian
   //of the mapping
   double J = shape_and_test(s,psif,testf);
   
   //Premultiply the weights and the Jacobian
   double W = w*J;
   
   //Need to find position to feed into flux function, initialise to zero
   Vector<double> interpolated_x(Dim,0.0);
   
   //Calculate coords
   for(unsigned l=0;l<n_node;l++) 
    {
     //Loop over velocity components
     for(unsigned i=0;i<Dim;i++)
      {
       interpolated_x[i] += nodal_position(l,i)*psif[l];
      }
    }
   
   // Get gradient of singular fct (incl. amplitude)
   Vector<double> grad_sing(Dim,0.0);
   if (Poisson_sing_el_pt!=0)
    {
     grad_sing=Poisson_sing_el_pt->gradient_of_singular_fct(interpolated_x);
    }

   // Compute outer unit normal at the specified local coordinate
   Vector<double> unit_normal(Dim);
   outer_unit_normal(s,unit_normal);
   
   // Get flux associated with singular fct
   double flux_sing=0.0;
   for (unsigned i=0;i<Dim;i++)
    {
     flux_sing+=unit_normal[i]*grad_sing[i];
    }

   //Get the imposed flux
   double flux=0.0;
   get_flux(interpolated_x,flux);
   

   // Subtract off the flux from the singular fct
   double flux_for_u_fe=flux-flux_sing;


   //Now add to the appropriate equations
   
   //Loop over the test functions
   for(unsigned l=0;l<n_node;l++)
    {
     local_eqn = nodal_local_eqn(l,u_index_poisson);
     /*IF it's not a boundary condition*/
     if(local_eqn >= 0)
      {
       //Add the prescribed flux terms
       residuals[local_eqn] -= flux_for_u_fe*testf[l]*W;
         
       // Imposed traction doesn't depend upon the solution, 
       // --> the Jacobian is always zero, so no Jacobian
       // terms are required
      }
    }
  }
}









}

#endif
