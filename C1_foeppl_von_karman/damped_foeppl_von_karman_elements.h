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
#ifndef OOMPH_THERMO_FVK_ELEMENTS_HEADER
#define OOMPH_THERMO_FVK_ELEMENTS_HEADER

#include<sstream>

//OOMPH-LIB headers
#include "foeppl_von_karman_elements.h"


namespace oomph
{
//=============================================================
/// A class for all subparametric elements that solve the 2D-
/// Biharmonic equations.
/// \f[
/// \frac{\partial^4 u}{\partial x_i^4} = f(x_j)
/// \f]
/// This contains the generic maths. Shape functions, geometric
/// mapping etc. must get implemented in derived class.
//=============================================================
class DampedFoepplVonKarmanEquations
 : public virtual FoepplVonKarmanEquations
{

public:

 /// \short Function pointer to an alphaDT type swelling function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*SwellingFctPt)(const Vector<double>& x, double& f);

 /// \short Function pointer to a scalar prestrain function fct(x,f(x)) --
 /// x is a Vector!
 typedef void (*PrestrainFctPt)(const Vector<double>& x, double& f);

 /// Constructor (must initialise the Pressure_fct_pt to null)
 DampedFoepplVonKarmanEquations() :
  FoepplVonKarmanEquations(),
  Prestrain_fct_pt(0),
  Swelling_fct_pt(0)
 {}
 
 /// Broken copy constructor
 DampedFoepplVonKarmanEquations(const DampedFoepplVonKarmanEquations& dummy)
  {
   BrokenCopy::broken_copy("DampedFoepplVonKarmanEquations");
  }

 /// Broken assignment operator
 void operator=(const DampedFoepplVonKarmanEquations&)
  {
   BrokenCopy::broken_assign("DampedFoepplVonKarmanEquations");
  }

 /// Output with default number of plot points
 void output(std::ostream &outfile)
  {
   const unsigned n_plot=5;
   DampedFoepplVonKarmanEquations::output(outfile,n_plot);
  }

 /// \short Output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(std::ostream &outfile, const unsigned &n_plot);

 /// C_style output with default number of plot points
 void output(FILE* file_pt)
  {
   const unsigned n_plot=5;
   DampedFoepplVonKarmanEquations::output(file_pt,n_plot);
  }

 /// \short C-style output FE representation of soln: x,y,u or x,y,z,u at
 /// n_plot^DIM plot points
 void output(FILE* file_pt, const unsigned &n_plot);

 
 /// Output: x, y, sigma_xx, sigma_xy, sigma_yy,
 ///         (sigma_1, sigma_2, sigma_1x, sigma1y, sigma2x, sigma2y)
 ///                                   (if principal_stresses==true)
 /// at the Gauss integration points to obtain a smooth point
 /// cloud (stress may be discontinuous across elements in general)
 void output_smooth_stress(std::ostream &outfile,
			   const bool &principal_stresses=false);


 
 /// Output exact soln: x,y,u_exact or x,y,z,u_exact at n_plot^DIM plot points
 void output_fct(std::ostream &outfile, const unsigned &n_plot,
                 FiniteElement::SteadyExactSolutionFctPt exact_soln_pt);

 /// \short Output exact soln: x,y,u_exact or x,y,z,u_exact at
 /// n_plot^DIM plot points (dummy time-dependent version to
 /// keep intel compiler happy)
 virtual void output_fct(std::ostream &outfile, const unsigned &n_plot,
                         const double& time,
                         FiniteElement::UnsteadyExactSolutionFctPt
                         exact_soln_pt)
  {
   throw OomphLibError(
    "There is no time-dependent output_fct() for these elements ",
    OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }

 
 /// Fill in the strain tensor
 void get_epsilon(DenseMatrix<double>& epsilon,
		  const DenseMatrix<double>& grad_u,
		  const DenseMatrix<double>& grad_w,
		  const double& c_swell,
		  const double& prestrain)const
 {
  // Truncated Green Lagrange strain tensor
  DenseMatrix<double> dummy_epsilon(this->dim(),this->dim(),0.0);
  for(unsigned alpha=0;alpha<this->dim();++alpha)
   {
    for(unsigned beta=0;beta<this->dim();++beta)
     {
      // Truncated Green Lagrange strain tensor
      dummy_epsilon(alpha,beta) += 0.5* grad_u(alpha,beta)
       + 0.5*grad_u(beta,alpha)
       + 0.5*grad_w(0,alpha)*grad_w(0,beta);
     }
    // Swelling & prestrain
    dummy_epsilon(alpha,alpha) += prestrain - c_swell;
    
   }
  epsilon=dummy_epsilon;
 }

 /// Fill in the stress tensor
 void get_sigma(DenseMatrix<double>& sigma,
		const DenseMatrix<double>& grad_u,
		const DenseMatrix<double>& grad_w,
		const double c_swell,
		const double prestrain)const
 {
  // Get the Poisson ratio
  double nu(get_nu());

  // Truncated Green Lagrange strain tensor
  DenseMatrix<double> epsilon(this->dim(),this->dim(),0.0);
  get_epsilon(epsilon, grad_u, grad_w, c_swell, prestrain);
   
  // Now construct the Stress
  for(unsigned alpha=0;alpha<this->dim();++alpha)
   {
    for(unsigned beta=0;beta<this->dim();++beta)
     {
      // The Laplacian term: Trace[ \epsilon ] I
      // \nu * \epsilon_{\alpha \beta} delta_{\gamma \gamma}
      sigma(alpha,alpha) += nu*epsilon(beta,beta)/(1-nu*nu);
      
      // The scalar transform term: \epsilon
      // (1-\nu) * \epsilon_{\alpha \beta}
      sigma(alpha,beta) += (1-nu)* epsilon(alpha,beta)/(1-nu*nu);
     }
   }
 }
 
 // TODO: Fix this function to use index loops
 /// Fill in the stress tensor using a precalculated epsilon
 void get_sigma_from_epsilon(DenseMatrix<double>& sigma,
			     const DenseMatrix<double>& epsilon)const
 {
  // Get the Poisson ratio
  double nu(get_nu());
   
  // Now construct the Stress
  sigma(0,0) = (epsilon(0,0) + nu*epsilon(1,1)) / (1.0 - nu*nu);
  sigma(1,1) = (epsilon(1,1) + nu*epsilon(0,0)) / (1.0 - nu*nu);
  sigma(0,1) = epsilon(0,1) / (1.0 + nu);
  sigma(1,0) = sigma(0,1);
 }
 
 /// Get the principal stresses from the stress tensor
 void get_principal_stresses(const DenseMatrix<double>& sigma,
		Vector<double>& eigenvals,
		DenseMatrix<double>& eigenvecs)const
 {
  // Ensure that our eigenvectors are the right size
  eigenvals.resize(2);
  eigenvecs.resize(2);
  
  // Store the axial and shear stresses
  double s00 = sigma(0,0);
  double s01 = sigma(0,1);
  double s11 = sigma(1,1);


  // Handle the shear free case (sorted so max is first and min is second)
  if(s01==0.0)
   {
    // If x-stress is larger make it the first principal stress.
    if(s00>=s11)
     {
      eigenvals[0] = s00;
      eigenvals[1] = s11;
      eigenvecs(0,0)=1.0;
      eigenvecs(1,0)=0.0;
      eigenvecs(0,1)=0.0;
      eigenvecs(1,1)=1.0;
     }
    // else make y-stress the first principal stress.
    else
     {
      eigenvals[0] = s11;
      eigenvals[1] = s00;
      eigenvecs(0,0)=0.0;
      eigenvecs(1,0)=1.0;
      eigenvecs(0,1)=1.0;
      eigenvecs(1,1)=0.0;      
     }
   }
  
  // If there is shear stress calculate the principal stresses.
  else
   {
#ifdef PARANOID
    // Check that the discriminant of the characteristic polynomial of stress is
    // non-negative to within some tolerance.
    double neg_tol = -1.0e-10;
    double discriminant = (s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01);

    // If discriminant is less than tolerance then throw an error
    if (discriminant<neg_tol)
     {
      std::ostringstream error_stream;
      error_stream << "When calculating the principal stresses (eigenvalues of"
		   << " the stress tensor), the discriminant of the"
		   << " characteristic polynomial was " << discriminant
		   << " which is negative." << std::endl;
      throw OomphLibError(error_stream.str(),
			  OOMPH_CURRENT_FUNCTION,
			  OOMPH_EXCEPTION_LOCATION);
     }
    // // Else if discriminant is greater than tolerance but still negative
    // // then throw a warning
    // else if (discriminant<0.0)
    //  {
    //   std::ostringstream error_stream;
    //   error_stream << "Careful, when calculating the principal stresses "
    //                << "(eigenvalues of the stress tensor), the discriminant "
    // 		   << "of the characteristic polynomial was " << discriminant
    // 		   << " which is negative." << std::endl;
    //   OomphLibWarning(error_stream.str(),
    // 		      OOMPH_CURRENT_FUNCTION,
    // 		      OOMPH_EXCEPTION_LOCATION);
    //  }
#endif
    // We take the fabs() "floating point absolute" before sqrt() to catch
    // small negative FP errors.
    eigenvals[0] =
     0.5 * ( (s00 + s11) + sqrt(fabs((s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01))) );
    eigenvals[1] =
     0.5 * ( (s00 + s11) - sqrt(fabs((s00+s11)*(s00+s11) - 4.0*(s00*s11-s01*s01))) );

    // For max eval we choose y-positive evecs (suited to Benjamins sheet
    // problem)
    double sign = (eigenvals[0]-s00<0.0) ? -1.0 : 1.0;
    // Calculate the normalised principal stress direction for eigenvals[0]
    eigenvecs(0,0) =
     sign * (s01 / sqrt(s01*s01 + (eigenvals[0]-s00)*(eigenvals[0]-s00)));
    eigenvecs(1,0) =
     sign * ((eigenvals[0]-s00) / sqrt(s01*s01 + (eigenvals[0]-s00)*(eigenvals[0]-s00)));

    // For min eval we choose x-positive evecs (suited to Benjamins sheet
    // problem)
    sign = (s01<0.0) ? -1.0 : 1.0;
    // Calculate the normalised principal stress direction for eigenvals[1]
    eigenvecs(0,1) =
     sign * (s01 / sqrt(s01*s01 + (eigenvals[1]-s00)*(eigenvals[1]-s00)));
    eigenvecs(1,1) =
     sign * ((eigenvals[1]-s00) / sqrt(s01*s01 + (eigenvals[1]-s00)*(eigenvals[1]-s00)));
   }
 }
 
 /// Access function: Pointer to prestrain function
 PrestrainFctPt& prestrain_fct_pt() {return Prestrain_fct_pt;}
 
 /// Access function: Pointer to prestrain function. Const version
 PrestrainFctPt prestrain_fct_pt() const {return Prestrain_fct_pt;}
 
 /// Access function: Pointer to swelling function
 SwellingFctPt& swelling_fct_pt() {return Swelling_fct_pt;}
 
 /// Access function: Pointer to swelling function. Const version
 SwellingFctPt swelling_fct_pt() const {return Swelling_fct_pt;}
 
 /// Get swelling at (Eulerian) position x. This function is
 /// virtual to allow overloading.
 inline virtual void get_prestrain(const unsigned& ipt,
				   const Vector<double>& x,
				   double& epsilon_0) const
  {
   //If no prestrain function has been set, return zero
   if(Prestrain_fct_pt==0)
    {
     epsilon_0 = 0.0;
    }
   else
    {
     // Get prestrain magnitude
     (*Prestrain_fct_pt)(x,epsilon_0);
    }
  }
 
 /// Get swelling at (Eulerian) position x. This function is
 /// virtual to allow overloading.
 inline virtual void get_swelling_foeppl_von_karman(const unsigned& ipt,
						    const Vector<double>& x,
						    double& swelling) const
  {
   //If no swelling function has been set, return zero
   if(Swelling_fct_pt==0)
    {
     swelling = 0.0;
    }
   else
    {
     // Get swelling magnitude
     (*Swelling_fct_pt)(x,swelling);
    }
  }

 
 /// Number of 'flux' terms for Z2 error estimation
 unsigned num_Z2_flux_terms()
 {
  return 3;
 }
 

 /// Get 'flux' for Z2 error recovery:   Upper triangular entries
 /// in stress tensor.
 void get_Z2_flux(const Vector<double>& s, Vector<double>& flux)
 {
  unsigned dim = this->dim();
#ifdef PARANOID
  unsigned num_entries = num_Z2_flux_terms();
  if (flux.size() != num_entries)
   {
    std::ostringstream error_message;
    error_message << "The flux vector has the wrong number of entries, "
		  << flux.size() << ", whereas it should be " << num_entries
		  << std::endl;
    throw OomphLibError(error_message.str(),
			OOMPH_CURRENT_FUNCTION,
			OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Interpolate unknowns to get the displacement gradients
  Vector<double> u = interpolated_u_foeppl_von_karman(s);
  DenseMatrix<double> duidxj(dim,dim);
  DenseMatrix<double> dwdxi(1,dim);
  
  // Copy out gradient entries to the containers to pass to get_sigma
  // [zdec] dont hard code this
  duidxj(0,0) = u[8];
  duidxj(0,1) = u[9];
  duidxj(1,0) = u[10];
  duidxj(1,1) = u[11];
  dwdxi(0,0) = u[1];
  dwdxi(0,1) = u[2];
  
  // Get global x
  Vector<double> x(dim);
  interpolated_x(s,x);
  
  // Get degree of swelling and prestrain at x
  double C;
  double e0;
  get_swelling_foeppl_von_karman(0, x, C);
  get_prestrain(0, x, e0);
  
  // Get stress matrix
  DenseMatrix<double> epsilon(dim);
  DenseMatrix<double> sigma(dim);
  this->get_epsilon(epsilon, duidxj, dwdxi, C, e0);
  this->get_sigma_from_epsilon(sigma, epsilon);
  
  flux[0] = sigma(0,0);
  flux[1] = sigma(0,1);
  flux[2] = sigma(1,1);
 }
 

 /// Order of recovery shape functions for Z2 error estimation:
 /// Cubic.
 unsigned nrecovery_order()
 {
  return 3;
 }
 
 
 /// Add the element's contribution to its residual vector (wrapper)
 void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   //Call the generic residuals function with flag set to 0
   //using a dummy matrix argument
   fill_in_generic_residual_contribution_damped_foeppl_von_karman(
    residuals,GeneralisedElement::Dummy_matrix,0);
  }

 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian)
  {
   //Call the generic routine with the flag set to 1
   fill_in_generic_residual_contribution_damped_foeppl_von_karman(residuals,jacobian,1);
  }

 /// Add the element's contribution to its residual vector and
 /// element Jacobian matrix (wrapper)
 void fill_in_contribution_to_jacobian_and_mass_matrix(Vector<double> &residuals,
                                   DenseMatrix<double> &jacobian,DenseMatrix<double> &mass_matrix)
  {
   //Call fill in Jacobian 
   fill_in_contribution_to_jacobian(residuals,jacobian);
   // There is no mass matrix: we will just want J w = 0

   // -- COPIED FROM DISPLACMENT FVK EQUATIONS --
   // Dummy diagonal (won't result in global unit matrix but
   // doesn't matter for zero eigenvalue/eigenvector
   unsigned ndof=mass_matrix.nrow();
   for (unsigned i=0;i<ndof;i++)
    {
     mass_matrix(i,i)+=1.0;
    }
  }
 
 /// \short Calculate the elastic energy of the element and return it as a
 /// double.
 virtual Vector<double> element_elastic_and_kinetic_energy();



protected:
 
 /// \short Compute element residual Vector only (if flag=and/or element
 /// Jacobian matrix
 virtual void fill_in_generic_residual_contribution_damped_foeppl_von_karman(
  Vector<double> &residuals, DenseMatrix<double> &jacobian,
  const unsigned& flag);
 
 /// Pointer to prestrain function:
 PrestrainFctPt Prestrain_fct_pt; 
 
 /// Pointer to swelling function:
 SwellingFctPt Swelling_fct_pt; 
};


} //end namespace oomph
#endif
