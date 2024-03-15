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
// Non--inline functions for the foeppl_von_karman equations
#include "damped_foeppl_von_karman_elements.h"

namespace oomph
{
 //======================================================================
 /// Fill in generic residual and jacobian contribution 
 //======================================================================
 void  DampedFoepplVonKarmanEquations::
 fill_in_generic_residual_contribution_damped_foeppl_von_karman(
								Vector<double> &residuals,
								DenseMatrix<double> &jacobian,
								const unsigned& flag)
 {
  //Find the dimension of the element
  const unsigned dim = this->dim();
  //Find out how many nodes there are
  const unsigned n_u_node = nnode_inplane();
  const unsigned n_w_node = nnode_outofplane(); 
  //Find out how many bubble nodes there are
  const unsigned n_b_node = nbubble_basis();
  //Find out how many nodes positional dofs there are
  unsigned n_basis_type = nnodal_basis_type();
  // Find the internal dofs
  const unsigned n_b_position_type = nbubble_basis_type();

  //Get the Poisson ratio of the plate
  const double nu=get_nu();

  //Get the virtual dampening coefficient
  const double mu=get_mu();
  
  //Set up memory for the shape and test functions
  Shape psi_u(n_u_node), test_u(n_u_node);
  DShape dpsi_udxi(n_u_node,dim), dtest_udxi(n_u_node,dim);

  //Local c1-shape funtion
  Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
   psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
 
  DShape dpsi_dxi(n_w_node,n_basis_type,dim),dtest_dxi(n_w_node,n_basis_type,dim),
   dpsi_b_dxi(n_b_node,n_b_position_type,dim),dtest_b_dxi(n_b_node,n_b_position_type,dim),
   d2psi_dxi2(n_w_node,n_basis_type,3), d2test_dxi2(n_w_node,n_basis_type,3),
   d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

  //Set the value of n_intpt
  const unsigned n_intpt = this->integral_pt()->nweight();
  //Integers to store the local equation and unknown numbers
  int local_eqn=0, local_unknown=0;

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);

    //Calculate values of unknown
    Vector<double> interpolated_w(1,0.0);
    Vector<double> interpolated_dwdt(1,0.0);
    DenseMatrix<double> interpolated_dwdxi(1,dim,0.0);
    DenseMatrix<double> interpolated_d2wdxi2(1,3,0.0);
   
    //Calculate values of unknown
    Vector<double> interpolated_u(2,0.0);
    DenseMatrix<double> interpolated_duidxj(2,dim,0.0);

    //Allocate and initialise to zero
    Vector<double> interp_x(dim,0.0);
    Vector<double> s(dim);
    s[0] = this->integral_pt()->knot(ipt,0);
    s[1] = this->integral_pt()->knot(ipt,1);

    interpolated_x(s,interp_x);
    //Call the derivatives of the shape and test functions for the unknown
    double J = d2shape_and_d2test_eulerian_foeppl_von_karman(s,
							     psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
							     test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);

    dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,psi_u,dpsi_udxi,test_u,dtest_udxi);
    //Premultiply the weights and the Jacobian
    double W = w*J;

    //Calculate function value and derivatives
    //-----------------------------------------
    // Loop over nodes
    for(unsigned l=0;l<n_w_node;l++)
     {
      TimeStepper* timestepper_pt = this->node_pt(l)->time_stepper_pt();
      unsigned n_time = 0;
      if( !(timestepper_pt->is_steady()) )
       {
	n_time=timestepper_pt->ntstorage();
       }
      for(unsigned k=0;k<n_basis_type;k++)
       {
	//Get the nodal value of the unknown
	double w_value = this->raw_nodal_value(l,k+2);
	interpolated_w[0] += w_value*psi(l,k);
	// Loop over directions
	for(unsigned j=0;j<dim;j++)
	 {
	  interpolated_dwdxi(0,j) += w_value*dpsi_dxi(l,k,j);
	 }
	for(unsigned j=0;j<3;j++)
	 {
	  interpolated_d2wdxi2(0,j) += w_value*d2psi_dxi2(l,k,j);
	 }

	// Add the contributions to the time derivative from node l, type k
	double dwdt_value = 0.0;
	for(unsigned t=0; t<n_time; t++)
	 {
	  dwdt_value +=
	   timestepper_pt->weight(1,t) * this->raw_nodal_value(t, l, k+2);
	 }
	interpolated_dwdt[0] += dwdt_value * psi(l,k);
       }
     }

    // Loop over internal dofs
    // oomph_info << "About to loop over the bubble bases" << std::endl;
    for(unsigned l=0;l<nbubble_basis();l++)
     {
      // oomph_info << "Basis number: " << l << std::endl;
      // THIS NEEDS TO BE FIXED, NOT GENERAL ENOUGH.
      // (1) should be changed to right index
      TimeStepper* timestepper_pt = internal_data_pt(1) // :(
       ->time_stepper_pt();
      unsigned n_time = 0;
      if( !(timestepper_pt->is_steady()) )
       {
	n_time=timestepper_pt->ntstorage();
       }
      
      for(unsigned k=0;k<n_b_position_type;k++)
       {
	//Get the nodal value of the unknown
	double w_value = get_w_bubble_dof(l,k);
	interpolated_w[0] += w_value*psi_b(l,k);
	// Add the contributions to the time derivative
	double dwdt_value = 0.0;
	for(unsigned t=0; t<n_time; t++)
	 {
	  // oomph_info << "At time " << t << ": ";
	  dwdt_value +=
	   timestepper_pt->weight(1,t) * get_w_bubble_dof(l,k,t);
	  // oomph_info << dwdt_value
	  // 	     << " = " << timestepper_pt
	  // 	     << " * " << get_w_bubble_dof(l,k,t) << std::endl;;
	 }
	interpolated_dwdt[0] += dwdt_value * psi_b(l,k);
	// Loop over directions
	for(unsigned j=0;j<dim;j++)
	 {
	  interpolated_dwdxi(0,j) += w_value*dpsi_b_dxi(l,k,j);
	 }
	for(unsigned j=0;j<3;j++)
	 {
	  interpolated_d2wdxi2(0,j) += w_value*d2psi_b_dxi2(l,k,j);
	 }
       }
     }

    // Loop over nodes
    for(unsigned l=0;l<n_u_node;l++)
     {
      for(unsigned alpha=0;alpha<2;alpha++)
       {
	double u_value = this->raw_nodal_value(l,alpha);
	interpolated_u[alpha] += u_value*psi_u(l);
  
        // Loop over directions
        for(unsigned beta=0; beta<dim; beta++)
         {
          interpolated_duidxj(alpha,beta) += u_value*dpsi_udxi(l,beta);
         }
       }
     }

     
    //Get the swelling and prestrain
    //--------------------
    double c_swell(0.0);
    double epsilon_0(0.0);
    
    get_swelling_foeppl_von_karman(ipt,interp_x,c_swell);
    get_prestrain(ipt,interp_x,epsilon_0);
    
    // Get the stress
    DenseMatrix<double> sigma(2,2,0.0);
    get_sigma(sigma, interpolated_duidxj, interpolated_dwdxi,
	      c_swell, epsilon_0); 
   
    // Check 
    //    Vector<double> interpolated_w(6,0.0)pli:
    //    interpolated_w=this-> interpolated_w_foeppl_von_karman(s);
    //   interpolated_w[0]=interpolated_w[0];
    //   interpolated_dwdxi(0,0)=interpolated_w[1];
    //   interpolated_dwdxi(0,1)=interpolated_w[2];
    //   interpolated_d2wdxi2(0,0)=interpolated_w[3];
    //   interpolated_d2wdxi2(0,1)=interpolated_w[4];
    //   interpolated_d2wdxi2(0,2)=interpolated_w[5];
    //Get pressure function
    //-------------------
    double pressure(0.0);
    Vector<double> pressure_gradient(2,0.0);
    get_pressure_foeppl_von_karman(ipt,interp_x,pressure);
    get_in_plane_forcing_foeppl_von_karman(ipt,interp_x,pressure_gradient);

    // Loop over the nodal test functions
    for(unsigned l=0;l<n_w_node;l++)
     {
      for(unsigned k=0;k<n_basis_type;k++)
       {
	//Get the local equation
	local_eqn = this->nodal_local_eqn(l,k+2);
	//IF it's not a boundary condition
	if(local_eqn >= 0)
	 {
	  // Add virtual time derivative term for dampened buckling
	  residuals[local_eqn] += mu*interpolated_dwdt[0]*test(l,k)*W;
	  // Add body force/pressure term here
	  residuals[local_eqn] -= pressure*test(l,k)*W;
	  for(unsigned alpha=0;alpha<2;++alpha)
	   {
	    for(unsigned beta=0; beta<2;++beta)
	     {
	      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
	      residuals[local_eqn] += (1-nu)
	       *interpolated_d2wdxi2(0,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
	      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
	      residuals[local_eqn] += (nu)
	       *interpolated_d2wdxi2(0,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
	      // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
	      residuals[local_eqn] += eta()*sigma(alpha,beta)
	       *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;
	     }
	   }
	  // Calculate the jacobian
	  //-----------------------
	  if(flag)
	   {
	    //Loop over the in--plane unknowns again
	    for(unsigned ll=0;ll<n_u_node;ll++)
	     {
	      // Loop over displacements
	      for(unsigned gamma=0;gamma<2;++gamma)
	       {
		local_unknown = this->nodal_local_eqn(ll,gamma);
		// If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Loop over displacements
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    // Loop over dimensions 
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Nonlinear Terms 
		      if(alpha==beta)
		       {
			jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)* 
			 interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(alpha==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
			 *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(beta==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
			 *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		     }
		   }
		 }
	       }
	     }
	    
	    //Loop over the test functions again
	    for(unsigned l2=0;l2<n_w_node;l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_basis_type;k2++)
	       {
		local_unknown = this->nodal_local_eqn(l2,k2+2);
		// If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // Linear Terms
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
		       *dpsi_dxi(l2,k2,alpha)*dtest_dxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
		     }
		   }
		 }
	       }
	     }
	    //Loop over the internal test functions
	    for(unsigned l2=0;l2<nbubble_basis();l2++)
	     {
	      for(unsigned k2=0;k2<n_b_position_type;k2++)
	       {
		local_unknown = local_w_bubble_equation(l2,k2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi_b(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Linear Terms 
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_dxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_dxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_dxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_dxi(l,k,beta)*W;//DGR: eta^2 surely
		     }
		   }
		 }
	       }
	     }
	   } // End of flag
	 }
       }
     }

    // Loop over the internal test functions
    for(unsigned l=0;l<nbubble_basis();l++)
     {
      for(unsigned k=0;k<n_b_position_type;k++)
       {
	//Get the local equation
	local_eqn = local_w_bubble_equation(l,k);
	//IF it's not a boundary condition
	if(local_eqn >= 0)
	 {
	  // Add virtual time derivative term for dampened buckling
	  residuals[local_eqn] += mu*interpolated_dwdt[0]*test_b(l,k)*W;     
	  // Add body force/pressure term here
	  residuals[local_eqn] -= pressure*test_b(l,k)*W;
	  for(unsigned alpha=0;alpha<2;++alpha)
	   {
	    for(unsigned beta=0; beta<2;++beta)
	     {
	      // The Linear Terms 
	      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
	      residuals[local_eqn] += (1-nu)
	       *interpolated_d2wdxi2(0,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
	      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
	      residuals[local_eqn] += (nu)
	       *interpolated_d2wdxi2(0,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
	      // sigma_{\alpha\beta} w_\alpha \delta w_{,\beta}
	      residuals[local_eqn] += eta()*sigma(alpha,beta)
	       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
	     }
	   }
	  // Calculate the jacobian
	  //-----------------------
	  if(flag)
	   {
	    //Loop over the in--plane unknowns again
	    for(unsigned ll=0;ll<n_u_node;ll++)
	     {
	      // Loop over displacements
	      for(unsigned gamma=0;gamma<2;++gamma)
	       {
		local_unknown = this->nodal_local_eqn(ll,gamma);
		// If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Loop over displacements
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    // Loop over dimensions 
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Nonlinear Terms 
		      if(alpha==beta)
		       {
			jacobian(local_eqn,local_unknown) += eta()*nu*dpsi_udxi(ll,gamma)*
			 interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(alpha==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,beta)
			 *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		      if(beta==gamma)
		       {
			jacobian(local_eqn,local_unknown) += eta()*(1-nu)/2.*dpsi_udxi(ll,alpha)
			 *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)/(1-nu*nu)*W;
		       }
		     }
		   }
		 }
	       }

	     }
	    //Loop over the test functions again
	    for(unsigned l2=0;l2<n_w_node;l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_basis_type;k2++)
	       {
		local_unknown = this->nodal_local_eqn(l2,k2+2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test_b(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Linear Terms 
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
		       *dpsi_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_dxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		     }
		   }
		 }
	       }
	     }
	    //Loop over the test functions again
	    for(unsigned l2=0;l2<nbubble_basis();l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_b_position_type;k2++)
	       {
		local_unknown = local_w_bubble_equation(l2,k2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Contribution from buckle stabilising drag
		  jacobian(local_eqn,local_unknown) +=
		   psi_b(l2,k2)*node_pt(l2)->time_stepper_pt()->weight(1,0)
		   *mu*test_b(l,k)*W;
		  // Loop over dimensions 
		  for(unsigned alpha=0;alpha<2;++alpha)
		   {
		    for(unsigned beta=0; beta<2;++beta)
		     {
		      // The Linear Terms 
		      // w_{,\alpha\beta} \delta \kappa_\alpha\beta
		      jacobian(local_eqn,local_unknown) += (1-nu)
		       *d2psi_b_dxi2(l2,k2,alpha+beta)*d2test_b_dxi2(l,k,alpha+beta)*W;
		      // w_{,\alpha\alpha} \delta \kappa_\beta\beta
		      jacobian(local_eqn,local_unknown) += nu
		       *d2psi_b_dxi2(l2,k2,beta+beta)*d2test_b_dxi2(l,k,alpha+alpha)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*sigma(alpha,beta) 
                       *dpsi_b_dxi(l2,k2,alpha)*dtest_b_dxi(l,k,beta)*W;
		      // Nonlinear terms
		      jacobian(local_eqn,local_unknown) += eta()*nu/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,beta)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,alpha)*dpsi_b_dxi(l2,k2,beta)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		      jacobian(local_eqn,local_unknown) += eta()/2.*(1-nu)/(1-nu*nu)*
		       interpolated_dwdxi(0,beta)*dpsi_b_dxi(l2,k2,alpha)
		       *interpolated_dwdxi(0,alpha)*dtest_b_dxi(l,k,beta)*W;
		     }
		   }
		 }
	       }
	     }
	   } // End of flag
	 }
       }
     } //end loop over nodes
    for(unsigned l=0;l<n_u_node;++l)
     {
      // Now loop over displacement equations
      for(unsigned alpha=0; alpha<2 ; ++alpha)
       {
	//Get the local equation
	local_eqn = this->nodal_local_eqn(l,alpha);
	//IF it's not a boundary condition
	if(local_eqn >= 0)
	 {
	  // Now add on the stress terms
	  for(unsigned beta=0; beta<2;++beta)
	   {
	    // Variations in u_alpha
	    residuals[local_eqn] += sigma(alpha,beta)*dtest_udxi(l,beta)*W;
	   }
	  // Add the forcing term
	  residuals[local_eqn] += pressure_gradient[alpha]*psi_u(l)*W;

	  // Now loop over Jacobian
	  if(flag)
	   {
	    for(unsigned ll=0;ll<n_u_node;++ll)
	     {
	      // Now loop over displacement equations
	      for(unsigned gamma=0; gamma<2 ; ++gamma)
	       {
		//Get the local equation
		local_unknown = this->nodal_local_eqn(ll,gamma);
		// If the unknown is a degree of freedom 
		if(local_unknown>=0)
		 {
		  // Now add on the stress terms
		  for(unsigned beta=0; beta<2;++beta)
		   {
		    // The Linear Terms 
		    if(alpha==beta)
		     {
		      jacobian(local_eqn,local_unknown) += nu*dpsi_udxi(ll,gamma)
		       *dtest_udxi(l,beta)/(1-nu*nu)*W;
		     }
		    if(alpha==gamma)
		     {
		      jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,beta)
		       *dtest_udxi(l,beta)/(1-nu*nu)*W;
		     }
		    if(beta==gamma)
		     {
		      jacobian(local_eqn,local_unknown) += (1-nu)/2.*dpsi_udxi(ll,alpha)
		       *dtest_udxi(l,beta)/(1-nu*nu)*W;
		     }
		   }// End loop beta
		 }// End if local unknown
	       }// End loop gamma dof
	     }// End loop nodal dofs
       
	    //Loop over the w nodal dofs
	    for(unsigned l2=0;l2<n_w_node;l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_basis_type;k2++)
	       {
		//If at a non-zero degree of freedom add in the entry
		//Get the local equation
		local_unknown = this->nodal_local_eqn(l2,k2+2);
		// If the unknown is a degree of freedom 
		if(local_unknown>=0)
		 {
		  // Now add on the stress terms
		  for(unsigned beta=0; beta<2;++beta)
		   {
		    // The NonLinear Terms 
		    jacobian(local_eqn,local_unknown) += nu*dpsi_dxi(l2,k2,beta)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_dxi(l2,k2,alpha)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_dxi(l2,k2,beta)
		     *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		   }// End loop beta
		 }// End if local unknown
	       }// End loop position type dof
	     }// End loop w nodal dofs
       
	    //Loop over the w internal dofs
	    for(unsigned l2=0;l2<nbubble_basis();l2++)
	     {
	      // Loop over position dofs
	      for(unsigned k2=0;k2<n_b_position_type;k2++)
	       {
		//If at a non-zero degree of freedom add in the entry
		//Get the local equation
		local_unknown = local_w_bubble_equation(l2,k2);
		//If at a non-zero degree of freedom add in the entry
		if(local_unknown >= 0)
		 {
		  // Now add on the stress terms
		  for(unsigned beta=0; beta<2;++beta)
		   {
		    // The NonLinear Terms 
		    jacobian(local_eqn,local_unknown) += nu*dpsi_b_dxi(l2,k2,beta)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,alpha)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_b_dxi(l2,k2,alpha)
		     *interpolated_dwdxi(0,beta)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		    jacobian(local_eqn,local_unknown) +=(1-nu)/2.*dpsi_b_dxi(l2,k2,beta)
		     *interpolated_dwdxi(0,alpha)*dtest_udxi(l,beta)/(1-nu*nu)*W;
		   } // End loop beta
		 } // End if local_unknown
	       } // End loop position type
	     } // End loop internal dofs
	   } // End if flag
	 } // End if local eqn
       } // End loop over alpha
     } // End loop over nodes
   } // End of loop over integration points
 }

 // Define the energy calculation function here
 Vector<double> DampedFoepplVonKarmanEquations::element_elastic_and_kinetic_energy()
 {
 
  double element_elastic_energy = 0.0;
  double element_kinetic_energy = 0.0;
  const unsigned dim = this->dim();
  //Find out how many nodes there are
  const unsigned n_u_node = nnode_inplane();
  const unsigned n_w_node = nnode_outofplane(); 
  //Find out how many bubble nodes there are
  const unsigned n_b_node = nbubble_basis();
  //Find out how many nodes positional dofs there are
  unsigned n_basis_type = nnodal_basis_type();
  // Find the internal dofs
  const unsigned n_b_position_type = nbubble_basis_type();

  //Get the Poisson ratio of the plate
  const double nu=get_nu();

  //Get the virtual dampening coefficient
  const double mu=get_mu();
  
  //Set up memory for the shape and test functions
  Shape psi_u(n_u_node), test_u(n_u_node);
  DShape dpsi_udxi(n_u_node,dim), dtest_udxi(n_u_node,dim);

  //Local c1-shape funtion
  Shape psi(n_w_node,n_basis_type),test(n_w_node,n_basis_type),
   psi_b(n_b_node,n_b_position_type),test_b(n_b_node,n_b_position_type);
 
  DShape dpsi_dxi(n_w_node,n_basis_type,dim),dtest_dxi(n_w_node,n_basis_type,dim),
   dpsi_b_dxi(n_b_node,n_b_position_type,dim),dtest_b_dxi(n_b_node,n_b_position_type,dim),
   d2psi_dxi2(n_w_node,n_basis_type,3), d2test_dxi2(n_w_node,n_basis_type,3),
   d2psi_b_dxi2(n_b_node,n_b_position_type,3), d2test_b_dxi2(n_b_node,n_b_position_type,3);

  //Set the value of n_intpt
  const unsigned n_intpt = this->integral_pt()->nweight();

  //Loop over the integration points
  for(unsigned ipt=0;ipt<n_intpt;ipt++)
   {
    //Get the integral weight
    double w = this->integral_pt()->weight(ipt);

    //Calculate values of unknown
    Vector<double> interpolated_w(1,0.0);
    Vector<double> interpolated_dwdt(1,0.0);
    DenseMatrix<double> interpolated_dwdxi(1,dim,0.0);
    DenseMatrix<double> interpolated_d2wdxi2(1,3,0.0);
   
    //Calculate values of unknown
    Vector<double> interpolated_u(2,0.0);
    DenseMatrix<double> interpolated_duidxj(2,dim,0.0);

    //Allocate and initialise to zero
    Vector<double> interp_x(dim,0.0);
    Vector<double> s(dim);
    s[0] = this->integral_pt()->knot(ipt,0);
    s[1] = this->integral_pt()->knot(ipt,1);
    interpolated_x(s,interp_x);
    //Call the derivatives of the shape and test functions for the unknown
    double J = d2shape_and_d2test_eulerian_foeppl_von_karman(s,
							     psi, psi_b, dpsi_dxi, dpsi_b_dxi, d2psi_dxi2, d2psi_b_dxi2,
							     test, test_b, dtest_dxi, dtest_b_dxi, d2test_dxi2, d2test_b_dxi2);

    dshape_u_and_dtest_u_eulerian_foeppl_von_karman(s,psi_u,dpsi_udxi,test_u,dtest_udxi);
    //Premultiply the weights and the Jacobian
    double W = w*J;


    //Calculate function value and derivatives
    //-----------------------------------------
    // Loop over nodes
    for(unsigned l=0;l<n_w_node;l++)
     {
      TimeStepper* timestepper_pt = this->node_pt(l)->time_stepper_pt();
      unsigned n_time = 0;
      if( !(timestepper_pt->is_steady()) )
       {
	n_time=timestepper_pt->ntstorage();
       }
      for(unsigned k=0;k<n_basis_type;k++)
       {
	//Get the nodal value of the unknown
	double w_value = this->raw_nodal_value(l,k+2);
	interpolated_w[0] += w_value*psi(l,k);
	// Loop over directions
	for(unsigned j=0;j<dim;j++)
	 {
	  interpolated_dwdxi(0,j) += w_value*dpsi_dxi(l,k,j);
	 }
	for(unsigned j=0;j<3;j++)
	 {
	  interpolated_d2wdxi2(0,j) += w_value*d2psi_dxi2(l,k,j);
	 }

	// Add the contributions to the time derivative from node l, type k
	double dwdt_value = 0.0;
	for(unsigned t=0; t<n_time; t++)
	 {
	  dwdt_value +=
	   timestepper_pt->weight(1,t) * this->raw_nodal_value(t, l, k+2);
	 }
	interpolated_dwdt[0] += dwdt_value * psi(l,k);
       }
     }

    // Loop over internal dofs
    // oomph_info << "About to loop over the bubble bases" << std::endl;
    for(unsigned l=0;l<nbubble_basis();l++)
     {
      // oomph_info << "Basis number: " << l << std::endl;
      // THIS NEEDS TO BE FIXED, NOT GENERAL ENOUGH.
      // (1) should be changed to right index
      TimeStepper* timestepper_pt = internal_data_pt(1) // :(
       ->time_stepper_pt();
      unsigned n_time = 0;
      if( !(timestepper_pt->is_steady()) )
       {
	n_time=timestepper_pt->ntstorage();
       }
      
      for(unsigned k=0;k<n_b_position_type;k++)
       {
	//Get the nodal value of the unknown
	double w_value = get_w_bubble_dof(l,k);
	interpolated_w[0] += w_value*psi_b(l,k);
	// Add the contributions to the time derivative
	double dwdt_value = 0.0;
	for(unsigned t=0; t<n_time; t++)
	 {
	  // oomph_info << "At time " << t << ": ";
	  dwdt_value +=
	   timestepper_pt->weight(1,t) * get_w_bubble_dof(l,k,t);
	  // oomph_info << dwdt_value
	  // 	     << " = " << timestepper_pt
	  // 	     << " * " << get_w_bubble_dof(l,k,t) << std::endl;;
	 }
	interpolated_dwdt[0] += dwdt_value * psi_b(l,k);
	// Loop over directions
	for(unsigned j=0;j<dim;j++)
	 {
	  interpolated_dwdxi(0,j) += w_value*dpsi_b_dxi(l,k,j);
	 }
	for(unsigned j=0;j<3;j++)
	 {
	  interpolated_d2wdxi2(0,j) += w_value*d2psi_b_dxi2(l,k,j);
	 }
       }
     }

    // Loop over nodes
    for(unsigned l=0;l<n_u_node;l++)
     {
      for(unsigned alpha=0;alpha<2;alpha++)
       {
	double u_value = this->raw_nodal_value(l,alpha);
	interpolated_u[alpha] += u_value*psi_u(l);
  
	// Loop over directions
	for(unsigned beta=0; beta<dim; beta++)
	 {
	  interpolated_duidxj(alpha,beta) += u_value*dpsi_udxi(l,beta);
	 }
       }
     }

    //----------------------------------------------------------------------
    // Add the elastic energy at ipt.


    //Get the swelling and prestrain
    //--------------------
    double c_swell(0.0);
    double epsilon_0(0.0);
    get_swelling_foeppl_von_karman(ipt,interp_x,c_swell);
    get_prestrain(ipt,interp_x,epsilon_0);
        
    // Truncated Green Lagrange strain tensor and 2Piola stress tensor
    DenseMatrix<double> epsilon(dim,dim,0.0);
    get_epsilon(epsilon, interpolated_duidxj, interpolated_dwdxi,
		c_swell, epsilon_0);
    DenseMatrix<double> sigma(dim,dim,0.0);
    get_sigma_from_epsilon(sigma, epsilon);

    // Add the plane strain energy contribution at ipt.
    for(unsigned alpha=0;alpha<dim;++alpha)
     {
      for(unsigned beta=0;beta<dim;++beta)
       {
	element_elastic_energy +=
	 0.5*eta()*sigma(alpha,beta)*epsilon(alpha,beta)*W;
       }
     }
    
    // Add the bending energy contribution at ipt.
    element_elastic_energy +=
     0.5 * (
	    (interpolated_d2wdxi2(0,0)+interpolated_d2wdxi2(0,2))
	    *(interpolated_d2wdxi2(0,0)+interpolated_d2wdxi2(0,2))
	    + 2.0*(1.0-nu) * (
			      interpolated_d2wdxi2(0,1)*interpolated_d2wdxi2(0,1)
			      - interpolated_d2wdxi2(0,0)*interpolated_d2wdxi2(0,2)
			      )
	    ) * W;

    //----------------------------------------------------------------------
    // Add the kinetic energy at ipt.
    element_kinetic_energy += 0.5*mu*interpolated_dwdt[0]*interpolated_dwdt[0] * W;

   }
  Vector<double> energies(2,0.0);
  // First entry is the elastic energy and second is kinetic.
  energies[0] = element_elastic_energy;
  energies[1] = element_kinetic_energy;
  
  return energies;
 }


//======================================================================
/// Output function:
///
///   x,y,u,(epsilon),(sigma),(sigma_evals,sigma_evecs)
///
/// nplot points in each coordinate direction
//======================================================================
void DampedFoepplVonKarmanEquations::output(std::ostream &outfile,
					    const unsigned &nplot)
{
 unsigned dim = this->dim();
 
 //Vector of local coordinates
 Vector<double> s(dim),x(dim);

 // Tecplot header info
 outfile << this->tecplot_zone_string(nplot);

 // Storage for variables
 double c_swell(0.0);
 double epsilon_0(0.0);
 Vector<double> u;
 DenseMatrix<double> interpolated_dwdxj(1,dim,0.0);
 DenseMatrix<double> interpolated_duidxj(2,dim,0.0);
 DenseMatrix<double> epsilon(2,2,0.0);
 DenseMatrix<double> sigma(2,2,0.0);
 DenseMatrix<double> sigma_eigenvecs(2,2,0.0);
 Vector<double> sigma_eigenvals(2,0.0);
 unsigned num_plot_points=this->nplot_points(nplot);
 Vector<double> r(3);
 
 // Loop over plot points
 for (unsigned iplot=0;iplot<num_plot_points;iplot++)
  {
   // Get local and global coordinates of plot point
   this->get_s_plot(iplot,nplot,s);
   interpolated_x(s,x);

   // Get interpolated unknowns
   u = interpolated_u_foeppl_von_karman(s);

   // Get prestrain and degree of swelling for use in the strain tensor
   this->get_swelling_foeppl_von_karman(iplot,x,c_swell);
   this->get_prestrain(iplot,x,epsilon_0);
   
   // TODO: make the output customisable with flags e.g. output_axial_strain,
   //       output_principal_stress
   // TODO: make the indexing below more general and not hard coded.

   // Copy gradients from u into interpolated gradient matrices...
   interpolated_dwdxj(0,0) = u[1]; //dwdx1
   interpolated_dwdxj(0,1) = u[2]; //dwdx2
   interpolated_duidxj(0,0)= u[8]; //du1dx1
   interpolated_duidxj(0,1)= u[9]; //du1dx2
   interpolated_duidxj(1,0)= u[10]; //du2dx1
   interpolated_duidxj(1,1)= u[11]; //du2dx2
   // ...which are used to retrieve the strain tensor epsilon.
   get_epsilon(epsilon, interpolated_duidxj,interpolated_dwdxj,
	       c_swell, epsilon_0);
   
   // Use epsilon to find the stress sigma.
   get_sigma_from_epsilon(sigma, epsilon);
   
   // Get the principal values and the corresponding directions of stress.
   get_principal_stresses(sigma, sigma_eigenvals, sigma_eigenvecs);
   
   
   for(unsigned i=0;i<this->dim();i++)
    {
     outfile << x[i] << " ";
    }
   
   // Loop for variables
   for(Vector<double>::iterator it=u.begin();it!=u.end();++it)
    {
     outfile << *it << " " ;
    }
   
   // Output axial strains
   outfile << epsilon(0,0) << " " << epsilon(0,1) << " " << epsilon(1,1) << " ";

   // Output axial stress
   outfile << sigma(0,0) << " " << sigma(0,1) << " " << sigma(1,1) << " ";
   
   // Output principal stresses
   outfile << sigma_eigenvals[0] << " " << sigma_eigenvals[1] << " ";
   
   // Output principal stress directions
   outfile << sigma_eigenvecs(0,0) << " " << sigma_eigenvecs(1,0) << " "
	   << sigma_eigenvecs(0,1) << " " << sigma_eigenvecs(1,1) << " ";
   
   // End output line
   outfile << std::endl;
  }
 // Write tecplot footer (e.g. FE connectivity lists)
 this->write_tecplot_zone_footer(outfile,nplot);
}



//======================================================================
/// C-style output function: // [TODO]
///
///   x,y,u   or    x,y,z,u
///
/// nplot points in each coordinate direction
//======================================================================
void DampedFoepplVonKarmanEquations::output(FILE* file_pt,
                                    const unsigned &nplot)
{
  throw OomphLibError("C-style output not implemented for "
                      "DampedFoepplVonKarmanEquations yet.\n"
                      "Please use c++ output.",
		      OOMPH_CURRENT_FUNCTION,
		      OOMPH_EXCEPTION_LOCATION);
}


//======================================================================
 /// Output exact solution // [TODO]
 ///
 /// Solution is provided via function pointer.
 /// Plot at a given number of plot points.
 ///
 ///   x,y,u_exact    or    x,y,z,u_exact
//======================================================================
void DampedFoepplVonKarmanEquations::output_fct(std::ostream &outfile,
                                       const unsigned &nplot,
                  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
{
  throw OomphLibError("Exact solution output not implemented for "
                      "DampedFoepplVonKarmanEquations yet.",
                      OOMPH_CURRENT_FUNCTION,
                      OOMPH_EXCEPTION_LOCATION);
}


 



//======================================================================
 /// Output: x, y, sigma_xx, sigma_xy, sigma_yy,
 ///         (sigma_1, sigma_2, sigma_1x, sigma1y, sigma2x, sigma2y)
 ///                                   (if principal_stresses==true)
 /// at the Gauss integration points to obtain a smooth point
 /// cloud (stress may be discontinuous across elements in general)
//======================================================================
void DampedFoepplVonKarmanEquations::output_smooth_stress(std::ostream &outfile,
							  const bool &principal_stresses)
 {
  unsigned dim = this->dim();

  // Vector of local coordinates
  Vector<double> s(dim);
  // Vector of global coordinates
  Vector<double> x(dim);

  // Storage for variables
  double c_swell(0.0);
  double epsilon_0(0.0);
  Vector<double> u;
  DenseMatrix<double> interpolated_dwdxj(1,dim,0.0);
  DenseMatrix<double> interpolated_duidxj(2,dim,0.0);
  DenseMatrix<double> sigma(2,2,0.0);
  Vector<double> sigma_eigenvals(2,0.0);
  DenseMatrix<double> sigma_eigenvecs(2,2,0.0);

  // The number of plot points is the number of integration (Gauss) points
  unsigned num_plot_points=this->integral_pt()->nweight();

  double iplot=0;
  
  // Loop over the plot points
  for(unsigned iplot=0; iplot<num_plot_points; iplot++)
   {
    // Get the local and global coordinates of the plot point
    s[0] = this->integral_pt()->knot(iplot,0);
    s[1] = this->integral_pt()->knot(iplot,1);
    //s[0] = 0.5;
    //s[1] = 0.5;
    interpolated_x(s,x);
    
      // Get interpolated unknowns
   u = interpolated_u_foeppl_von_karman(s);

   // Get prestrain and degree of swelling for use in the strain tensor
   this->get_swelling_foeppl_von_karman(iplot,x,c_swell);
   this->get_prestrain(iplot,x,epsilon_0);
   
   // Copy gradients from u into interpolated gradient matrices...
   interpolated_dwdxj(0,0) = u[1]; //dwdx1
   interpolated_dwdxj(0,1) = u[2]; //dwdx2
   interpolated_duidxj(0,0)= u[8]; //du1dx1
   interpolated_duidxj(0,1)= u[9]; //du1dx2
   interpolated_duidxj(1,0)= u[10]; //du2dx1
   interpolated_duidxj(1,1)= u[11]; //du2dx2

   DenseMatrix<double> epsilon(2,2,0.0);
   get_epsilon(epsilon, interpolated_duidxj,interpolated_dwdxj,
	       c_swell, epsilon_0);

   // Use epsilon to find the stress sigma.
   get_sigma_from_epsilon(sigma, epsilon);

   // Output the global coordinates of the plot point
   for(unsigned i=0;i<this->dim();i++)
    {
     outfile << x[i] << " ";
    }
   
   // Output axial stresses [zdec] here also w for debugging
   outfile << sigma(0,0) << " " << sigma(0,1) << " " << sigma(1,1) << " ";
   
   // If we want the principal stresses, calculate them and add them to the
   // output
   if(principal_stresses)
    {
     // Get the principal values and the corresponding directions of stress.
     get_principal_stresses(sigma, sigma_eigenvals, sigma_eigenvecs);
     
     // Output principal stress magnitudes
     outfile << sigma_eigenvals[0] << " " << sigma_eigenvals[1] << " ";
     
     // Output principal stress directions
     outfile << sigma_eigenvecs(0,0) << " " << sigma_eigenvecs(1,0) << " "
	     << sigma_eigenvecs(0,1) << " " << sigma_eigenvecs(1,1) << " ";
    }
   
   // End output line
   outfile << std::endl;
  }
 }

} 


