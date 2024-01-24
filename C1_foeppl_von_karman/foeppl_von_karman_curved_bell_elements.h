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

#ifndef OOMPH_BIHARMONIC_CURVABLE_BELL_ELEMENTS_HEADER
#define OOMPH_BIHARMONIC_CURVABLE_BELL_ELEMENTS_HEADER

#include "foeppl_von_karman_equations.h"
#include "../C1_basis/Bell_element_basis.h"
#include "../C1_basis/C1_curved_elements.h"
#include "../C1_basis/my_geom_object.h"
#include "../C1_basis/SubparametricTElement.h"

namespace oomph
{
 //=============================================================================
 /// FoepplVonKarmanC1CurvableBellElement elements are a subparametric scheme
 /// with  linear Lagrange interpolation for approximating the geometry and
 /// the C1-functions for approximating variables.
 /// These elements use TElement<2,NNODE_1D> with Lagrangian bases to interpolate
 /// the in-plane unknowns u_x, u_y although with sub-parametric simplex shape
 /// interpolation.
 ///
 /// The Bell Hermite basis is used to interpolate the out-of-plane unknown w for
 /// straight sided elements and is upgraded to Bernadou basis if the element is
 /// on a curved boundary. This upgrade similarly corrects the shape
 /// interpolation.
 ///
 /// As a result, all nodes contain the two in-plane Lagrangian dofs, vertex
 /// nodes (i.e. nodes 0,1,2) each contain an additional six out-of-plane Hermite
 /// dofs.
 ///
 ///  e.g FeopplVonKarmanC1CurvableBellElement<4> (unupgraded)
 ///                        | 0 | 1 | 2 | 3  | 4  |  5   |  6    |  7   |
 ///          o <---------- |u_x|u_y| w |dwdx|dwdy|d2wdx2|d2wdxdy|d2wdy2|
 ///         / \            .   
 ///        x   x <-------- |u_x|u_y|
 ///       /     \          .
 ///      x   x   x         .
 ///     /         \        .
 ///    o---x---x---o       .
 ///                  
 ///  e.g FeopplVonKarmanC1CurvableBellElement<3> (upgraded)
 ///                         | 0 | 1 | 2 | 3  | 4  |  5   |  6    |  7   |
 ///          o_ <---------- |u_x|u_y| w |dwdn|dwdt|d2wdn2|d2wdndt|d2wdt2|
 ///         /  \            .
 ///        /    |           .
 ///       x     x <-------- |u_x|u_y|
 ///      /      |           .
 ///     /        \          .
 ///    o-----x----o         .
 ///                    
 //=============================================================================
 
 template <unsigned NNODE_1D>
 class FoepplVonKarmanC1CurvableBellElement :
  public virtual CurvableBellElement<NNODE_1D>,
  public virtual FoepplVonKarmanEquations
 {
 public:
  
  /// Constructor: Call constructors for C1CurvableBellElement and
  /// FoepplVonKarman equations
  FoepplVonKarmanC1CurvableBellElement() :
   CurvableBellElement<NNODE_1D>(Nfield,Field_is_bell_interpolated),
   FoepplVonKarmanEquations(),
   Rotated_basis_fct_pt(0),
   Nnodes_to_rotate(0)
  {
    // Use the higher order integration scheme
   TGauss<2,4>* new_integral_pt = new TGauss<2,4>;
   this->set_integration_scheme(new_integral_pt); 

   Association_matrix_pt=0;
  }

  /// Destructor: clean up alloacations
  ~FoepplVonKarmanC1CurvableBellElement()
  {
  }

  /// Broken copy constructor
  FoepplVonKarmanC1CurvableBellElement(const
				     FoepplVonKarmanC1CurvableBellElement<NNODE_1D>& dummy)
  {
   BrokenCopy::broken_copy("FoepplVonKarmanC1CurvableBellElement");
  }

  /// Broken assignment operator
  void operator=(const FoepplVonKarmanC1CurvableBellElement<NNODE_1D>&)
  {
   BrokenCopy::broken_assign("FoepplVonKarmanC1CurvableBellElement");
  }

  /// Function pointer to basis vectors function which sets  basis vectors
  /// b1 and b2 (which are in general functions of x)
  typedef void (*BasisVectorsFctPt) (const Vector<double>& x,
				     Vector<double>& b1,
				     Vector<double>& b2,
				     DenseMatrix<double>& Db1,
				     DenseMatrix<double>& Db2);

  /// Get the zeta coordinate
  inline void interpolated_zeta(const Vector<double> &s,
				Vector<double> &zeta) const
  {
   // If there is a macro element use it
   if(this->Macro_elem_pt!=0) {this->get_x_from_macro_element(s,zeta);}
   //Otherwise interpolate zeta_nodal using the shape functions
   else {interpolated_x(s,zeta);}
  }

  /// Lagrange interpolated shape for in--plane displacements
  void shape_u(const Vector<double> &s, Shape &psi) const;
  
  /// Lagrange interpolated d_shape for in--plane displacements
  void dshape_u_local(const Vector<double> &s,
		      Shape &psi, DShape &dpsids) const;

  /// Add the element's contribution to its residual vector (wrapper) with cached
  /// association matrix
  void fill_in_contribution_to_residuals(Vector<double> &residuals)
  {
   // Store the expensive-to-construct matrix
   this->store_association_matrix();
   //Call the generic routine with the flag set to 1
   FoepplVonKarmanEquations::fill_in_contribution_to_residuals(residuals);
   // Remove the expensive-to-construct matrix
   this->delete_association_matrix();
  }
  
  /// Add the element's contribution to its residual vector and
  /// element Jacobian matrix (wrapper) with caching of association matrix
  void fill_in_contribution_to_jacobian(Vector<double> &residuals,
					DenseMatrix<double> &jacobian)
  {
   // Store the expensive-to-construct matrix
   this->store_association_matrix();
   //Call the generic routine with the flag set to 1
   FoepplVonKarmanEquations::fill_in_contribution_to_jacobian(residuals,jacobian);
   // Remove the expensive-to-construct matrix
   this->delete_association_matrix();
  }

  //----------------------------------------------------------------------------
  // Interface to FoepplVonKarmanEquations (can this all be (static) data?)

  /// Field indices for u
  virtual Vector<unsigned> u_field_indices() const
  { return {0,1}; }
  
  /// Field indices for u
  virtual unsigned u_alpha_field_index(const unsigned& alpha) const
  { return alpha; }
  
  /// Field index for w
  virtual unsigned w_field_index() const
  { return 2; }
  
  /// Interface to return the number of nodes used by u
  virtual unsigned nu_node() const
  { return CurvableBellElement<NNODE_1D>::nnode_for_field(u_alpha_field_index(0)); }

  /// Interface to return the number of nodes used by w
  virtual unsigned nw_node() const
  { return CurvableBellElement<NNODE_1D>::nnode_for_field(w_field_index()); }

  
  /// Interface to get the local indices of the nodes used by u
  virtual Vector<unsigned> get_u_node_indices() const
  { return CurvableBellElement<NNODE_1D>::nodal_indices_for_field(u_alpha_field_index(0)); }
  
  /// Interface to get the local indices of the nodes used by w
  virtual Vector<unsigned> get_w_node_indices() const
  { return CurvableBellElement<NNODE_1D>::nodal_indices_for_field(w_field_index()); }

  
  /// Interface to get the number of basis types for u at node j
  virtual unsigned nu_type_at_each_node() const
  { return CurvableBellElement<NNODE_1D>::nnodal_basis_type_for_field(u_alpha_field_index(0)); }

  /// Interface to get the number of basis types for w at node j
  virtual unsigned nw_type_at_each_node() const
  { return CurvableBellElement<NNODE_1D>::nnodal_basis_type_for_field(w_field_index()); }

  
  /// Interface to retrieve the value of u_alpha at node j of
  /// type k
  virtual double get_u_alpha_value_at_node_of_type(const unsigned& alpha,
						   const unsigned& j_node,
						   const unsigned& k_type) const
  {
   unsigned n_u_types = nu_type_at_each_node();
   unsigned nodal_type_index = alpha*n_u_types + k_type;
   return raw_nodal_value(j_node, nodal_type_index);
  }

  /// Interface to retrieve the t-th history value of
  /// u_alpha at node j of type k
  virtual double get_u_alpha_value_at_node_of_type(const unsigned& t_time,
						   const unsigned& alpha,
						   const unsigned& j_node,
						   const unsigned& k_type) const
  {
   unsigned n_u_types = nu_type_at_each_node();
   unsigned nodal_type_index = alpha*n_u_types + k_type;
   return raw_nodal_value(t_time, j_node, nodal_type_index);
  }
  
  /// Interface to retrieve the value of w at node j of type k
  virtual double get_w_value_at_node_of_type(const unsigned& j_node,
					     const unsigned& k_type) const
  {
   unsigned n_u_fields = 2;
   unsigned n_u_types = nu_type_at_each_node();
   unsigned nodal_type_index = n_u_fields * n_u_types + k_type;
   return raw_nodal_value(j_node, nodal_type_index);
  }
  
  /// Interface to retrieve the t-th history value of w at node j
  /// of type k
  virtual double get_w_value_at_node_of_type(const unsigned& t_time,
					     const unsigned& j_node,
					     const unsigned& k_type) const
  {
   unsigned n_u_fields = 2;
   unsigned n_u_types = nu_type_at_each_node();
   unsigned nodal_type_index = n_u_fields * n_u_types + k_type;
   return raw_nodal_value(t_time, j_node, nodal_type_index);
  }

  // [IN-PLANE-INTERNAL] 
  // /// (pure virtual) interface to get the pointer to the internal data used to
  // /// interpolate u (NOTE: assumes each u field has exactly one internal data)
  // virtual Vector<Data*> u_internal_data_pts()
  //  {
  //   // Write me
  //  }

  /// Interface to get the pointer to the internal data used to
  /// interpolate w (NOTE: assumes w field has exactly one internal data)
  virtual Data* w_internal_data_pt() const
  {
   unsigned i_field = w_field_index();
   unsigned index = CurvableBellElement<NNODE_1D>::index_of_internal_data_for_field(i_field);
   return internal_data_pt(index);
  }

  
  /// Interface to get the number of internal types for the u
  /// fields (left in tact to ensure that this always returns zero)
  virtual unsigned nu_type_internal() const
  { return CurvableBellElement<NNODE_1D>::ninternal_basis_type_for_field(u_field_indices()[0]); }

  /// Interface to get the number of internal types for the w
  /// fields
  virtual unsigned nw_type_internal() const
  { return CurvableBellElement<NNODE_1D>::ninternal_basis_type_for_field(w_field_index()); }

  
  // [IN-PLANE-INTERNAL] 
  // /// (pure virtual) interface to retrieve the value of u_alpha of internal
  // /// type k
  // virtual double get_u_alpha_internal_value_of_type(const unsigned& alpha,
  // 						    const unsigned& k_type) const = 0;

  // /// (pure virtual) interface to retrieve the t-th history value of u_alpha of
  // /// internal type k
  // virtual double get_u_alpha_internal_value_of_type(const unsigned& time,
  // 						    const unsigned& alpha,
  // 						    const unsigned& k_type) const = 0;

  /// Interface to retrieve the value of w of internal type k
  virtual double get_w_internal_value_of_type(const unsigned& k_type) const
  {
   unsigned index = w_field_index();
   return CurvableBellElement<NNODE_1D>::internal_value_for_field_of_type(index, k_type);
  }
  
  /// Interface to retrieve the t-th history value of w of
  /// internal type k
  virtual double get_w_internal_value_of_type(const unsigned& t_time,
					      const unsigned& k_type) const
  {
   unsigned index = w_field_index();
   return CurvableBellElement<NNODE_1D>::internal_value_for_field_of_type(t_time, index, k_type);
  }
  

 protected:
    
  /// In-plane basis functions and derivatives w.r.t. global
  /// coords at local coordinate s; return det(Jacobian of mapping)
  virtual double dbasis_u_eulerian_foeppl_von_karman(const Vector<double> &s,
						     Shape& psi_n,
						     DShape& dpsi_n_dx) const;
    
  /// In-plane basis/test functions at and derivatives w.r.t
  /// global coords at local coordinate s; return det(Jacobian of mapping)
  virtual double dbasis_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double> &s,
							       Shape& psi_n,
							       DShape& dpsi_n_dx,
							       Shape& test_n,
							       DShape& dtest_n_dx) const;

  /// Out-of-plane basis/test functions at local coordinate s
  virtual void basis_and_test_w_foeppl_von_karman(const Vector<double> &s,
						  Shape& psi_n,
						  Shape& psi_i,
						  Shape& test_n,
						  Shape& test_i) const;
   
  /// Out-of-plane basis/test functions and first derivs w.r.t.
  /// to global coords at local coordinate s; return det(Jacobian of mapping)
  virtual double dbasis_and_dtest_w_eulerian_foeppl_von_karman(const Vector<double> &s,
							       Shape& psi_n,
							       Shape& psi_i,
							       DShape& dpsi_n_dx,
							       DShape& dpsi_i_dx,
							       Shape& test_n,
							       Shape& test_i,
							       DShape& dtest_n_dx,
							       DShape& dtest_i_dx) const;

  /// Out-of-plane basis functions and derivs w.r.t. global
  /// coords at local coordinate s; return det(Jacobian of mapping)
  virtual double d2basis_w_eulerian_foeppl_von_karman(const Vector<double>& s,
						      Shape& psi_n,
						      Shape& psi_i,
						      DShape& dpsi_n_dx,
						      DShape& dpsi_i_dx,
						      DShape& d2psi_n_dx2,
						      DShape& d2psi_i_dx2) const;
  
  /// Out-of-plane basis/test functions and first/second derivs
  /// w.r.t. to global coords at local coordinate s;
  /// return det(Jacobian of mapping)
  virtual double d2basis_and_d2test_w_eulerian_foeppl_von_karman(const Vector<double>& s,
								 Shape& psi_n,
								 Shape& psi_i,
								 DShape& dpsi_n_dx,
								 DShape& dpsi_i_dx,
								 DShape& d2psi_n_dx2,
								 DShape& d2psi_i_dx2,
								 Shape& test_n,
								 Shape& test_i,
								 DShape& dtest_n_dx,
								 DShape& dtest_i_dx,
								 DShape& d2test_n_dx2,
								 DShape& d2test_i_dx2) const;

  // End of FoepplVonKarmanEquations interface functions
  //----------------------------------------------------------------------------

 public:

  /// Function to pin all deflection dofs
  void pin_all_deflection_dofs() const;

  /// Function to pin particular in--plane displacement dofs
  void fix_in_plane_displacement_dof(const unsigned& dof_number,
				     const unsigned& b,
				     const ScalarFctPt& u);

  /// Function to pin particular out--of--plane displacement dofs
  void fix_out_of_plane_displacement_dof(const unsigned& dof_number,
					 const unsigned& b,
					 const ScalarFctPt& w);

  // Is this element curved?
  bool element_is_curved() const
  { return CurvableBellElement<NNODE_1D>::element_is_curved(); }
  
 
 protected:
  /// Get rotation matrices that change the degrees of freedom to the basis set
  /// by Rotated_basis_fct_pt
  inline void rotation_matrix_at_node(const unsigned& inode, DenseDoubleMatrix&
				      rotation_matrix) const;
 
  /// Transform the shape functions so that they correspond to
  /// the new rotated dofs
  inline void rotate_shape(Shape& shape) const;

  /// Transform the shape functions and first derivatives so that they 
  /// correspond to the new rotated dofs
  inline void rotate_shape(Shape& shape, DShape& dshape) const;
  
  /// Transform the shape functions, first and second   derivatives so 
  /// that they correspond to the new rotated dofs
  inline void rotate_shape(Shape& shape, DShape& dshape,
			   DShape& d2shape) const;
  
 public:

  /// Set up the rotated degrees of freedom
  inline void set_up_rotated_dofs(const unsigned& nnodes_to_rotate,
				  const Vector<unsigned>& nodes_to_rotate,
				  const BasisVectorsFctPt& basis_vectors_fct_pt);

  /// Upgrade the Bell element to a curved Bernadou element
  virtual void upgrade_element_to_curved(const MyC1CurvedElements::Edge& curved_edge, const double& s_ubar,
					 const double& s_obar,  CurvilineGeomObject* parametric_edge, 
					 const unsigned& boundary_order)
  {
   CurvableBellElement<NNODE_1D>::
    upgrade_element_to_curved(curved_edge, s_ubar, s_obar,
			      parametric_edge, boundary_order);
  }
  
  /// Required  # of `values' (pinned or dofs)
  /// at node n
  inline unsigned required_nvalue(const unsigned &n) const
  {
   return Initial_Nvalue[n];
  }

  /// Output function:
  ///  x, y, ux, uy, w
  void output(std::ostream &outfile)
  {
   FoepplVonKarmanEquations::output(outfile);
  }
 
  /// Full output function with a rich set of unknowns:
  ///  x, y, ux, uy, w, dw, ddw, du, strain, stress, principal stress
  void full_output(std::ostream &outfile)
  {
   FoepplVonKarmanEquations::full_output(outfile);
  }

  /// Output function:
  ///   x,y,u   or    x,y,z,u at n_plot*(n_plot+1)/2 plot points
  void output(std::ostream &outfile, const unsigned &n_plot)
  {FoepplVonKarmanEquations::output(outfile,n_plot);}

  /// Output function:
  ///   x,y,u   or    x,y,z,u at n_plot*(n_plot+1)/2plot points
  void output_interpolated_exact_soln(std::ostream &outfile,
				      FiniteElement::SteadyExactSolutionFctPt exact_soln_pt, const unsigned& n_plot);

  /// C-style output function:
  ///  x,y,u   or    x,y,z,u
  void output(FILE* file_pt)
  {FoepplVonKarmanEquations::output(file_pt);}

  /// C-style output function:
  ///   x,y,u   or    x,y,z,u at n_plot*(n_plot+1)/2 plot points
  void output(FILE* file_pt, const unsigned &n_plot)
  {FoepplVonKarmanEquations::output(file_pt,n_plot);}


  /// Output function for an exact solution:
  ///  x,y,u_exact   or    x,y,z,u_exact at n_plot*(n_plot+1)/2 plot points
  void output_fct(std::ostream &outfile, const unsigned &n_plot,
		  FiniteElement::SteadyExactSolutionFctPt exact_soln_pt)
  {
   FoepplVonKarmanEquations::output_fct(outfile,n_plot,
					exact_soln_pt);
  }

  /// Output function for a time-dependent exact solution.
  ///  x,y,u_exact   or    x,y,z,u_exact at n_plot*(n_plot+1)/2 points
  /// (Calls the steady version)
  void output_fct(std::ostream &outfile, const unsigned &n_plot,
		  const double& time,
		  FiniteElement::UnsteadyExactSolutionFctPt exact_soln_pt)
  {
   FoepplVonKarmanEquations::output_fct(outfile,n_plot,time,exact_soln_pt);
  }


  
  // Private Data Members
  
 private:
  
  /// Static number of fields (is always 3)
  static const unsigned Nfield;

  /// Static bool vector with the Bell interpolation of the fields
  /// (always only w)
  static const std::vector<bool> Field_is_bell_interpolated;
  
  /// Static int that holds the number of variables at
  /// nodes: always the same
  static const unsigned Initial_Nvalue[];

  /// Enum to store which edge is curved set to none when element has no curved
  /// edges
  MyC1CurvedElements::Edge Curved_edge;

  /// Pointer to Stored Association matrix
  DenseMatrix<double>* Association_matrix_pt;

  /// A Pointer to the function that sets up the rotated basis at point x
  BasisVectorsFctPt Rotated_basis_fct_pt;

  /// Which nodes are we rotating
  Vector<unsigned> Nodes_to_rotate;

  /// Number of nodes to rotate
  unsigned Nnodes_to_rotate;


 };




 ////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////
 ////////////////////////////////////////////////////////////////////


 //==============================================================================
 /// Face geometry for the FoepplVonKarmanC1CurvableBellElement elements: The
 /// spatial dimension of the face elements is one lower than that of the bulk
 /// element but they have the same number of points along their 1D edges.
 //==============================================================================
 template < unsigned NNODE_1D>
 class FaceGeometry<FoepplVonKarmanC1CurvableBellElement<NNODE_1D> >:
  public virtual TElement<1,NNODE_1D>
 {

 public:

  /// \short Constructor: Call the constructor for the
  /// appropriate lower-dimensional TElement
  FaceGeometry() : TElement<1,NNODE_1D>() {}

 };



 //==============================================================================
 /// Set up the rotated degrees of freedom: includes a check for the number of
 /// rotation nodes being greater than three.
 //==============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::set_up_rotated_dofs(const unsigned&
									nnodes_to_rotate, const Vector<unsigned>& nodes_to_rotate, const
									BasisVectorsFctPt& basis_vectors_fct_pt)
 {
  // Change the member Nnode_to_rotate
  Nnodes_to_rotate=nnodes_to_rotate;
#ifdef PARANOID
  // Check that the number of nodes is smaller than 3
  if( nnodes_to_rotate > 3 )
   {
    throw OomphLibError(
			"There are only three nodes per element, so we cannot rotate more than\
 three ", OOMPH_CURRENT_FUNCTION,OOMPH_EXCEPTION_LOCATION);
   }
#endif

  Nodes_to_rotate = nodes_to_rotate;

  // Point to the basis vectors function
  Rotated_basis_fct_pt = basis_vectors_fct_pt;
 }


 
 //==============================================================================
 /// Rotate the shape functions according to
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //==============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 rotation_matrix_at_node (const unsigned& inode, DenseDoubleMatrix&
			  rotation_matrix) const
 {
  // Initialise x normal and tangent
  Vector<double> x(2,0.0);

  // Get the node pointer
  Node* nod_pt=this->node_pt(inode);

  // Get the position of the vertex
  x[0]=nod_pt->x(0);
  x[1]=nod_pt->x(1);

  // Initialise the two basis vectors
  Vector<Vector<double> > bi(2,Vector<double>(2,0.0));
  Vector<DenseMatrix<double> > Dbi(2,DenseMatrix<double>(2,2,0.0));

  // Now find the two new basis vectors
  // Get the normal - overload if not continuous
  (*Rotated_basis_fct_pt)(x,bi[0],bi[1], Dbi[0], Dbi[1]);

  // Rotation matrix, B
  DenseMatrix<double> b1(2,2,0.0),b22(3,3,0.0), b21(3,2,0.0);

  // Fill in the submatrices
  for(unsigned alpha=0; alpha<2; ++alpha)
   {
    for(unsigned beta=0; beta<2; ++beta)
     {
      // Fill in b1 - the Jacobian
      // Fill in the rotation of the first derivatives
      b1(alpha,beta) = bi[beta][alpha];
 
      // Avoid double counting the cross derivative
      if(alpha<=beta)
       {
	// Define row index
	const unsigned row = alpha + beta;
	for(unsigned gamma=0; gamma<2; ++gamma)
	 {
	  // Fill in b21 - the non affine part of the Jacobian derivative
	  // Define column index
	  unsigned col_b21 = gamma;
	  // Fill in the non-affine part of the rotation of the second derivatives
	  // if( beta>= alpha) ?
	  b21(row,col_b21) += Dbi[gamma](alpha,beta);
	  for(unsigned delta=0; delta<2; ++delta)
	   {
      
	    // Fill in b22 - the Affine part of the Jacobian derivative
	    // Redefine column index for the next submatrix
	    unsigned col_b22 = gamma + delta;
	    // Fill in the affine part of the rotation of the second derivatives
	    // if( beta>= alpha) ?
	    b22(row,col_b22) += bi[gamma][alpha] * bi[delta][beta];
	   }
	 }
       }
     }
   }


  // Fill in the submatrices to the full (6x6) matrix - we need to right multiply
  // this matrix so we need the transpose of the Jacobian
  // W dof remains the same
  rotation_matrix(0,0)=1.0;
  // Fill in b1
  for(unsigned i=0; i<2; ++i){
   for(unsigned j=0; j<2; ++j)
    { rotation_matrix(1+j,1+i) = b1(i,j); }
  }
  // Fill in b21
  for(unsigned i=0; i<3; ++i){
   for(unsigned j=0; j<2; ++j)
    { rotation_matrix(1+j,3+i) = b21(i,j); }
  }
  // Fill in b22
  for(unsigned i=0; i<3; ++i){
   for(unsigned j=0; j<3; ++j)
    { rotation_matrix(3+j,3+i) = b22(i,j); }
  }
 }


 //============================================================================
 /// Fetch the in-plane basis functions and their derivatives w.r.t. global
 /// coordinates at s and return Jacobian of mapping.
 //=============================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 dbasis_u_eulerian_foeppl_von_karman(const Vector<double>& s,
				     Shape& psi_n,
				     DShape& dpsi_n_dx) const
 { 
  // Initialise and get dpsi w.r.t local coord
  const unsigned n_node = this->nnode();
  DShape dpsids(n_node,this->dim());
  dshape_u_local(s,psi_n,dpsids);
  double J;

  // Get the Jacobian of the mapping
  DenseMatrix<double> jacobian(this->dim(),this->dim()),
   inverse_jacobian(this->dim(),this->dim());
  J = CurvableBellElement<NNODE_1D>::local_to_eulerian_mapping(s,jacobian,inverse_jacobian);

  // Now find the global derivatives
  for(unsigned l=0;l<n_node;++l)
   { 
    for(unsigned i=0; i<this->dim(); ++i)
     {
      // Initialise to zero
      dpsi_n_dx(l,i)=0.0;
      for(unsigned j=0; j<this->dim(); ++j)
       {
	// Convert to local coordinates
	dpsi_n_dx(l,i)+=inverse_jacobian(i,j)*dpsids(l,j);
       }
     }
   }
  return J;
 }

 
 //============================================================================
 /// Fetch the in-plane basis functions and test functions and derivatives
 /// w.r.t. global coordinates at s and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //=============================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 dbasis_and_dtest_u_eulerian_foeppl_von_karman(const Vector<double>& s,
					       Shape& psi_n,
					       DShape& dpsi_n_dx,
					       Shape& test_n,
					       DShape& dtest_n_dx) const
 {
  // Initialise and get dpsi w.r.t local coord
  const unsigned n_node = this->nnode();
  DShape dpsids(n_node,this->dim());
  dshape_u_local(s,psi_n,dpsids);
  double J;

  // Get the Jacobian of the mapping
  DenseMatrix<double> jacobian(this->dim(),this->dim()),
   inverse_jacobian(this->dim(),this->dim());
  J = CurvableBellElement<NNODE_1D>::local_to_eulerian_mapping(s,jacobian,inverse_jacobian);

  // Now find the global derivatives
  for(unsigned l=0;l<n_node;++l)
   { 
    for(unsigned i=0; i<this->dim(); ++i)
     {
      // Initialise to zero
      dpsi_n_dx(l,i)=0.0;
      for(unsigned j=0; j<this->dim(); ++j)
       {
	// Convert to local coordinates
	dpsi_n_dx(l,i)+=inverse_jacobian(i,j)*dpsids(l,j);
       }
     }
   }
  test_n = psi_n;
  dtest_n_dx = dpsi_n_dx;
  return J;
 }

 
 //======================================================================
 /// Define the shape functions and test functions and derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //======================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 basis_and_test_w_foeppl_von_karman(const Vector<double> &s,
				    Shape& psi_n,
				    Shape& psi_i,
				    Shape& test_n,
				    Shape& test_i) const
 {
  throw OomphLibError(
		      "This still needs testing for curved elements.",
		      "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::\
shape_and_test_foeppl_von_karman(...)", OOMPH_EXCEPTION_LOCATION);

  this->c1_basis(s,psi_n,psi_i); 
 
  // Rotate the degrees of freedom
  rotate_shape(psi_n);
 
  // Galerkin
  // (Shallow) copy the basis functions
  test_n = psi_n;
  test_i = psi_i;
 }


 //======================================================================
 /// Fetch the basis functions and test functions and first derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //======================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 dbasis_and_dtest_w_eulerian_foeppl_von_karman(const Vector<double> &s,
					       Shape& psi_n,
					       Shape& psi_i,
					       DShape& dpsi_n_dx,
					       DShape& dpsi_i_dx,
					       Shape& test_n,
					       Shape& test_i,
					       DShape& dtest_n_dx,
					       DShape& dtest_i_dx) const
 {
  // Throw if called 
  throw OomphLibError(
		      "This still needs testing for curved elements.",
		      "void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::\
dshape_and_dtest_foeppl_von_karman(...)",
		      OOMPH_EXCEPTION_LOCATION);// HERE

  // Get the basis 
  double J=this->d_c1_basis_eulerian(s,psi_n,psi_i,dpsi_n_dx,dpsi_i_dx);
 
  // Rotate the degrees of freedom
  rotate_shape(psi_n, dpsi_n_dx);
  // Galerkin
  // (Shallow) copy the basis functions
  test_n = psi_n;
  dtest_n_dx= dpsi_n_dx;
  test_i = psi_i;
  dtest_i_dx= dpsi_i_dx;
  
  return J;
 }



 //======================================================================
 /// Fetch the basis functions and first/second derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //======================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 d2basis_w_eulerian_foeppl_von_karman(const Vector<double>& s,
						 Shape& psi_n,
						 Shape& psi_i,
						 DShape& dpsi_n_dx,
						 DShape& dpsi_i_dx,
						 DShape& d2psi_n_dx,
						 DShape& d2psi_i_dx) const
 {
  //Call the geometrical shape functions and derivatives
  double J=CurvableBellElement<NNODE_1D>::d2_c1_basis_eulerian(s,psi_n,psi_i,dpsi_n_dx,dpsi_i_dx,d2psi_n_dx,d2psi_i_dx);
  // Rotate the dofs
  rotate_shape(psi_n, dpsi_n_dx, d2psi_n_dx);
 
  //Return the jacobian
  return J;
 }

 

 
 //======================================================================
 /// Fetch the basis functions and test functions and first/second derivatives
 /// w.r.t. global coordinates and return Jacobian of mapping.
 ///
 /// Galerkin: Test functions = shape functions
 //======================================================================
 template < unsigned NNODE_1D>
 double FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 d2basis_and_d2test_w_eulerian_foeppl_von_karman(const Vector<double>& s,
						 Shape& psi_n,
						 Shape& psi_i,
						 DShape& dpsi_n_dx,
						 DShape& dpsi_i_dx,
						 DShape& d2psi_n_dx,
						 DShape& d2psi_i_dx,
						 Shape& test_n,
						 Shape& test_i,
						 DShape& dtest_n_dx,
						 DShape& dtest_i_dx,
						 DShape& d2test_n_dx,
						 DShape& d2test_i_dx) const
  {
  //Call the geometrical shape functions and derivatives
  double J=CurvableBellElement<NNODE_1D>::d2_c1_basis_eulerian(s,
							       psi_n,psi_i,
							       dpsi_n_dx,dpsi_i_dx,
							       d2psi_n_dx,d2psi_i_dx);
  // Rotate the dofs
  rotate_shape(psi_n, dpsi_n_dx, d2psi_n_dx);
 
  // Galerkin
  //Set the test functions equal to the shape functions (this is a shallow copy)
  test_n = psi_n;
  dtest_n_dx = dpsi_n_dx;
  d2test_n_dx = d2psi_n_dx;
  test_i = psi_i;
  dtest_i_dx = dpsi_i_dx;
  d2test_i_dx = d2psi_i_dx;

  //Return the jacobian
  return J;
 }


 
 //======================================================================
 /// Rotate the shape functions according to the specified basis on the 
 /// boundary. This function does a DenseDoubleMatrix solve to determine 
 /// new basis - which could be speeded up by caching the matrices higher
 /// up and performing the LU decomposition only once
 //======================================================================
 template < unsigned NNODE_1D>
 inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 rotate_shape(Shape& psi) const
 {
  const unsigned n_dof_types  = nw_type_at_each_node();
  // Loop over the nodes with rotated dofs
  for(unsigned i=0; i<Nnodes_to_rotate; ++i)
   {
    // Get the nodes
    unsigned inode = Nodes_to_rotate[i];

    // Construct the vectors to hold the shape functions
    Vector<double> psi_vector(n_dof_types);

    // Get the rotation matrix
    DenseDoubleMatrix  rotation_matrix(n_dof_types,n_dof_types,0.0);
    this->rotation_matrix_at_node(inode,rotation_matrix);

    // Copy to the vectors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy to the vectors
      for(unsigned k=0;k<n_dof_types;++k)
       {
	// Copy over shape functions
	// psi_vector[l]=psi(inode,l);
	psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
       }
     }

    // Copy back to shape the rotated vetcors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
     }
   }
 }

 
 //======================================================================
 /// Rotate the shape functions according to the specified basis on the 
 /// boundary. This function does a DenseDoubleMatrix solve to determine 
 /// new basis - which could be speeded up by caching the matrices higher
 /// up and performing the LU decomposition only once
 //======================================================================
 template < unsigned NNODE_1D>
 inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 rotate_shape(Shape& psi, DShape& dpsidx) const
 {
  const unsigned n_dof_types  = nw_type_at_each_node();
  const unsigned n_dim  = this->dim();
  // Loop over the nodes with rotated dofs
  for(unsigned i=0; i<Nnodes_to_rotate; ++i)
   {
    // Get the nodes
    unsigned inode = Nodes_to_rotate[i];

    // Construct the vectors to hold the shape functions
    Vector<double> psi_vector(n_dof_types);
    Vector<Vector<double> > dpsi_vector_dxi(n_dim,Vector<double>(n_dof_types));

    // Get the rotation matrix
    DenseDoubleMatrix  rotation_matrix(n_dof_types,n_dof_types,0.0);
    this->rotation_matrix_at_node(inode,rotation_matrix);

    // Copy to the vectors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy to the vectors
      for(unsigned k=0;k<n_dof_types;++k)
       {
	// Copy over shape functions
	// psi_vector[l]=psi(inode,l);
	psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
	// Copy over first derivatives
	for(unsigned i=0;i<n_dim; ++i)
	 {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
       }
     }

    // Copy back to shape the rotated vetcors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<n_dim; ++i)
       {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
     }
   }
 }

 //======================================================================
 /// Rotate the shape functions according to the specified basis on the 
 /// boundary. This function does a DenseDoubleMatrix solve to determine 
 /// new basis - which could be speeded up by caching the matrices higher
 /// up and performing the LU decomposition only once
 //======================================================================
 template < unsigned NNODE_1D>
 inline void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 rotate_shape(Shape& psi, DShape& dpsidx, DShape& d2psidx) const
 {
  const unsigned n_dof_types  = nw_type_at_each_node();
  const unsigned n_dim  = this->dim();
  // n_dimth triangle number
  const unsigned n_2ndderiv = ((n_dim+1)*(n_dim))/2; // Guaranteed to be positive
  // Loop over the nodes with rotated dofs
  for(unsigned i=0; i<Nnodes_to_rotate; ++i)
   {
    // Get the nodes
    unsigned inode = Nodes_to_rotate[i];

    // Construct the vectors to hold the shape functions
    Vector<double> psi_vector(n_dof_types);
    Vector<Vector<double> > dpsi_vector_dxi(n_dim,Vector<double>(n_dof_types));
    Vector<Vector<double> > d2psi_vector_dxidxj(n_2ndderiv,Vector<double>(n_dof_types));

    // Get the rotation matrix
    DenseDoubleMatrix  rotation_matrix(n_dof_types,n_dof_types,0.0);
    this->rotation_matrix_at_node(inode,rotation_matrix);

    // Copy to the vectors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy to the vectors
      for(unsigned k=0;k<n_dof_types;++k)
       {
	// Copy over shape functions
	// psi_vector[l]=psi(inode,l);
	psi_vector[l]+=psi(inode,k)*rotation_matrix(l,k); 
	// Copy over first derivatives
	for(unsigned i=0;i<n_dim; ++i)
	 {dpsi_vector_dxi[i][l]+=dpsidx(inode,k,i)*rotation_matrix(l,k);}
	for(unsigned i=0;i<n_2ndderiv; ++i)
	 {d2psi_vector_dxidxj[i][l]+=d2psidx(inode,k,i)*rotation_matrix(l,k);}
       }
     }

    // Copy back to shape the rotated vetcors
    for(unsigned l=0;l<n_dof_types;++l)
     {
      // Copy over shape functions
      psi(inode,l)=psi_vector[l];
      // Copy over first derivatives
      for(unsigned i=0;i<n_dim; ++i)
       {dpsidx(inode,l,i)=dpsi_vector_dxi[i][l];}
      // Copy over second derivatives
      for(unsigned i=0;i<n_2ndderiv; ++i)
       {d2psidx(inode,l,i)=d2psi_vector_dxidxj[i][l];}
     }
   }
 }


 
 //=============================================================================
 /// Function to pin all deflection dofs
 //=============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::pin_all_deflection_dofs() const
 {
  const unsigned w_index = w_field_index();
    
  // Get nodes
  const unsigned n_node = nw_node();
  const Vector<unsigned> nodes = get_w_node_indices();

  // Curved Bell elements only have deflection dofs at vertices
  for(unsigned j_nodei=0; j_nodei<n_node; j_nodei++)
   {
    // Get the j_nodei-th node used by i_field
    unsigned j_node = nodes[j_nodei];
    Node* nod_pt=this->node_pt(j_node);

    // Get the number of types at the current node
    unsigned n_type = nw_type_at_each_node();
    
    // Check if it is on the boundary
    for(unsigned k_type=0; k_type<n_type; k_type++)
     {
      // Pin and set the value
      nod_pt->pin(2+k_type);
      nod_pt->set_value(2+k_type,0.0);
     }
   }

  // Now fix internal dofs
  unsigned n_internal = nw_type_internal();
  for(unsigned k_type=0; k_type<n_internal; k_type++)
   {
    // Get node
    // Pin and set the value
    this->internal_data_for_field_pt(w_index)->pin(k_type);
   }
 }



 //=============================================================================
 /// Function to pin particular out-of-plane displacement dofs
 //=============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 fix_out_of_plane_displacement_dof(const unsigned& j_type,
				   const unsigned& b_boundary,
				   const ScalarFctPt& specified_deflection_fct_pt)
 {
  const unsigned w_index = w_field_index();
  const unsigned first_nodal_type_index = this->first_nodal_type_index_for_field(w_index);
  const unsigned n_vertices = nw_node();
  const unsigned n_type = nw_type_at_each_node();

#ifdef PARANOID
  // Check that the dof number is a sensible value
  if(j_type >= n_type)
   {
    throw OomphLibError("Foppl von Karman elements only have 6 Hermite deflection degrees\
of freedom at internal points. They are {w ; w,x ; w,y ; w,xx ; w,xy ; w,yy}",
                        "FoepplVonKarmanC1CurvableBellElement:fix_out_of_plane_dof()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
  // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<n_vertices; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(b_boundary);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      Vector<double> x(2);
      x[0]=nod_pt->x(0);
       x[1]=nod_pt->x(1);
      // Get value
      double value;
      specified_deflection_fct_pt(x,value);
      // Pin and set the value
      nod_pt->pin(first_nodal_type_index+j_type);
      nod_pt->set_value(first_nodal_type_index+j_type,value);
     }
   }
 }



 
 //=============================================================================
 /// Function to pin particular out-of-plane displacement dofs
 //=============================================================================
 template < unsigned NNODE_1D>
 void FoepplVonKarmanC1CurvableBellElement<NNODE_1D>::
 fix_in_plane_displacement_dof(const unsigned& alpha,
			       const unsigned& b,
			       const ScalarFctPt& specified_displacement_fct_pt)
 {
  // Initialise constants that we use in this function
  const unsigned field_index = u_alpha_field_index(alpha);
  const unsigned nodal_type_index = this->first_nodal_type_index_for_field(field_index);
  const unsigned n_node = nu_node();
  const unsigned n_type = 2*nu_type_at_each_node();
  const unsigned dim = this->dim();

#ifdef PARANOID
  // Check that the dof number is a sensible value
  if(alpha >= n_type)
   {
    throw OomphLibError("Foppl von Karman elements only have 2 in-plane displacement degrees\
of freedom at internal points. They are {ux, uy}",
                        "FoepplVonKarmanC1CurvableBellElement:fix_out_of_plane_dof()",
                        OOMPH_EXCEPTION_LOCATION);
   }
#endif
  
   // Bell elements only have deflection dofs at vertices
  for(unsigned n=0; n<n_node; ++n)
   {
    // Get node
    Node* nod_pt=this->node_pt(n);
    // Check if it is on the boundary
    bool is_boundary_node=nod_pt->is_on_boundary(b);
    if(is_boundary_node)
     {
      // Extract nodal coordinates from node:
      // Since the element isn't necessarily isoparametric the nodes 'position'
      // is not necessarily correct?
      Vector<double> x(dim),s(dim);
      this->local_coordinate_of_node(n,s);
      interpolated_x(s, x); 
      // Fill in value
      double value;
      specified_displacement_fct_pt(x,value);
      // Pin and set the value
      nod_pt->pin(alpha);
      nod_pt->set_value(nodal_type_index, value);
     }
   }
 }
}


#endif
