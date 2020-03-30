// =================================================================================
// Set the attributes of the primary field variables
// =================================================================================
// This function sets attributes for each variable/equation in the app. The
// attributes are set via standardized function calls. The first parameter for each
// function call is the variable index (starting at zero). The first set of
// variable/equation attributes are the variable name (any string), the variable
// type (SCALAR/VECTOR), and the equation type (EXPLICIT_TIME_DEPENDENT/
// TIME_INDEPENDENT/AUXILIARY). The next set of attributes describe the
// dependencies for the governing equation on the values and derivatives of the
// other variables for the value term and gradient term of the RHS and the LHS.
// The final pair of attributes determine whether a variable represents a field
// that can nucleate and whether the value of the field is needed for nucleation
// rate calculations.


// Interpolation function for the parameters of the individual phases
scalarvalueType interp_function(scalarvalueType phi)
{
	return phi * phi * phi * (6 * phi * phi - 15 * phi + 10);
}

scalarvalueType derivative_interp_function(scalarvalueType phi)
{
	return 30 * phi * phi * (phi - 1) * (phi - 1);
}

void variableAttributeLoader::loadVariableAttributes(){
	// Variable 0
	set_variable_name				(0,"c");
	set_variable_type				(0,SCALAR);
	set_variable_equation_type		(0,EXPLICIT_TIME_DEPENDENT);

    set_dependencies_value_term_RHS(0, "c");
    set_dependencies_gradient_term_RHS(0, "grad(mu)");

	// Variable 1
	set_variable_name				(1,"mu");
	set_variable_type				(1,SCALAR);
	set_variable_equation_type		(1,AUXILIARY);

    set_dependencies_value_term_RHS(1, "c, grad(u)");
    set_dependencies_gradient_term_RHS(1, "grad(c)");

	// Variable 2
	set_variable_name				(2,"u");
	set_variable_type				(2,VECTOR);
	set_variable_equation_type		(2,TIME_INDEPENDENT);

	set_dependencies_value_term_RHS(4, "");
    set_dependencies_gradient_term_RHS(4, "c, grad(u)");
    set_dependencies_value_term_LHS(4, "");
    set_dependencies_gradient_term_LHS(4, "c, grad(change(u))");
}

// =============================================================================================
// explicitEquationRHS (needed only if one or more equation is explict time dependent)
// =============================================================================================
// This function calculates the right-hand-side of the explicit time-dependent
// equations for each variable. It takes "variable_list" as an input, which is a list
// of the value and derivatives of each of the variables at a specific quadrature
// point. The (x,y,z) location of that quadrature point is given by "q_point_loc".
// The function outputs two terms to variable_list -- one proportional to the test
// function and one proportional to the gradient of the test function. The index for
// each variable in this list corresponds to the index given at the top of this file.

template <int dim, int degree>
void customPDE<dim,degree>::explicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

// --- Getting the values and derivatives of the model variables ---
scalarvalueType c = variable_list.get_scalar_value(0);
scalargradType mux = variable_list.get_scalar_gradient(1);

// --- Setting the expressions for the terms in the governing equations ---
scalarvalueType eq_c = c;
scalargradType eqx_c = constV(-McV*userInputs.dtValue)*mux;

// --- Submitting the terms for the governing equations ---
variable_list.set_scalar_value_term_RHS(0,eq_c);
variable_list.set_scalar_gradient_term_RHS(0,eqx_c);

}

// =============================================================================================
// nonExplicitEquationRHS (needed only if one or more equation is time independent or auxiliary)
// =============================================================================================
// This function calculates the right-hand-side of all of the equations that are not
// explicit time-dependent equations. It takes "variable_list" as an input, which is
// a list of the value and derivatives of each of the variables at a specific
// quadrature point. The (x,y,z) location of that quadrature point is given by
// "q_point_loc". The function outputs two terms to variable_list -- one proportional
// to the test function and one proportional to the gradient of the test function. The
// index for each variable in this list corresponds to the index given at the top of
// this file.

template <int dim, int degree>
void customPDE<dim,degree>::nonExplicitEquationRHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
				 dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

 // The part of chemical potential
 // --- Getting the values and derivatives of the model variables ---

 scalarvalueType c = variable_list.get_scalar_value(0);
 scalargradType cx = variable_list.get_scalar_gradient(0);

 // --- Setting the expressions for the terms in the governing equations ---

 // The derivative of the local free energy
 scalarvalueType fcV = 0;

 // The array of the local free energy parameters
 constV a [9] = {8.072789087, -81.24549382, 408.0297321, -1244.129167,
 	2444.046270, -3120.635139, 2506.663551, -1151.003178, 230.2006355};

// calculation of the polynomial of the local free energy
for (unsigned int i=0, scalarvalueType mult = 1.0; i < 10; i++)
{
	for (unsigned int j=1; j <= i+2; j++)
	{
		mult = mult * c;
	}
	fcV += mult * a[i];
}

// iterpolation function and it's derivative calculation
scalarvaluetype h_eta = interp_function(c);
scalarvaluetype dh_deta = derivative_interp_function(c);

// difference of stiffnesses
dealii::VectorizedArray<double> delta_CIJ[CIJ_tensor_size][CIJ_tensor_size];
for (unsigned int i=0; i<CIJ_tensor_size; i++){
	  for (unsigned int j=0; j<CIJ_tensor_size; j++){
		  delta_CIJ[i][j] = (CIJ_part[i][j] - CIJ_matr[i][j])*dh_deta;
	  }
}

// effective stiffness in the point
dealii::VectorizedArray<double> CIJ_eff[CIJ_tensor_size][CIJ_tensor_size];
for (unsigned int i=0; i<CIJ_tensor_size; i++){
	  for (unsigned int j=0; j<CIJ_tensor_size; j++){
		  CIJ_eff[i][j] = CIJ_part[i][j] * h_eta + CIJ_matr[i][j]*sum_hV*(1.0 - h_eta);
	  }
}

// calculate transformation strain
dealii::VectorizedArray<double> epsilon_transform;
for (unsigned int i=0; i < dim, i++)
{
	for (unsigned int j=0, j < dim, j++)
	{
		epsilon_transform[i][j] = h_eta * epsilon_0[i][j];
	}
}

// calculate derivative of transformation strain
dealii::VectorizedArray<double> epsilon_transform_derivative;
for (unsigned int i=0; i < dim, i++)
{
	for (unsigned int j=0; j < dim, j++)
	{
		epsilon_transform_derivative[i][j] = dh_deta * epsilon_0[i][j];
	}
}

// calculate the strain from displacement field
dealii::VectorizedArray<double> epsilon[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		epsilon[i][j]= constV(0.5)*(ux[i][j]+ux[j][i]);
	}
}

// compute elastic strain and derivative elastic strain
dealii::VectorizedArray<double> epsilon_el[dim][dim],epsilon_el_der[dim][dim];
for (unsigned int i=0; i<dim; i++)
{
	for (unsigned int j=0; j<dim; j++)
	{
		epsilon_el[i][j] = epsilon[i][j] - epsilon_transform[i][j];
		epsilon_el_der[i]][j] = epsilon[i][j] - epsilon_transform_derivative[i][j];
	}
}

// compute the stress values needed
dealii::VectorizedArray<double> sigma[dim][dim], delta_sigma[dim][dim], sigma_derivative[dim][dim];
computeStress<dim>(delta_CIJ, epsilon_el, delta_sigma);
computeStress<dim>(CIJ_eff, epsilon_el, sigma);
computeStress<dim>(CIJ_eff, epsilon_el_der, sigma_derivative);

// Compute elastic contribution to chemical potential
scalarvalueType fel = 0.0;

for (unsigned int i = 0, i < dim, i++)
{
	for (unsigned int j = 0, j < dim, j++)
	{
		fel += delta_sigma[i][j]*epsilon_el[i][j] + sigma_derivative[i][j]*epsilon_el[i][j] +
		sigma[i][j]*epsilon_el_der[i][j];
	}
}

fel *= 0.5;
 // The terms for the governing equations
 scalarvalueType eq_mu = fcV+fel;
 scalargradType eqx_mu = constV(KcV)*cx;

 // --- Submitting the terms for the governing equations ---

 variable_list.set_scalar_value_term_RHS(1,eq_mu);
 variable_list.set_scalar_gradient_term_RHS(1,eqx_mu);



// elasticity equation
vectorgradType eqx_u;

//compute the term in the equation
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		eqx_u[i][j] = -stress[i][j];
	}
}

// --- Submitting the terms for the governing equations ---

variable_list.set_vector_gradient_term_RHS(0,eqx_u);

}

// =============================================================================================
// equationLHS (needed only if at least one equation is time independent)
// =============================================================================================
// This function calculates the left-hand-side of time-independent equations. It
// takes "variable_list" as an input, which is a list of the value and derivatives of
// each of the variables at a specific quadrature point. The (x,y,z) location of that
// quadrature point is given by "q_point_loc". The function outputs two terms to
// variable_list -- one proportional to the test function and one proportional to the
// gradient of the test function -- for the left-hand-side of the equation. The index
// for each variable in this list corresponds to the index given at the top of this
// file. If there are multiple elliptic equations, conditional statements should be
// sed to ensure that the correct residual is being submitted. The index of the field
// being solved can be accessed by "this->currentFieldIndex".

template <int dim, int degree>
void customPDE<dim,degree>::equationLHS(variableContainer<dim,degree,dealii::VectorizedArray<double> > & variable_list,
		dealii::Point<dim, dealii::VectorizedArray<double> > q_point_loc) const {

dealii::VectorizedArray<double> CIJ_eff[CIJ_tensor_size][CIJ_tensor_size];

// iterpolation function and it's derivative calculation
scalarvaluetype h_eta = interp_function(c);
scalarvaluetype dh_deta = derivative_interp_function(c);

// Calculation of effective stiffness
dealii::VectorizedArray<double> CIJ_eff[CIJ_tensor_size][CIJ_tensor_size];
for (unsigned int i=0; i<CIJ_tensor_size; i++){
	  for (unsigned int j=0; j<CIJ_tensor_size; j++){
		  CIJ_eff[i][j] = CIJ_part[i][j] * h_eta + CIJ_matr[i][j]*sum_hV*(constV(1.0) - h_eta);
	  }
}

// Get change in variable
vectorgradType Dux = variable_list.get_change_in_vector_gradient(2);

//compute strain tensor
dealii::VectorizedArray<double> depsilon[dim][dim];
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		depsilon[i][j]= constV(0.5)*(Dux[i][j]+Dux[j][i]);
	}
}

// compute stress tensor
dealii::VectorizedArray<double> dsigma[dim][dim];
computeStress<dim>(CIJ_eff, depsilon, dsigma);

//compute the term in the governing equation
vectorgradType eqx_Du;
for (unsigned int i=0; i<dim; i++){
	for (unsigned int j=0; j<dim; j++){
		eqx_Du[i][j] = dsigma[i][j];
	}
}

 // --- Submitting the terms for the governing equations ---

variable_list.set_vector_gradient_term_LHS(2,eqx_Du);

}