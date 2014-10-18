#include "BonStdCInterface.h"
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <limits.h>

/* Function Declarations */
Bool eval_f(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data);

Bool eval_grad_f(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data);

Bool eval_g(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data);

Bool eval_jac_g(Index n, Number *x, Bool new_x, Index m, Index nnz_jac_g, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data);


Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor, Index m, Number *lambda, Bool new_lambda,
	Index nnz_h_lag, Index *iRow, Index *jCol, Number *values, UserDataPtr user_data);


Bool intermediate_cb(Index alg_mod, Index iter_count, Number obj_value,
	Number inf_pr, Number inf_du, Number mu, Number d_norm,
	Number regularization_size, Number alpha_du,
	Number alpha_pr, Index ls_trials, UserDataPtr user_data);

struct MyUserData
{
	// nothing
};

/* Main Program */
int main()
{
	struct MyUserData user_data;
	Index n = 4;
	Index m = 3
	VariableTypeC* var_types = ( VariableTypeC* )malloc( sizeof( VariableTypeC )*n );
	var_types[0] = BINARY;
	var_types[1] = CONTINUOUS;
	var_types[2] = CONTINUOUS;
	var_types[3] = INTEGER;

	LinearityTypeC* var_linearity_types = ( LinearityTypeC* )malloc( sizeof( LinearityTypeC )*n );
	var_linearity_types[0] = LINEAR;
	var_linearity_types[1] = NON_LINEAR;
	var_linearity_types[2] = NON_LINEAR;
	var_linearity_types[3] = LINEAR;

	LinearityTypeC* constraint_linearity_types = ( LinearityTypeC* )malloc( sizeof( LinearityTypeC )*m );
	constraint_linearity_types[0] = NON_LINEAR;
	constraint_linearity_types[1] = LINEAR;
	constraint_linearity_types[2] = LINEAR;

	Index nnz_jac_g = 7;
	Index nnz_h_lag = 2;

	Number* x_L = (Number*)malloc(sizeof(Number)*n);
	Number* x_U = (Number*)malloc(sizeof(Number)*n);
	Number* g_L = (Number*)malloc(sizeof(Number)*m);
	Number* g_U = (Number*)malloc(sizeof(Number)*m);

	x_l[0] = 0.;
	x_u[0] = 1.;

	x_l[1] = 0.;
	x_u[1] = DBL_MAX;

	x_l[2] =0.;
	x_u[2] = DBL_MAX;

	x_l[3] = 0;
	x_u[3] = 5;

	g_l[0] = -DBL_MAX;
	g_u[0] = 1./4.;

	g_l[1] = -DBL_MAX;
	g_u[1] = 0;

	g_l[2] = -DBL_MAX;
	g_u[2] = 2;

	Number* x = (Number*)malloc(sizeof(Number)*n);
	x[0] = 0.0;
	x[1] = 0.0;
	x[2] = 0.0;
	x[3] = 0.0;
	BonminSosInfo* sos_info = NULL;
	BonminBrancingInfo* branch_info = NULL;
	BonminProblem minlp = CreateBonminProblem(  n, x_L, x_U,
		m, g_L, g_U,
		nnz_jac_g, nnz_h_lag,
		&eval_f,
		&eval_g,
		&eval_grad_f,
		&eval_jac_g,
		&eval_h,
		var_types,
		var_linearity_types,
		constraint_linearity_types,
		sos_info,
		branch_info);

	AddBonminNumOption( minlp, "bonmin.time_limit", 5 );
	AddBonminStrOption( minlp, "mu_oracle", "loqo" );
	ReadBonminOptFile( minlp, "MyBonmin.opt" );
	ReadBonminOptFile( minlp );
	AddBonminStrOption( minlp, "bonmin.algorithm B-BB\n");

	free(x_L);
	free(x_U);
	free(g_L);
	free(g_U);
	free( var_types );
	free( var_linearity_types );
	free( constraint_linearity_types );

	Number final_obj_value;
	Int status = BonminSolve( minlp, x, NULL, &final_obj_value, mult_g, mult_x_L, mult_x_U, &user_data );

	FreeBonminProblem( minlp );
	free( x );
	free(mult_g);
	free(mult_x_L);
	free(mult_x_U);
	return 0;
}

Bool eval_f(Index n, Number* x, Bool new_x, Number* obj_value, UserDataPtr user_data)
{
	assert(n == 4);

	*obj_value = - x[0] - x[1] - x[2];

	return TRUE;
}

Bool eval_grad_f(Index n, Number* x, Bool new_x, Number* grad_f, UserDataPtr user_data)
{
	assert(n == 4);

	grad_f[0] = -1.;
	grad_f[1] = -1.;
	grad_f[2] = -1.;
	grad_f[3] = 0.;

	return TRUE;
}

Bool eval_g(Index n, Number* x, Bool new_x, Index m, Number* g, UserDataPtr user_data)
{
	struct MyUserData* my_data = user_data;

	assert(n == 4);
	assert(m == 3);

	g[0] = (x[1] - 1./2.)*(x[1] - 1./2.) + (x[2] - 1./2.)*(x[2] - 1./2.);
	g[1] = x[0] - x[1];
	g[2] = x[0] + x[2] + x[3];

	return TRUE;
}

Bool eval_jac_g(Index n, Number *x, Bool new_x, Index m, Index nnz_jac_g,
	Index *iRow, Index *jCol, Number *values, UserDataPtr user_data)
{
	assert(n==4);
	assert(nnz_jac == 7);
	if(values == NULL) {
		iRow[0] = 2;
		jCol[0] = 1;

		iRow[1] = 3;
		jCol[1] = 1;

		iRow[2] = 1;
		jCol[2] = 2;

		iRow[3] = 2;
		jCol[3] = 2;

		iRow[4] = 1;
		jCol[4] = 3;

		iRow[5] = 3;
		jCol[5] = 3;

		iRow[6] = 3;
		jCol[6] = 4;
	}
	else
	{
		values[0] = 1.;
		values[1] = 1;

		values[2] = 2*x[1] - 1;
		values[3] = -1.;

		values[4] = 2*x[2] - 1;
		values[5] = 1.;

		values[6] = 1.;
	}

	return TRUE;
}

Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
	Index m, Number *lambda, Bool new_lambda,
	Index nnz_h_lag, Index *iRow, Index *jCol,
	Number *values, UserDataPtr user_data)
{
	assert (n==4);
	assert (m==3);
	assert(nnz_h_lag==2);
	if(values==NULL)
	{
		iRow[0] = 2;
		jCol[0] = 2;

		iRow[1] = 3;
		jCol[1] = 3;
	}
	else {
		values[0] = 2*lambda[0];
		values[1] = 2*lambda[0];
	}
	return TRUE;
}