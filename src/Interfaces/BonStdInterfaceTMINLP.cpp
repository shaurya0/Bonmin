#include "BonStdInterfaceTMINLP.hpp"

using namespace Ipopt;
namespace Bonmin
{
	StdInterfaceTMINLP::StdInterfaceTMINLP(Index n_var,
	                const Number* x_L, const Number* x_U,
	                Index n_con,
	                const Number* g_L, const Number* g_U,
	                Index nele_jac,
	                Index nele_hess,
	                Index index_style,
	                const Number* start_x,
	                const Number* start_lam,
	                const Number* start_z_L,
	                const Number* start_z_U,
	                Eval_F_CB eval_f,
	                Eval_G_CB eval_g,
	                Eval_Grad_F_CB eval_grad_f,
	                Eval_Jac_G_CB eval_jac_g,
	                Eval_H_CB eval_h,
	                VariableType* var_types,
	                TNLP::LinearityType* var_linearity_types,
	                TNLP::LinearityType* constraint_linearity_types,
	                // Intermediate_CB intermediate_cb,
	                Number* x_sol,
	                Number* z_L_sol,
	                Number* z_U_sol,
	                Number* g_sol,
	                Number* lam_sol,
	                Number* obj_sol,
	                UserDataPtr user_data,
	                Number obj_scaling,
	                const Number* x_scaling,
	                const Number* g_scaling)
					:
					TMINLP(),
					n_var_(n_var),
					n_con_(n_con),
					x_L_(x_L),
					x_U_(x_U),
					g_L_(g_L),
					g_U_(g_U),
					nele_jac_(nele_jac),
					nele_hess_(nele_hess),
					index_style_(index_style),
					start_x_(start_x),
					start_lam_(start_lam),
					start_z_L_(start_z_L),
					start_z_U_(start_z_U),
					eval_f_(eval_f),
					eval_g_(eval_g),
					eval_grad_f_(eval_grad_f),
					eval_jac_g_(eval_jac_g),
					eval_h_(eval_h),
					var_types_(var_types),
					var_linearity_types_(var_linearity_types),
					constraint_linearity_types_(constraint_linearity_types),
					sos_info_(NULL),
					branch_info_(NULL),
					// intermediate_cb_(intermediate_cb),
					user_data_(user_data),
					obj_scaling_(obj_scaling),
					x_scaling_(NULL),
					g_scaling_(NULL),
					non_const_x_(NULL),
					x_sol_(x_sol),
					z_L_sol_(z_L_sol),
					z_U_sol_(z_U_sol),
					g_sol_(g_sol),
					lambda_sol_(lam_sol),
					obj_sol_(obj_sol)
					{
						if (x_scaling)
						{
							Number* tmp = new Number[n_var_];
							for (Index i=0; i<n_var_; i++)
							{
								tmp[i] = x_scaling[i];
							}
							x_scaling_ = tmp;
						}
						if (g_scaling)
						{
							Number* tmp = new Number[n_con_];
							for (Index i=0; i<n_con_; i++)
							{
								tmp[i] = g_scaling[i];
							}
							g_scaling_ = tmp;
						}
					}

	StdInterfaceTMINLP::~StdInterfaceTMINLP()
	{
		delete [] non_const_x_;
		delete [] x_scaling_;
		delete [] g_scaling_;
	}



	bool StdInterfaceTMINLP::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g, Index& nnz_h_lag, TNLP::IndexStyleEnum& index_style)
	{
		n = n_var_; // # of variables (variable types have been asserted in the constructor
		m = n_con_; // # of constraints
		nnz_jac_g = nele_jac_; // # of non-zeros in the jacobian
		nnz_h_lag = nele_hess_; // # of non-zeros in the hessian

		index_style = (index_style_ == 0) ? TNLP::C_STYLE : TNLP::FORTRAN_STYLE;

		return true;
	}


	bool StdInterfaceTMINLP::get_bounds_info(Index n, Number* x_l, Number* x_u, Index m, Number* g_l, Number* g_u)
	{
		DBG_ASSERT(n == n_var_);
		DBG_ASSERT(m == n_con_);

		for (Index i=0; i<n; i++) {
			x_l[i] = x_L_[i];
			x_u[i] = x_U_[i];
		}

		for (Index i=0; i<m; i++) {
			g_l[i] = g_L_[i];
			g_u[i] = g_U_[i];
		}

		return true;
	}

	bool StdInterfaceTMINLP::get_scaling_parameters(
		Number& obj_scaling,
		bool& use_x_scaling, Index n,
		Number* x_scaling,
		bool& use_g_scaling, Index m,
		Number* g_scaling)
	{
		DBG_ASSERT(n==n_var_);
		DBG_ASSERT(m==n_con_);
		obj_scaling = obj_scaling_;
		if (x_scaling_) {
			use_x_scaling = true;
			for (Index i=0; i<n_var_; i++) {
				x_scaling[i] = x_scaling_[i];
			}
		}
		else {
			use_x_scaling = false;
		}
		if (g_scaling_) {
			use_g_scaling = true;
			for (Index i=0; i<n_con_; i++) {
				g_scaling[i] = g_scaling_[i];
			}
		}
		else {
			use_g_scaling = false;
		}
		return true;
	}

	bool StdInterfaceTMINLP::get_variables_types(Index n, VariableType* var_types)
	{
		DBG_ASSERT(n==n_var_);
		memcpy( var_types, var_types_, n*sizeof( VariableType ) );
		return true;
	}

	bool StdInterfaceTMINLP::get_variables_linearity(Index n, TNLP::LinearityType* var_linearity_types)
	{
		DBG_ASSERT(n==n_var_);
		memcpy( var_linearity_types, var_linearity_types_, n*sizeof( TNLP::LinearityType ) );
		return true;
	}

	bool StdInterfaceTMINLP::get_constraints_linearity(Index m, TNLP::LinearityType* constraint_linearity_types)
	{
		DBG_ASSERT(m==n_con_);
		memcpy( constraint_linearity_types, constraint_linearity_types_, m*sizeof( TNLP::LinearityType ) );
		return true;
	}

	bool StdInterfaceTMINLP::get_starting_point(Index n, bool init_x, Number* x, bool init_z, Number* z_L, Number* z_U, Index m, bool init_lambda, Number* lambda)
	{
		bool retval=true;

		DBG_ASSERT(n == n_var_);
		DBG_ASSERT(m == n_con_);

		if (init_x) {
			for (Index i=0; i<n; i++) {
				x[i] = start_x_[i];
			}
		}

		if (init_z) {
			if (start_z_L_==NULL) {
				retval = false;
			}
			else {
				for (Index i=0; i<n; i++) {
					z_L[i] = start_z_L_[i];
				}
			}
			if (start_z_U_==NULL) {
				retval = false;
			}
			else {
				for (Index i=0; i<n; i++) {
					z_U[i] = start_z_U_[i];
				}
			}
		}

		if (init_lambda) {
			if (start_lam_==NULL) {
				retval = false;
			}
			else {
				for (Index i=0; i<m; i++) {
					lambda[i] = start_lam_[i];
				}
			}
		}

		return retval;
	}


	bool StdInterfaceTMINLP::eval_f(Index n, const Number* x, bool new_x,
		Number& obj_value)
	{
		DBG_ASSERT(n==n_var_);

		apply_new_x(new_x, n, x);

		Bool retval = (*eval_f_)(n, non_const_x_, (Bool)new_x,
			&obj_value, user_data_);
		return (retval!=0);
	}

	bool StdInterfaceTMINLP::eval_grad_f(Index n, const Number* x, bool new_x,
		Number* grad_f)
	{
		DBG_ASSERT(n==n_var_);

		apply_new_x(new_x, n, x);

		Bool retval = (*eval_grad_f_)(n, non_const_x_, (Bool)new_x, grad_f,
			user_data_);
		return (retval!=0);
	}

	bool StdInterfaceTMINLP::eval_g(Index n, const Number* x, bool new_x,
		Index m, Number* g)
	{
		DBG_ASSERT(n==n_var_);
		DBG_ASSERT(m==n_con_);

		apply_new_x(new_x, n, x);

		Bool retval = (*eval_g_)(n, non_const_x_, (Bool)new_x, m, g, user_data_);

		return (retval!=0);
	}

	bool StdInterfaceTMINLP::eval_jac_g(Index n, const Number* x, bool new_x,
		Index m, Index nele_jac, Index* iRow,
		Index *jCol, Number* values)
	{
		DBG_ASSERT(n==n_var_);
		DBG_ASSERT(nele_jac==nele_jac_);

		Bool retval=1;

		if ( (iRow && jCol && !values) || (!iRow && !jCol && values) ) {
			apply_new_x(new_x, n, x);
			retval = (*eval_jac_g_)(n, non_const_x_, (Bool)new_x, m, nele_jac,
				iRow, jCol, values, user_data_);
		}
		else {
			DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
		}
		return (retval!=0);
	}

	bool StdInterfaceTMINLP::eval_h(Index n, const Number* x, bool new_x,
		Number obj_factor, Index m,
		const Number* lambda, bool new_lambda,
		Index nele_hess, Index* iRow, Index* jCol,
		Number* values)
	{
		DBG_ASSERT(n==n_var_);
		DBG_ASSERT(m==n_con_);
		DBG_ASSERT(nele_hess==nele_hess_);

		Bool retval=1;

		if ( (iRow && jCol && !values) || (!iRow && !jCol && values) ) {
			apply_new_x(new_x, n, x);
			Number* non_const_lambda = new Number[m];
			if (lambda) {
				for (Index i=0; i<m; i++) {
					non_const_lambda[i] = lambda[i];
				}
			}

			retval = (*eval_h_)(n, non_const_x_, (Bool)new_x, obj_factor, m,
				non_const_lambda, (Bool)new_lambda, nele_hess,
				iRow, jCol, values, user_data_);
			delete [] non_const_lambda;
		}
		else {
			DBG_ASSERT(false && "Invalid combination of iRow, jCol, and values pointers");
		}
		return (retval!=0);
	}

	void StdInterfaceTMINLP::apply_new_x(bool new_x, Index n, const Number* x)
	{
		if (new_x)
		{
			//copy the data to the non_const_x_
			if (!non_const_x_)
			{
				non_const_x_ = new Number[n];
			}

			DBG_ASSERT(x && "x is NULL");
			for (Index i=0; i<n; i++)
			{
				non_const_x_[i] = x[i];
			}
		}
	}

	// bool StdInterfaceTMINLP::intermediate_callback(AlgorithmMode mode,
	// 	Index iter, Number obj_value,
	// 	Number inf_pr, Number inf_du,
	// 	Number mu, Number d_norm,
	// 	Number regularization_size,
	// 	Number alpha_du, Number alpha_pr,
	// 	Index ls_trials,
	// 	const IpoptData* ip_data,
	// 	IpoptCalculatedQuantities* ip_cq)
	// {
	// 	Bool retval = 1;
	// 	if (intermediate_cb_) {
	// 		retval = (*intermediate_cb_)((Index)mode, iter, obj_value, inf_pr, inf_du,
	// 			mu, d_norm, regularization_size, alpha_du,
	// 			alpha_pr, ls_trials, user_data_);
	// 	}
	// 	return (retval!=0);
	// }
}
