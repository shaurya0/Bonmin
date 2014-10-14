#include "BonStdCInterface.h"
#include "BonStdInterfaceTMINLP.hpp"
#include "IpOptionsList.hpp"
#include "IpIpoptApplication.hpp"
#include "BonBonminSetup.hpp"
#include "BonCbc.hpp"

using namespace Bonmin;

struct BonminProblemInfo
{
    Index n;
    Number* x_L;
    Number* x_U;
    Index m;
    Number* g_L;
    Number* g_U;
    Index nele_jac;
    Index nele_hess;
    Index index_style;
    Eval_F_CB eval_f;
    Eval_G_CB eval_g;
    Eval_Grad_F_CB eval_grad_f;
    Eval_Jac_G_CB eval_jac_g;
    Eval_H_CB eval_h;
    VariableTypeC* var_types;
    LinearityTypeC* var_linearity_types;
    LinearityTypeC* constraint_linearity_types;
    // BranchingInfo branch,
    // SosInfo sos,
    Intermediate_CB intermediate_cb;
    BonminSetup bonmin_setup;
    Number obj_scaling;
    Number* x_scaling;
    Number* g_scaling;
};

BonminProblem CreateBonminProblem(
				Index n
				, Number* x_L
				, Number* x_U
				, Index m
				, Number* g_L
				, Number* g_U
				, Index nele_jac
				, Index nele_hess
				, Index index_style
				, Eval_F_CB eval_f
				, Eval_G_CB eval_g
				, Eval_Grad_F_CB eval_grad_f
				, Eval_Jac_G_CB eval_jac_g
				, Eval_H_CB eval_h
				, VariableTypeC* var_types
				, LinearityTypeC* var_linearity_types
				, LinearityTypeC* constraint_linearity_types )
{
	if ( n<1 || m<0 || !x_L || !x_U || (m>0 && (!g_L || !g_U)) ||
	        (m==0 && nele_jac != 0) || (m>0 && nele_jac < 1) || nele_hess < 0 ||
	        !eval_f || !eval_grad_f || (m>0 && (!eval_g || !eval_jac_g)) ||
	        !var_types || !var_linearity_types || (m>0 && !constraint_linearity_types) )
	{
	    return NULL;
	}
    BonminProblem retval = new BonminProblemInfo;

    retval->n = n;
    retval->x_L = new Number[n];
    for (Index i=0; i<n; i++) {
        retval->x_L[i] = x_L[i];
    }
    retval->x_U = new Number[n];
    for (Index i=0; i<n; i++) {
        retval->x_U[i] = x_U[i];
    }

    retval->m = m;
    if (m>0) {
        retval->g_L = new Number[m];
        for (Index i=0; i<m; i++) {
            retval->g_L[i] = g_L[i];
        }
        retval->g_U = new Number[m];
        for (Index i=0; i<m; i++) {
            retval->g_U[i] = g_U[i];
        }
    }
    else {
        retval->g_L = NULL;
        retval->g_U = NULL;
    }

    retval->nele_jac = nele_jac;
    retval->nele_hess = nele_hess;
    retval->index_style = index_style;
    retval->eval_f = eval_f;
    retval->eval_g = eval_g;
    retval->eval_grad_f = eval_grad_f;
    retval->eval_jac_g = eval_jac_g;
    retval->eval_h = eval_h;
    retval->var_types = var_types;
    retval->var_linearity_types = var_linearity_types;
    retval->constraint_linearity_types = constraint_linearity_types;
    retval->intermediate_cb = NULL;

    retval->obj_scaling = 1;
    retval->x_scaling = NULL;
    retval->g_scaling = NULL;

    retval->bonmin_setup.initializeOptionsAndJournalist();

    return retval;
}

void FreeBonminProblem(BonminProblem bonmin_problem)
{
	delete [] bonmin_problem->x_L;
	delete [] bonmin_problem->x_U;
	delete [] bonmin_problem->var_types;
	delete [] bonmin_problem->var_linearity_types;
	delete [] bonmin_problem->constraint_linearity_types;
	if (bonmin_problem->m>0) {
	    delete [] bonmin_problem->g_L;
	    delete [] bonmin_problem->g_U;
	}

	delete [] bonmin_problem->x_scaling;
	delete [] bonmin_problem->g_scaling;

	delete bonmin_problem;
}

Bool AddBonminStrOption(BonminProblem bonmin_problem, char* keyword, char* val)
{
	std::string tag( keyword );
	std::string value( val );
	return (Bool) bonmin_problem->bonmin_setup.options()->SetStringValue( tag, value );
}

Bool AddBonminNumOption(BonminProblem bonmin_problem, char* keyword, Number val)
{
	std::string tag(keyword);
	Ipopt::Number value=val;
	return (Bool) bonmin_problem->bonmin_setup.options()->SetNumericValue( tag, value );
}

Bool AddBonminIntOption(BonminProblem bonmin_problem, char* keyword, Int val)
{
	std::string tag(keyword);
	Ipopt::Index value=val;
	return (Bool) bonmin_problem->bonmin_setup.options()->SetIntegerValue( tag, value );
}

Bool OpenBonminOutputFile(BonminProblem bonmin_problem, char* file_name, Int print_level)
{
    std::string name(file_name);
    Ipopt::EJournalLevel level = Ipopt::EJournalLevel(print_level);
    // return (Bool) bonmin_problem->bonmin_setup.options()->OpenOutputFile(name, level);
    return ( Bool )true;
}

Bool SetBonminProblemScaling(BonminProblem bonmin_problem,
		      Number obj_scaling,
		      Number* x_scaling,
		      Number* g_scaling)
{
	bonmin_problem->obj_scaling = obj_scaling;
	if (x_scaling)
	{
	    if (!bonmin_problem->x_scaling)
	    {
	        bonmin_problem->x_scaling = new Number[bonmin_problem->n];
	    }
	    for (::Index i=0; i<bonmin_problem->n; i++)
	    {
	        bonmin_problem->x_scaling[i] = x_scaling[i];
	    }
	}
	else
	{
	    delete [] bonmin_problem->x_scaling;
	    bonmin_problem->x_scaling = NULL;
	}
	if (g_scaling)
	{
	    if (!bonmin_problem->g_scaling)
	    {
	        bonmin_problem->g_scaling = new Number[bonmin_problem->m];
	    }
	    for (::Index i=0; i<bonmin_problem->m; i++)
	    {
	        bonmin_problem->g_scaling[i] = g_scaling[i];
	    }
	}
	else
	{
	    delete [] bonmin_problem->g_scaling;
	    bonmin_problem->g_scaling = NULL;
	}

	return (Bool)true;
}

Int BonminSolve(
  BonminProblem bonmin_problem
, Number* x
, Number* g
, Number* obj_val
, Number* mult_g
, Number* mult_x_L
, Number* mult_x_U
, UserDataPtr user_data)
{
	using namespace Ipopt;

	// Initialize and process options
	// Ipopt::ApplicationReturnStatus retval = ipopt_problem->app->Initialize();
	// if (retval!=Ipopt::Solve_Succeeded) {
	    // return (::ApplicationReturnStatus) retval;
	// }

	if (!x) {
		return 0;
	    // ipopt_problem->app->Jnlst()->Printf(J_ERROR, J_MAIN,
	                                        // "Error: Array x with starting point information is NULL.");
	    // return (::ApplicationReturnStatus) Ipopt::Invalid_Problem_Definition;
	}

	// Copy the starting point information
	::Number* start_x = new ::Number[bonmin_problem->n];
	for (::Index i=0; i<bonmin_problem->n; i++)
	{
	    start_x[i] = x[i];
	}
	::Number* start_lam = NULL;
	if (mult_g)
	{
	    start_lam = new ::Number[bonmin_problem->m];
	    for (::Index i=0; i<bonmin_problem->m; i++)
	    {
	        start_lam[i] = mult_g[i];
	    }
	}
	::Number* start_z_L = NULL;
	if (mult_x_L)
	{
	    start_z_L = new ::Number[bonmin_problem->n];
	    for (::Index i=0; i<bonmin_problem->n; i++)
	    {
	        start_z_L[i] = mult_x_L[i];
	    }
	}
	::Number* start_z_U = NULL;
	if (mult_x_U)
	{
	    start_z_U = new ::Number[bonmin_problem->n];
	    for (::Index i=0; i<bonmin_problem->n; i++)
	    {
	        start_z_U[i] = mult_x_U[i];
	    }
	}


	SmartPtr<TMINLP> tminlp;
	try
	{

		tminlp = new Bonmin::StdInterfaceTMINLP(bonmin_problem->n, bonmin_problem->x_L,
                                    	bonmin_problem->x_U, bonmin_problem->m,
                                    	bonmin_problem->g_L, bonmin_problem->g_U,
                                    	bonmin_problem->nele_jac,
                                    	bonmin_problem->nele_hess,
                                    	bonmin_problem->index_style,
                                    	start_x, start_lam, start_z_L, start_z_U,
                                    	bonmin_problem->eval_f, bonmin_problem->eval_g,
                                    	bonmin_problem->eval_grad_f,
                                    	bonmin_problem->eval_jac_g,
                                    	bonmin_problem->eval_h,
                                    	bonmin_problem->var_types,
                                    	bonmin_problem->var_linearity_types,
                                    	bonmin_problem->constraint_linearity_types,
                                    	// bonmin_problem->intermediate_cb,
                                    	x, mult_x_L, mult_x_U, g, mult_g,
                                    	obj_val, user_data,
                                    	bonmin_problem->obj_scaling,
                                    	bonmin_problem->x_scaling,
                                    	bonmin_problem->g_scaling);

	bonmin_problem->bonmin_setup.initialize( GetRawPtr( tminlp ) );
	Bonmin::Bab branch_and_bound;
	branch_and_bound( bonmin_problem->bonmin_setup );
	}
	catch(TNLPSolver::UnsolvedError *E)
	{
		//There has been a failure to solve a problem with Ipopt.
		std::cerr<<"Ipopt has failed to solve a problem"<<std::endl;
	}
	catch(OsiTMINLPInterface::SimpleError &E)
	{
		std::cerr<<E.className()<<"::"<<E.methodName()
		<<std::endl
		<<E.message()<<std::endl;
	}
	catch(CoinError &E)
	{
		std::cerr<<E.className()<<"::"<<E.methodName()
		<<std::endl
		<<E.message()<<std::endl;
	}
}