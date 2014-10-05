#ifndef BONSTDINTERFACETMINLP_HPP
#define BONSTDINTERFACETMINLP_HPP

#include "BonTMINLP.hpp"


namespace Bonmin
{
    class StdInterfaceTMINLP : public TMINLP
    {
	public:
    	virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g,
        Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style){ return false; }

        virtual bool get_scaling_parameters(Ipopt::Number& obj_scaling,
                                        bool& use_x_scaling, Ipopt::Index n,
                                        Ipopt::Number* x_scaling,
                                        bool& use_g_scaling, Ipopt::Index m,
                                        Ipopt::Number* g_scaling){ return false; }

        virtual bool get_variables_types(Ipopt::Index n, VariableType* var_types){ return false; }

        virtual bool get_variables_linearity(Ipopt::Index n, Ipopt::TNLP::LinearityType* var_types){ return false; }

		virtual bool get_constraints_linearity(Ipopt::Index m, Ipopt::TNLP::LinearityType* const_types){ return false; }

	    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u,Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u){ return false; }

	    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x,bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U,
	    							    Ipopt::Index m, bool init_lambda,
	    							    Ipopt::Number* lambda){ return false; }

	    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value){ return false; }

	    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f){ return false; }

	    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g){ return false; }

	    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
						        Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow,
						        Ipopt::Index *jCol, Ipopt::Number* values){ return false; }

        virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x,
            				Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda,
            				bool new_lambda, Ipopt::Index nele_hess,
            				Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values){ return false; }

        virtual bool eval_gi(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index i, Ipopt::Number& gi){ return false; }

        virtual bool eval_grad_gi(	Ipopt::Index n, const Ipopt::Number* x, bool new_x,
    			      				Ipopt::Index i, Ipopt::Index& nele_grad_gi, Ipopt::Index* jCol,
    			      				Ipopt::Number* values){ return false; }

        virtual void finalize_solution(TMINLP::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, Ipopt::Number obj_value){}

    };
}


#endif // BONSTDINTERFACETMINLP_HPP
