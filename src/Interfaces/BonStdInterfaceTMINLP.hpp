#ifndef BONSTDINTERFACETMINLP_HPP
#define BONSTDINTERFACETMINLP_HPP

#include "BonTMINLP.hpp"
#include "BonStdCInterface.h"

namespace Bonmin
{
    class StdInterfaceTMINLP : public TMINLP
    {
    public:
        StdInterfaceTMINLP(Index n_var,
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
                         LinearityType* var_linearity_types,
                         ConstraintsLinearityType* constraint_linearity_types,
                         BranchingInfo branch,
                         SosInfo sos,
                         PerturbInfo perturb_info,
                         Intermediate_CB intermediate_cb,
                         Number* x_sol,
                         Number* z_L_sol,
                         Number* z_U_sol,
                         Number* g_sol,
                         Number* lam_sol,
                         Number* obj_sol,
                         UserDataPtr user_data,
                         Number obj_scaling=1,
                         const Number* x_scaling = NULL,
                         const Number* g_scaling = NULL);



        virtual ~StdInterfaceTMINLP();
        /** Default Constructor */
        StdInterfaceTMINLP();

        /** Copy Constructor */
        StdInterfaceTMINLP(const StdInterfaceTMINLP&);

        /** Overloaded Equals Operator */
        void operator=(const StdInterfaceTMINLP&);
        //@}

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

        virtual void finalize_solution(TMINLP::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, Ipopt::Number obj_value){}

        virtual const BranchingInfo * branchingInfo() const
        {
            return &branch_;
        }

        virtual const SosInfo * sosConstraints() const
        {
            return &sos_;
        }

        virtual const PerturbInfo* perturbInfo() const
        {
            return &perturb_info_;
        }


    private:
        /** Journlist */
        SmartPtr<const Journalist> jnlst_;

        /** @name Information about the problem */
        //@{
        /** Number of variables */
        const Index n_var_;
        /** Number of constraints */
        const Index n_con_;
        /** Pointer to Number array containing lower bounds for variables */
        const Number* x_L_;
        /** Pointer to Number array containing upper bounds for variables */
        const Number* x_U_;
        /** Pointer to Number array containing lower bounds for constraints */
        const Number* g_L_;
        /** Pointer to Number array containing upper bounds for constraints */
        const Number* g_U_;
        /** Number of non-zero elements in the constraint Jacobian */
        const Index nele_jac_;
        /** Number of non-zero elements in the Hessian */
        const Index nele_hess_;
        /** Starting value of the iRow and jCol parameters for matrices */
        const Index index_style_;
        /** Pointer to Number array containing starting point for variables */
        const Number* start_x_;
        /** Poitner to Number array containing starting values for
        *  constraint multipliers */
        const Number* start_lam_;
        /** Pointer to Number array containing starting values for lower
        *  bound multipliers */
        const Number* start_z_L_;
        /** Pointer to Number array containing starting values for upper
        *  bound multipliers */
        const Number* start_z_U_;
        /** Pointer to callback function evaluating value of objective function */
        Eval_F_CB eval_f_;
        /**  Pointer to callback function evaluating value of constraints */
        Eval_G_CB eval_g_;
        /** Pointer to callback function evaluating gradient of objective
        *  function */
        Eval_Grad_F_CB eval_grad_f_;
        /** Pointer to callback function evaluating Jacobian of constraints */
        Eval_Jac_G_CB eval_jac_g_;

        /** Pointer to callback function evaluating Hessian of Lagrangian */
        Eval_H_CB eval_h_;

        VariableType* var_types_;
        LinearityType* var_linearity_types_;
        LinearityType* constraint_linearity_types_;

        /** Pointer to intermediate callback function giving control to user */
        Intermediate_CB intermediate_cb_;

        /** Storage of branching priorities information.*/
        BranchingInfo branch_;

        /** Storage of sos constraints */
        SosInfo sos_;

        /** Storage for perturbation radii */
        PerturbInfo perturb_info_;

        /** Pointer to user data */
        UserDataPtr user_data_;
        /** Objective scaling factor */
        Number obj_scaling_;
        /** Scaling factors for variables (if not NULL) */
        const Number* x_scaling_;
        /** Scaling factors for constraints (if not NULL) */
        const Number* g_scaling_;
        //@}


        /** A non-const copy of x - this is kept up-to-date in apply_new_x */
        Number* non_const_x_;

        /** Pointers to the user provided vectors for solution */
        Number* x_sol_;
        Number* z_L_sol_;
        Number* z_U_sol_;
        Number* g_sol_;
        Number* lambda_sol_;
        Number* obj_sol_;
    };
}


#endif // BONSTDINTERFACETMINLP_HPP
