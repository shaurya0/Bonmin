#ifndef BONHEURISTICOFP_HPP
#define BONHEURISTICOFP_HPP

#include "BonBonminSetup.hpp"
#include "CbcHeuristic.hpp"
#include <algorithm>

namespace Bonmin
{
  class  HeuristicOFP : public CbcHeuristic
  {
  public:
    /// Default constructor
    HeuristicOFP();

    /// Constructor with setup
    HeuristicOFP(BonminSetup * setup);

    /// Copy constructor
    HeuristicOFP(const HeuristicOFP &copy);

    /// Destructor
    ~HeuristicOFP() {}

    /// Assignment operator
    HeuristicOFP & operator=(const HeuristicOFP & rhs);

    /** Virtual constructor.*/
    virtual CbcHeuristic * clone() const{
      return new HeuristicOFP(*this);
    }

    /// Resets stuff if model changes
    virtual void resetModel(CbcModel * model){
      setModel(model);
    }

    /** Change setup used for heuristic.*/
    void setSetup(BonminSetup * setup){
      setup_ = setup;
      Initialize(setup_->options());
    }

    // Compute l1_distance between continuous variable and rounded point
    double l1_distance(const double *x, const double* x_tilde, int numberIntegerColumns, std::vector<int> integerColumns);

    // Compute number of fractional variables
    int num_fractional_vars(const double *newSolution ,int numberIntegerColumns, std::vector<int> integerColumns, double integerTolerance);


    // Performs heuristic
    virtual int solution(double &solutionValue, double *betterSolution);

    /// Performs heuristic with add cust
    virtual int solution(double &solutionValue, double *betterSolution, OsiCuts & cs)
    {
      return solution(solutionValue, betterSolution);
    }

    /** Register the options for this heuristic */
    static void registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions);

    /** Initiaize using passed options.*/
    void Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options);

  private:
    /** Setup to use for local searches (will make copies).*/
    BonminSetup * setup_;

    /** Norm of the objective function - either 1 or 2 */
    int objective_norm_;

    // objective function user weight (u_2), value is read through options
    double obj_user_weight_;

    // distance function user weight (u_1), value is read through options
    double dist_user_weight_;

    // geometric reduction criterion for cycle handling to start
    double delta_alpha_;

  };

  class RoundingOFP
  {
  public:
    /// Default constructor
    RoundingOFP(TMINLP2TNLP* minlp, const double primalTolerance_, const double integerTolerance_);

    /// Destructor
    ~RoundingOFP();

    /// Rounds the solution
    void round(double* solution);

    // Computes ||max{0,g(x)}||, norm of constraint violation for all constraints
    double constraint_violation(const double *x);

    // Computes ||max{0,g(x)}||, i.e. norm of constraint violation for nonlinear constraints
    double nonlinear_constraint_violation(const double *x);

    // For a given variable index, checks if it belongs to a SOS constraint, inlined
    bool is_sos_var(int idx);

    // Returns the constraint row corresponding to a variable that is SOS1
    int sos_row(int idx);

    // Checks if a given point is SOS1 feasible
    bool is_sos_feasible(const double *x);

    // Returns true if the MINLP has atleast one SOS1 constraint
    bool has_sos_constraints() const;

    // Perform simple rounding on x[idx], inlined
    double simple_round(const double x_idx);

    // Constraint directed rounding on an index idx
    double constraint_round(const double *x, int idx);

    // Given indices that are being flipped as a result of stalling, this function makes sure
    // that the resulting flipped variables are SOS1 feasible. This includes modification of
    // the values
    void stall_sos_feasibility( double *rounded_solution, std::vector<int> &sos_row_inds );

    // Takes a row with a SOS1 constraint and mutates x to make that row feasible
    void sos_round_row_constraint(int iRow, double *rounded_solution);


  private:
    /// gutsOfConstructor
    void gutsOfConstructor();

    /// Pointer to problem
    TMINLP2TNLP* minlp_;

    /// Number of rows
    int numberRows_;

    /// Number of columns
    int numberColumns_;

    // Upper and lower bounds of the original constraints
    const double *x_l, *x_u, *g_l, *g_u;

    double *g_temp;

    // Store variable types
    const Bonmin::TMINLP::VariableType* variableType;

    // Maps from row number to variable indices which appear in SOS1 constraint in that row
    std::map< int,std::vector<int> > sos_constraints;

    // Key: variable index. Value : row index of SOS1 constraint in which the variable appears
    std::map<int, int> sos_var_row;

    // Sorted indices of variables which appear in SOS1 constraints
    std::vector<int> sos_columns;

    // Sorted indices of variables which are integer but do not appear in SOS1 constraints
    std::vector<int> int_not_sos;


    const double integerTolerance, primalTolerance;

  };

}



#endif // BONHEURISTICOFP_HPP
