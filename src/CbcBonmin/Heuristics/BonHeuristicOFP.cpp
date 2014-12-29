#include "BonHeuristicOFP.hpp"
#include "CoinHelperFunctions.hpp"
#include "CbcModel.hpp"

#include "OsiAuxInfo.hpp"

#include "CoinTime.hpp"

#include <fstream>
#include <iomanip>

using namespace std;

//#define DEBUG_BON_HEURISTIC_FPUMP

namespace Bonmin
{
class score_sorter {
public:
    //! Constructor
    score_sorter(const vector<double>& score):
        score_(score) {}

    bool operator() (const int x, const int y) const {
        return score_[x]>score_[y];
    }

private:
    const vector<double>& score_;
};


HeuristicOFP::HeuristicOFP()
    :
      CbcHeuristic(),
      setup_(NULL),
      objective_norm_(1),
      obj_user_weight_(1.0),
      dist_user_weight_(1.0),
      delta_alpha_(0.005)
{}

HeuristicOFP::HeuristicOFP(BonminSetup * setup)
    :
      CbcHeuristic(),
      setup_(setup),
      objective_norm_(1)
{
    Initialize(setup->options());
}

HeuristicOFP::HeuristicOFP(const HeuristicOFP &copy)
    :
      CbcHeuristic(copy),
      setup_(copy.setup_),
      objective_norm_(copy.objective_norm_),
      obj_user_weight_(copy.obj_user_weight_),
      dist_user_weight_(copy.dist_user_weight_),
      delta_alpha_(copy.delta_alpha_)
{
}

HeuristicOFP &
HeuristicOFP::operator=(const HeuristicOFP & rhs)
{
    if(this != &rhs) {
        CbcHeuristic::operator=(rhs);
        setup_ = rhs.setup_;
        objective_norm_ = rhs.objective_norm_;
        obj_user_weight_ = rhs.obj_user_weight_;
        dist_user_weight_ = rhs.dist_user_weight_;
        delta_alpha_ = rhs.delta_alpha_;
    }
    return *this;
}




double HeuristicOFP::l1_distance(const double *x, const double* rounded_solution, int numberIntegerColumns, std::vector<int> integerColumns){
    double l1_d = 0;
    for (int i=0; i < numberIntegerColumns; i++)
    {
        int j = integerColumns[i];
        l1_d += abs(x[j]-rounded_solution[j]);
    }
    return l1_d;
}

int HeuristicOFP::num_fractional_vars(const double *newSolution ,int numberIntegerColumns, std::vector<int> integerColumns, double integerTolerance){
    int numberFractionalVariables = 0;
    for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++) {
        int iColumn = integerColumns[iIntCol];
        double value=newSolution[iColumn];
        if (fabs(floor(value+0.5)-value)>integerTolerance)
            numberFractionalVariables++;
    }
    return numberFractionalVariables;
}



int
HeuristicOFP::solution(double &solutionValue, double *betterSolution)
{
    //measure time
    clock_t start = clock();

    // need to call setSetup twice for some reason
    setSetup(setup_);
    if(model_->getNodeCount() || model_->getCurrentPassNumber() > 1) return 0;
    bool integerSolutionAlreadyExists = false;

    if(model_->getSolutionCount()) {
        //      bestSolutionValue = model_->getObjValue();
        integerSolutionAlreadyExists = true;
        assert(solutionValue < 1.0e50);
    }
    const int maxNumberIterations = 200;
    int returnCode = 0; // 0 means it didn't find a feasible solution
    OsiTMINLPInterface * nlp = NULL;
    if(setup_->getAlgorithm() == B_BB)
        nlp = dynamic_cast<OsiTMINLPInterface *>(model_->solver()->clone());
    else
        nlp = dynamic_cast<OsiTMINLPInterface *>(setup_->nonlinearSolver()->clone());

    TMINLP2TNLP* minlp = nlp->problem();
    //set tolerances
    double integerTolerance = model_->getDblParam(CbcModel::CbcIntegerTolerance);
    double primalTolerance;
#if 0
    OsiSolverInterface * solver = model_->solver();
    solver->getDblParam(OsiPrimalTolerance,primalTolerance);
#endif
    primalTolerance=1.0e-6;
    //used as a termination criterion
    double obj_nlp_ofp = 1e20; //just some large value

    int numberColumns;
    int numberRows;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp->get_nlp_info(numberColumns, numberRows, nnz_jac_g,
                        nnz_h_lag, index_style);

    const Bonmin::TMINLP::VariableType* variableType = minlp->var_types();
    const double* x_sol = minlp->x_sol();
    const double* x_l = minlp->x_l();
    const double* x_u = minlp->x_u();
    double* x_ofp = new double[numberColumns];
    memcpy(x_ofp,x_sol,numberColumns*sizeof(double));
    const double nlp_utopia = minlp->obj_value();

    // exit if the current NLP solution is infeasible
    // infeasibility is determined by Ipopt
    if(minlp->optimization_status() != Ipopt::SUCCESS){
        delete nlp;
        return returnCode;
    }


    // create a set with the indices of the fractional variables

    vector<int> integerColumns; // stores the integer variables
    int numberFractionalVariables = 0;
    for (int iColumn=0;iColumn<numberColumns;iColumn++) {
        if (variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS)
        {
            integerColumns.push_back(iColumn);
            double value=x_ofp[iColumn];
            if (fabs(floor(value+0.5)-value)>integerTolerance)
            {
                numberFractionalVariables++;
            }
        }
    }
    int numberIntegerColumns = (int) integerColumns.size();
    // cout << "numberIntegerColumns == " << numberIntegerColumns << endl;
    // for_each(integerColumns.begin(), integerColumns.end(),print_);
    // cout << endl;

    // create space to store old solutions in order to prevent cycling
    const int numberOldSolutionsStored = 4;

    //store corresponding alpha values
    double *alpha_t = new double [numberOldSolutionsStored];
    double ** oldSolution = new double * [numberOldSolutionsStored];
    for (int j=0;j<numberOldSolutionsStored;j++){
        oldSolution[j]= new double[numberIntegerColumns];
        alpha_t[j] = COIN_DBL_MAX;
        for (int i=0;i<numberIntegerColumns;i++)
            oldSolution[j][i]=-COIN_DBL_MAX;
    }

    RoundingOFP roundObj(minlp, integerTolerance, primalTolerance);

    double* x_tilde = new double[numberIntegerColumns];
    int* indexes_x_tilde = new int[numberIntegerColumns];

    double* rounded_solution = new double[numberColumns];
    double *newSolution = new double[numberColumns];
    int ofp_iteration = 0;
    int stall_counter = 0;
    int cycle_counter = 0;

    //parameters for objective FP
    double alpha = 1.0;
    double phi = 0.9;

    double objectiveScalingFactor, distanceScalingFactor; //objectiveScalingFactor and distanceScalingFactor: normalizations for the objectives
    memcpy(rounded_solution, x_ofp, numberColumns*sizeof(double));
    roundObj.round(rounded_solution);

    for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++)
    {
        int iColumn = integerColumns[iIntCol];
        x_tilde[iIntCol] = rounded_solution[iColumn];
        indexes_x_tilde[iIntCol] = iColumn;
    }

    double temp = nlp->solveFeasibilityProblem(numberIntegerColumns,
                                               x_tilde,indexes_x_tilde,objective_norm_);

    double fp_utopia, fp_nadir, nlp_nadir;
    fp_nadir = l1_distance( x_ofp, rounded_solution, numberIntegerColumns, integerColumns );
    fp_utopia = l1_distance( x_sol, rounded_solution, numberIntegerColumns, integerColumns );
    minlp->eval_f( numberColumns, x_sol, true, nlp_nadir );

    // some small value
    static const double epsilon = 1e-12;

    //  check for divide by zero possibility
    double temp_scaling_value = fp_nadir - fp_utopia;
    if( std::abs(temp_scaling_value) <= epsilon )
    {
        temp_scaling_value = 1;
    }

    distanceScalingFactor = dist_user_weight_/temp_scaling_value;
    temp_scaling_value = nlp_nadir - nlp_utopia;
    if( std::abs(temp_scaling_value) <= epsilon )
    {
        temp_scaling_value = 1;
    }
    objectiveScalingFactor = obj_user_weight_/temp_scaling_value;



    bool sufficient_alpha_decrease;
    while( numberFractionalVariables ) { //ofp argument: obj_nlp_ofp > toleranceObjectiveFP || numberFractionalVariables
        ofp_iteration++;

        if( ofp_iteration > maxNumberIterations )
            break;
        alpha = phi*alpha;
        sufficient_alpha_decrease = false;

        obj_nlp_ofp = nlp->solveObjectiveFeasibilityProblem(numberIntegerColumns
                                                    		, x_tilde
                                                    		, indexes_x_tilde
                                                    		, alpha
                                                    		, objectiveScalingFactor
                                                    		, distanceScalingFactor);


        alpha_t[0] = alpha;

        memcpy( x_ofp,x_sol,numberColumns*sizeof(double) );
        memcpy( rounded_solution,x_ofp,numberColumns*sizeof(double) );
        roundObj.round( rounded_solution );

        // double ofp_obj;
        // minlp->eval_f(numberColumns, x_ofp, true, ofp_obj);
        // double l1_dist = l1_distance(x_ofp, rounded_solution, numberIntegerColumns, integerColumns);
        // std::cout << ofp_obj << "\t" << l1_dist << std::endl;

        numberFractionalVariables = num_fractional_vars(x_ofp,numberIntegerColumns, integerColumns, integerTolerance);
        if (numberFractionalVariables == 0)
            break;

        bool flip = true;
        for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++)
        {
            int iColumn = integerColumns[iIntCol];
            double value=rounded_solution[iColumn];
            x_tilde[iIntCol]=value;
            if(flip && fabs(x_tilde[iIntCol]-oldSolution[0][iIntCol])>integerTolerance)
                flip = false;
        }

        if(flip)
        {
            stall_counter++;
            vector<int> sortedIntegerColumns(numberIntegerColumns);
            vector<double> score(numberIntegerColumns);
            for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++)
            {
                int iColumn = integerColumns[iIntCol];
                sortedIntegerColumns[iIntCol] = iIntCol;
                double value=x_ofp[iColumn];
                score[iIntCol] = fabs(value-oldSolution[0][iIntCol]);
            }

            // Sort according to fractionality
            sort(sortedIntegerColumns.begin(),sortedIntegerColumns.end(),
                 score_sorter(score));

            int maxNumberToMove = int(10*CoinDrand48()) + 1;
            int numberMoved = 0;
            std::vector<int> indices_flipped;
            for(int i=0; i<numberIntegerColumns; i++)
            {
                if( numberMoved >= maxNumberToMove)
                    break;

                int iIntCol = sortedIntegerColumns[i];
                if(score[iIntCol] > 0.00)
                {
                    int iColumn = integerColumns[iIntCol];
                    double value=x_ofp[iColumn];
                    indices_flipped.push_back( iColumn );
                    if(value-oldSolution[0][iIntCol]>0.0)
                        value = oldSolution[0][iIntCol]+1.0;
                    else
                        value = oldSolution[0][iIntCol]-1.0;
                    if(value < x_l[iColumn]-primalTolerance)
                        value++;
                    else if(value > x_u[iColumn]+primalTolerance)
                        value--;
                    assert(fabs(floor(value+0.5)-value)<=integerTolerance);
                    x_tilde[iIntCol] = value;
                    rounded_solution[iColumn] = value;
                    numberMoved++;
                }
                else
                    break;
            }
            if ( roundObj.has_sos_constraints() )
            {
                if ( numberMoved > 0 && !roundObj.is_sos_feasible( rounded_solution ) )
                {
                    // Store the row index, and the indices of the
                    // variables that appear in that row and are also
                    // to be flipped
                    std::vector<int> sos_row_inds;
                    for (int j = 0; j < indices_flipped.size(); j++)
                    {
                        int iColumn = indices_flipped[j];
                        if( roundObj.is_sos_var( iColumn ) )
                        {
                            int iRow = roundObj.sos_row( iColumn );
                            sos_row_inds.push_back( iRow );
                        }
                    }
                    roundObj.stall_sos_feasibility(rounded_solution, sos_row_inds);
                    for(int iIntCol = 0; iIntCol <numberIntegerColumns; iIntCol++)
                    {
                        int iColumn = integerColumns[iIntCol];
                        x_tilde[iIntCol] = rounded_solution[iColumn];
                    }
                }
            }


            // Cycle detection
            bool matched = false;
            for (int k = numberOldSolutionsStored-1; k > 0; k--)
            {
                if (alpha_t[k] - alpha <= delta_alpha_ && ofp_iteration > numberOldSolutionsStored )
                    sufficient_alpha_decrease = true;
                else
                    continue;
                double * b = oldSolution[k];
                matched = true;
                for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++)
                {
                    if (fabs(x_tilde[iIntCol]-b[iIntCol])>integerTolerance)
                    {
                        matched=false;
                        break;
                    }
                }

                if (matched)
                    break;
            }
            if (matched && sufficient_alpha_decrease)
            {
                cycle_counter++;
                // perturbation
                for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++)
                {
                    int iColumn = integerColumns[iIntCol];
                    double value = x_ofp[iColumn];
                    double random = max(0.0,CoinDrand48()-0.3);
                    double difference = fabs(value-oldSolution[0][iIntCol]);
                    if(difference+random>0.5)
                    {
                        if(value-oldSolution[0][iIntCol]>0.0)
                            value = oldSolution[0][iIntCol]+1.0;
                        else
                            value = oldSolution[0][iIntCol]-1.0;
                        // make sure that the new value is within bounds
                        if(value < x_l[iColumn]-primalTolerance)
                            value++;
                        else if(value > x_u[iColumn]+primalTolerance)
                            value--;
                        assert(fabs(floor(value+0.5)-value)<=integerTolerance);
                    }
                    else
                    {
                        value = oldSolution[0][iIntCol];
                    }
                    x_tilde[iIntCol]=value;
                }
            }
        }
        // store the new solution and remove the oldest one
        for (int j=numberOldSolutionsStored-1;j>0;j--) {
            alpha_t[j] = alpha_t[j-1];
            for (int i = 0; i < numberIntegerColumns; i++){
                oldSolution[j][i]=oldSolution[j-1][i];

            }
        }
        for (int j = 0; j < numberIntegerColumns; j++){
            oldSolution[0][j] = x_tilde[j];
        }

    }//end main while loop

    // For plotting purposes
    // double ofp_dist = l1_distance(x_ofp, x_tilde, numberIntegerColumns,integerColumns);
    // double ofp_obj_;
    // minlp->eval_f(numberColumns, x_ofp, true, ofp_obj_);
    // cout << ofp_dist << "\t" << ofp_obj_ << endl;

    double ofp_obj;
    minlp->eval_f(numberColumns, x_ofp, true, ofp_obj);

    // fix the integer variables and solve the NLP
    for(int iIntCol=0; iIntCol<numberIntegerColumns; iIntCol++)
    {
        int iColumn = integerColumns[iIntCol];
        double value=floor(x_ofp[iColumn]+0.5);
        minlp->SetVariableUpperBound(iColumn, floor(value));
        minlp->SetVariableLowerBound(iColumn, ceil(value));
    }
    clock_t end = clock();
    float seconds = (float)(end-start) / CLOCKS_PER_SEC;
    std::cout << "------------ OFP results ------------" << std::endl;
    std::cout << "Total time == " << seconds << std::endl;
    std::cout << "Iterations == " << ofp_iteration << std::endl;
//    std::cout << "Stall count == " << stall_counter << std::endl;
//    std::cout << "Cycle count == " << cycle_counter << std::endl;
    nlp->initialSolve();

    bool feasible = true;
    if( minlp->optimization_status() != Ipopt::SUCCESS || ofp_iteration == maxNumberIterations )
    {
        feasible = false;
    }

    memcpy(newSolution,x_sol,numberColumns*sizeof(double));
    if( feasible )
    {
        double newSolutionValue;
        minlp->eval_f(numberColumns, newSolution, true, newSolutionValue);
        if( newSolutionValue < solutionValue )
        {
            memcpy( betterSolution, newSolution, numberColumns*sizeof(double) );
            solutionValue = newSolutionValue;
            returnCode = 1;
        }
    }

#if 0
    delete [] indexRow;
    delete [] indexCol;
    delete [] row;
    delete [] columnStart;
    delete [] columnLength;
#endif
    // delete [] rounded_solution;
    delete [] x_ofp;
    delete [] alpha_t;
    for (int j=0;j<numberOldSolutionsStored;j++)
        delete [] oldSolution[j];
    delete [] oldSolution;
    delete [] x_tilde;
    delete [] indexes_x_tilde;
    delete nlp;
    return returnCode;
}

void
HeuristicOFP::registerOptions(Ipopt::SmartPtr<Bonmin::RegisteredOptions> roptions){
    roptions->SetRegisteringCategory("MINLP Heuristics", RegisteredOptions::BonminCategory);

    roptions->AddStringOption2("heuristic_objective_feasibility_pump", "whether the heuristic objective feasibility pump should be used",
                               "no", "no", "don't use it", "yes", "use it", "");
    roptions->AddLowerBoundedNumberOption("ofp_objective_weight","user defined weight for the original objective function",
    							0.00001,true,1,"");
    roptions->AddLowerBoundedNumberOption("ofp_distance_weight","user defined weight for the distance function",
    							0.00001,true,1,"");
    roptions->AddBoundedIntegerOption("ofp_cycle_length","minimum number of iterations needed for cycle handling",
    								1,200, 30, "");


    roptions->setOptionExtraInfo("heuristic_objective_feasibility_pump", 63);
    roptions->setOptionExtraInfo("ofp_objective_weight", 63);
    roptions->setOptionExtraInfo("ofp_distance_weight", 63);
    roptions->setOptionExtraInfo("ofp_cycle_length", 63);
}

void
HeuristicOFP::Initialize(Ipopt::SmartPtr<Ipopt::OptionsList> options){
    options->GetNumericValue("ofp_objective_weight", obj_user_weight_, "bonmin.");
    options->GetNumericValue("ofp_distance_weight", dist_user_weight_, "bonmin.");

    int cycle_length;
    options->GetIntegerValue("ofp_cycle_length", cycle_length, "bonmin.");
    delta_alpha_ = 0.1*std::pow(0.9,cycle_length-1);
}

RoundingOFP::RoundingOFP(TMINLP2TNLP* minlp,
                             const double primalTolerance_,
                             const double integerTolerance_)
    :
      minlp_(minlp),
      primalTolerance(primalTolerance_),
      integerTolerance(integerTolerance_)
{
    gutsOfConstructor();
}

RoundingOFP::~RoundingOFP()
{
    delete [] g_temp;
    // delete [] col_and_jac_g_;
    //delete x_l, x_u, g_l, g_u;
}

void
RoundingOFP::gutsOfConstructor()
{
    std::vector<std::pair<int, int> > *col_and_jac_g_;
    int nnz_jac_g;
    int nnz_h_lag;
    Ipopt::TNLP::IndexStyleEnum index_style;
    minlp_->get_nlp_info(numberColumns_, numberRows_, nnz_jac_g,
                         nnz_h_lag, index_style);
    //
    const double* x_sol = minlp_->x_sol();

    // double *xsol = new double[numberColumns_];
    // std::copy(x_sol,x_sol+numberColumns_, xsol);

    x_l = minlp_->x_l();
    x_u = minlp_->x_u();
    g_l = minlp_->g_l();
    g_u = minlp_->g_u();

    variableType = minlp_->var_types();

    // Get the indicies of the jacobian
    // This is also a way of knowing which variables are
    // used in each row
    int* indexRow = new int[nnz_jac_g];
    int* indexCol = new int[nnz_jac_g];
    minlp_->eval_jac_g(numberColumns_, x_sol, false,
                       numberRows_, nnz_jac_g,
                       indexRow, indexCol, 0);

    g_temp = new double [numberRows_];

    // get the jacobian for the solution with zeros
    double* jac_g = new double [nnz_jac_g];
    double* zero_sol = new double [numberColumns_];
    minlp_->get_starting_point(numberColumns_, 1, zero_sol, 0, NULL, NULL, numberRows_, 0, NULL);
    //memset(zero_sol, 0, numberColumns_ * sizeof(double));
    minlp_->eval_jac_g(numberColumns_, zero_sol, true,
                       numberRows_, nnz_jac_g,
                       0, 0, jac_g);

    col_and_jac_g_ = new vector<pair<int, int> >[numberRows_];

    int indexCorrection = (index_style == Ipopt::TNLP::C_STYLE) ? 0 : 1;
    for(int i=0; i<nnz_jac_g; i++)
    {
        int thisIndexRow = indexRow[i]-indexCorrection;
        int thisIndexCol = indexCol[i]-indexCorrection;
        pair<int, int> value(thisIndexCol, static_cast<int>(jac_g[i]));
        col_and_jac_g_[thisIndexRow].push_back(value);
    }


    // Temporary set to store variables which appear in sos constraints
    std::set<int> sos_columns_;
    for(int iRow=0; iRow<numberRows_; iRow++)
    {
        vector<int> cols;
        bool sosConstraint = false;
        if(g_l[iRow] == 1.0 && g_u[iRow] == 1.0)
        {
            sosConstraint = true;
            vector<pair<int, int> > jac_g = col_and_jac_g_[iRow];
            for (unsigned int j=0; j<jac_g.size(); j++)
            {
                int iColumn = jac_g[j].first;
                if (x_sol[iColumn]>=1.0-integerTolerance || jac_g[j].second != 1.0 || variableType[iColumn] == Bonmin::TMINLP::CONTINUOUS) {
                    sosConstraint = false;
                    cols.clear();
                    break;
                }
                else
                {
                    sos_columns_.insert(iColumn); // Placeholder, to be inserted into vector
                    cols.push_back(iColumn);
                }
            }
        }
        if (sosConstraint)
        {
            sos_constraints[iRow] = cols;
            for (int i = 0; i < cols.size(); ++i)
            {
                sos_var_row[cols[i]] = iRow;
            }
        }

    }

    for (int iColumn=0; iColumn < numberColumns_ ; iColumn++)
    {
        std::set<int>::iterator it = sos_columns_.find(iColumn);
        if ( variableType[iColumn] != Bonmin::TMINLP::CONTINUOUS && it == sos_columns_.end() ) // Integer var and not in a SOS1 constraint
        {
            int_not_sos.push_back(iColumn);
        }
    }

    // Insert the SOS1 variable indices from the set into the member vector
    std::copy(sos_columns_.begin(), sos_columns_.end(), std::back_inserter(sos_columns));

    delete [] indexRow;
    delete [] indexCol;
    delete [] jac_g;
    delete [] zero_sol;
}

void
RoundingOFP::round(double* solution)
{
    // Performs SOS rounding for the appropriate variables
    if ( !sos_constraints.empty() )
    {
        for (std::map< int, std::vector<int> >::iterator it = sos_constraints.begin(); it != sos_constraints.end(); ++it)
        {
            double weighted_sum = 0.0;
            int counter = 1;
            for (std::vector<int>::iterator vit = it->second.begin(); vit != it->second.end(); ++vit)
            {
                weighted_sum += counter * solution[*vit];
                counter++;
            }
            double fl = floor(weighted_sum + 0.5); // simple rounding or quadratic rounding
            int indexColumnSelected = static_cast<int>(fl) - 1;
            if(indexColumnSelected < 0)
            {
                continue;
            }
            assert(indexColumnSelected < it->second.size());
            int j = 0;
            for (std::vector<int>::iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
                int iColumn = *vit;
                if (j == indexColumnSelected)
                    solution[iColumn] = 1.0;
                else
                    solution[iColumn] = 0.0;
                j++;
            }
        }
    }

    // Perform rounding for the remaining integer variables
    if (!int_not_sos.empty())
    {
        std::vector<int>::iterator it;
        for (it = int_not_sos.begin(); it != int_not_sos.end(); it++)
        {
            int iColumn = *it;
            double value=solution[iColumn];
            // cout << "value before round == " << value << endl;
            if (fabs(floor(value+0.5)-value)>integerTolerance) // Make sure that the current value is not integer
            {
                // Simple round
                value = simple_round(value);

                // Constraint rounding
                // value = constraint_round(solution, iColumn);

                // make sure that the new value is within bounds
                if(value < x_l[iColumn]-primalTolerance)
                    value++;
                else if(value > x_u[iColumn]+primalTolerance)
                    value--;
                solution[iColumn] = value;
            }
        }
    }

}

double RoundingOFP::constraint_violation(const double *x)
{
    minlp_->eval_g(numberColumns_, x, false, numberRows_, g_temp);
    double violation = 0.0;
    for (int iRow = 0; iRow < numberRows_; iRow++)
    {
        double upper_violation = 0.0;
        double lower_violation = 0.0;
        if(g_u[iRow] < 1e10)
        {
            upper_violation = max(0.0,g_temp[iRow] - g_u[iRow]);
            if (upper_violation < primalTolerance)
            {
                upper_violation = 0.0;
            }
        }
        if(g_l[iRow] > -1e10)
        {
            lower_violation = max(0.0,g_l[iRow] - g_temp[iRow]);
            if (lower_violation < primalTolerance)
            {
                lower_violation = 0.0;
            }
        }
        violation += max(upper_violation,lower_violation);
    }

    if (violation < primalTolerance)
    {
        violation = 0.0;
    }
    return max(0.0,violation);
}

bool RoundingOFP::is_sos_feasible(const double *x)
{
    if (!sos_constraints.empty()){
        for (std::map< int, std::vector<int> >::iterator it = sos_constraints.begin(); it != sos_constraints.end(); ++it){
            double sum = 0.0;
            for (std::vector<int>::iterator vit = it->second.begin(); vit != it->second.end(); ++vit){
                // Check that each index is integer
                double value = x[*vit];
                if ( fabs(floor(value+0.5)-value)>integerTolerance )
                    return false;
                else
                {
                    sum += value;
                }
            }
            // Check that feasibility is within tolerable bounds, else return false
            if ( sum > 1.0 + primalTolerance || sum < 1.0 - primalTolerance)
                return false;
        } // End map loop
    }
    return true;
}


// Too short function to be worth anything, will insert it directly when needed instead
inline bool RoundingOFP::is_sos_var(int idx)
{
    // Binary search returns true if the index is in the sos_columns vector
    return std::binary_search(sos_columns.begin(), sos_columns.end(), idx );
}

inline double RoundingOFP::simple_round(const double x_idx)
{
    return floor(x_idx + 0.5);
}

double RoundingOFP::constraint_round(const double *x, int idx)
{
    double c0,c1;
    double value = x[idx];
    if(x[idx] <= integerTolerance && x[idx] >= -integerTolerance)
        return value;
    else if(x[idx] <= 1.0 + integerTolerance && x[idx] >= 1.0 - integerTolerance)
        return value;
    else
    {
        double * x_copy = new double[numberColumns_];
        std::copy(x,x+numberColumns_, x_copy);
        x_copy[idx] = 0.0;
        c0 = constraint_violation(x_copy);
        x_copy[idx] = 1.0;
        c1 = constraint_violation(x_copy);
        if ( c0 == c1 )
        {
            value = simple_round( x_copy[idx] );
        }
        else if (c0 < c1)
            value = 0.0;
        else
            value = 1.0;
        return value;
    }
}

inline int RoundingOFP::sos_row(int idx)
{
    return sos_var_row[idx];
}

void RoundingOFP::stall_sos_feasibility(double *rounded_solution, std::vector<int> &sos_row_inds)
{
    for (std::vector<int>::iterator it = sos_row_inds.begin(); it != sos_row_inds.end(); ++it)
    {
        int iRow = *it;
        std::vector<int> sos_col_inds = sos_constraints[iRow];
        int sum = 0;
        for (std::vector<int>::iterator vit = sos_col_inds.begin(); vit != sos_col_inds.end(); vit++)
        {
            int iColumn = *vit;
            sum += rounded_solution[iColumn];
        }

        if( sum != 1) //SOS1 violation
            sos_round_row_constraint( iRow, rounded_solution );
    }
}

void RoundingOFP::sos_round_row_constraint(int iRow, double *rounded_solution)
{
    std::vector<int> sos_inds = sos_constraints[iRow];
    int inds_size = sos_inds.size();

    double min_violation = COIN_DBL_MAX;
    std::vector<double> best_sol( inds_size );
    for( std::vector<int>::iterator it = sos_inds.begin(); it != sos_inds.end(); it++ )
    {
        double violation;
        int iColumn = *it;
        rounded_solution[iColumn] = 1.0;
        for( std::vector<int>::iterator vit = sos_inds.begin(); vit != sos_inds.end(); vit++ )
        {
            int jColumn = *vit;
            if( jColumn != iColumn)
                rounded_solution[jColumn] = 0.0;
        }
        violation = constraint_violation( rounded_solution );
        if( violation < min_violation)
        {
            min_violation = violation;
            int i = 0;
            for( std::vector<int>::iterator vit = sos_inds.begin(); vit != sos_inds.end(); vit++ )
            {
                int jColumn = *vit;
                best_sol[i] = rounded_solution[jColumn];
                i++;
            }
        }
    }
    int i = 0;
    for( std::vector<int>::iterator it = sos_inds.begin(); it != sos_inds.end(); it++ )
    {
        int iColumn = *it;
        rounded_solution[iColumn] = best_sol[i];
        i++;
    }
}


bool RoundingOFP::has_sos_constraints() const
{
    if( sos_columns.size() >= 1)
        return true;
    return false;
}

}





