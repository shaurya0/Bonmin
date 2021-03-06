// (C) Copyright International Business Machines Corporation, 2006
// All Rights Reserved.
// This code is published under the Common Public License.
//
// Authors :
// Pierre Bonami, International Business Machines Corporation
//
// Date : 11/21/2006




#ifndef BonFilterWarmStart_H
#define BonFilterWarmStart_H

#include "CoinWarmStartBasis.hpp"
#include "CoinWarmStartPrimalDual.hpp"
#include "BonFilterTypes.hpp" /* for types */
#include "CoinSmartPtr.hpp"

#include <vector>

namespace Bonmin
{

  /** Warm start for filter interface.
   * Warm start for filter constists of a (possibly huge) array of integers.
   * This class inherits from CoinWarmStartPrimalDual, because that's what
   * this warmstart really is. <br>
   * For practical reason (integration in Cbc) this class also inherits from
   * CoinWarmStartBasis. <br>
  */
  class FilterWarmStart :
    public virtual CoinWarmStartPrimalDual, public virtual CoinWarmStartBasis,
    public Coin::ReferencedObject
  {
    typedef FilterTypes::fint fint;
    typedef FilterTypes::real real;

  public:
    /** Default values for istat */
    static fint def_istat[14];
    /** Constructor */
    FilterWarmStart(const fint xSize = 0,
        const real* xArray = NULL,
        const fint lamSize = 0,
        const real* lamArray = NULL,
        const fint lwsSize = 0,
        const fint *lwsArray = NULL,
        const fint istat[14] = def_istat);

    /** Copy constructor */
    FilterWarmStart(const FilterWarmStart & other);

    /** virtual copy */
    virtual CoinWarmStart * clone() const
    {
      return new FilterWarmStart(*this);
    }

    /** Destructor. */
    virtual ~FilterWarmStart();

    /** Generate differences.*/
    virtual CoinWarmStartDiff* generateDiff(const CoinWarmStart * const other) const;

    /** Apply differences. */
    virtual void applyDiff(const CoinWarmStartDiff * const cswDiff);

    /** Access to lws array */
    const fint *lwsArray() const
    {
      return lwsArray_;
    }

    /** Access to lws size. */
    fint lwsSize() const
    {
      return lwsSize_;
    }

    const fint* istat()const
    {
      return istat_;
    }

    /// flush the starting point
    void flushPoint();

    ///Is this an empty warm start?
    bool empty() const
    {
      return empty_;
    }
  private:
    /** Size of fint lws array store. */
    fint lwsSize_;

    /** fint lws array to store */
    fint* lwsArray_;

    /** Filter's istat (AW: I think we only need first entry) */
    fint istat_[14];

    ///Say if warm start is empty
    bool empty_;
  };

  class FilterWarmStartDiff : public CoinWarmStartPrimalDualDiff
  {
    typedef FilterTypes::fint fint;
    typedef FilterTypes::real real;

    friend class FilterWarmStart;

  public:
    FilterWarmStartDiff(CoinWarmStartPrimalDualDiff * diff,
			fint capacity);

    virtual ~FilterWarmStartDiff() {}

    virtual CoinWarmStartDiff * clone() const
    {
      return new FilterWarmStartDiff(*this);
    }

    void flushPoint();
  private:
    /** One difference is two integers (indice and difference). */
    typedef std::pair<fint, fint> OneDiff;
    /** Vector of all the differences.*/
    std::vector<OneDiff> differences;
    /** istat */
    fint istat_[14];
  };

} /* end namespace Bonmin */
#endif

