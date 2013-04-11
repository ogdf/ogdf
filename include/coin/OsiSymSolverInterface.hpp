/*===========================================================================*/
/*                                                                           */
/* This file is part of the SYMPHONY Branch, Cut, and Price Callable         */
/* Library.                                                                  */
/*                                                                           */
/* SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and    */
/* Laci Ladanyi (ladanyi@us.ibm.com).                                        */
/*                                                                           */
/* (c) Copyright 2004-2006 Ted Ralphs and Lehigh University.                 */
/* All Rights Reserved.                                                      */
/*                                                                           */
/* The authors of this file are Menal Guzelsoy and Ted Ralphs                */
/*                                                                           */
/* This software is licensed under the Eclipse Public License. Please see    */
/* accompanying file for terms.                                              */
/*                                                                           */
/*===========================================================================*/

#ifndef OsiSymSolverInterface_hpp
#define OsiSymSolverInterface_hpp

#include "OsiSolverInterface.hpp"
#include "OsiSymSolverParameters.hpp"
#include "SymWarmStart.hpp"

#include <string>

typedef struct SYM_ENVIRONMENT sym_environment;

//#############################################################################

/** OSI Solver Interface for SYMPHONY

  Many OsiSolverInterface query methods return a const pointer to the
  requested read-only data. If the model data is changed or the solver
  is called, these pointers may no longer be valid and should be 
  refreshed by invoking the member function to obtain an updated copy
  of the pointer.
  For example:
  \code
      OsiSolverInterface solverInterfacePtr ;
      const double * ruBnds = solverInterfacePtr->getRowUpper();
      solverInterfacePtr->applyCuts(someSetOfCuts);
      // ruBnds is no longer a valid pointer and must be refreshed
      ruBnds = solverInterfacePtr->getRowUpper();
  \endcode

  Querying a problem that has no data associated with it will result in
  zeros for the number of rows and columns, and NULL pointers from
  the methods that return vectors.
*/

class OsiSymSolverInterface : virtual public OsiSolverInterface {
   friend void OsiSymSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir);
   
public:
   ///@name Solve methods 
   //@{
   /// Solve initial LP relaxation
   virtual void initialSolve();

   /// Resolve an IP problem modification
   virtual void resolve();

   /// Invoke solver's built-in enumeration algorithm
   virtual void branchAndBound();
   
   /// Invoke solver's multi-criteria enumeration algorithm
   virtual void multiCriteriaBranchAndBound();

   /// Get a lower bound for the new rhs problem using the warm start tree.
   virtual double getLbForNewRhs(int cnt, int *index, 
				 double * value);
   /// Get an upper bound for the new rhs problem using the warm start tree.
   virtual double getUbForNewRhs(int cnt, int *index, 
				 double * value);
#if 0
   /// Get a lower bound for the new obj problem using the warm start tree.
   virtual double getLbForNewObj(int cnt, int *index, 
				 double * value);
   /// Get an upper bound for the new obj problem using the warm start tree.
#endif
   virtual double getUbForNewObj(int cnt, int *index, 
				 double * value);

  //@}

  //---------------------------------------------------------------------------
  /**@name Parameter set/get methods

     The set methods return true if the parameter was set to the given value,
     false otherwise. There can be various reasons for failure: the given
     parameter is not applicable for the solver (e.g., refactorization
     frequency for the volume algorithm), the parameter is not yet implemented
     for the solver or simply the value of the parameter is out of the range
     the solver accepts. If a parameter setting call returns false check the
     details of your solver.

     The get methods return true if the given parameter is applicable for the
     solver and is implemented. In this case the value of the parameter is
     returned in the second argument. Otherwise they return false.

  */
  //@{
    // Set an integer parameter
    virtual bool setIntParam(OsiIntParam key, int value);

    // Set SYMPHONY int parameter
    virtual bool setSymParam(OsiSymIntParam key, int value);

    // Set SYMPHONY int parameter directly by the C interface parameter name
    virtual bool setSymParam(const std::string key, int value);


    // Set an double parameter
    virtual bool setDblParam(OsiDblParam key, double value);

    // Set SYMPHONY double parameter
    virtual bool setSymParam(OsiSymDblParam key, double value);

    // Set SYMPHONY double parameter directly by the C interface parameter name
    virtual bool setSymParam(const std::string key, double value);



    // Set a string parameter
    virtual bool setStrParam(OsiStrParam key, const std::string & value);

    // Set SYMPHONY string parameter
    virtual bool setSymParam(OsiSymStrParam key, const std::string & value);

    // Set SYMPHONY string parameter directly by the C interface parameter name
    virtual bool setSymParam(const std::string key, const std::string value);



    // Get an integer parameter
    virtual bool getIntParam(OsiIntParam key, int& value) const;

    // Get SYMPHONY int parameter
    virtual bool getSymParam(OsiSymIntParam key, int& value) const;

    // Get SYMPHONY int parameter directly by the C interface parameter name
    virtual bool getSymParam(const std::string key, int& value) const;


    // Get an double parameter
    virtual bool getDblParam(OsiDblParam key, double& value) const;

    // Get SYMPHONY double parameter
    virtual bool getSymParam(OsiSymDblParam key, double& value) const;

    // Get SYMPHONY double parameter directly by the C interface parameter name
    virtual bool getSymParam(const std::string key, double& value) const;



    // Get a string parameter
    virtual bool getStrParam(OsiStrParam key, std::string& value) const;

    // Get SYMPHONY string parameter
    virtual bool getSymParam(OsiSymStrParam key, std::string& value) const;

    // Get SYMPHONY string parameter directly by the C interface parameter name
    virtual bool getSymParam(const std::string key, std::string& value) const;

  //@}

  //---------------------------------------------------------------------------
  ///@name Methods returning info on how the solution process terminated
  //@{
    /// Are there numerical difficulties?
   virtual bool isAbandoned() const;

    /// Is optimality proven?
   virtual bool isProvenOptimal() const;

    /// Is primal infeasiblity proven?
   virtual bool isProvenPrimalInfeasible() const;

    /// Is dual infeasiblity proven?
   virtual bool isProvenDualInfeasible() const {
      throw CoinError("Error: Function not implemented",
		      "isProvenDualInfeasible", "OsiSymSolverInterface");
   }
    /// Is the given primal objective limit reached?
   //virtual bool isPrimalObjectiveLimitReached() const;

    /// Is the given dual objective limit reached?
    //virtual bool isDualObjectiveLimitReached() const{
    //   throw CoinError("Error: Function not implemented",
		//       "isDualObjectiveLimitReached", "OsiSymSolverInterface");
    //}
    /// Iteration limit reached?
   virtual bool isIterationLimitReached() const;

   /// Time limit reached?
   virtual bool isTimeLimitReached() const;

   /// Target gap achieved?
   virtual bool isTargetGapReached() const;

  //@}

  //---------------------------------------------------------------------------
  /**@name Warm start methods */
  //@{
    /*! \brief Get an empty warm start object
      
      This routine returns an empty warm start object. Its purpose is
      to provide a way to give a client a warm start object of the
      appropriate type, which can resized and modified as desired.
    */

    virtual CoinWarmStart *getEmptyWarmStart () const{
       throw CoinError("Error: Function not implemented",
		       "getEmptyWarmStart", "OsiSymSolverInterface");
    }

    /** Get warm start information.

      If there is no valid solution, an empty warm start object (0 rows, 0
      columns) wil be returned.
    */

   /* 
      virtual CoinWarmStart* getWarmStart(bool keepTreeInSymEnv = false) const;
   */

    virtual CoinWarmStart* getWarmStart() const;

    /** Set warm start information.
    
      Return true/false depending on whether the warm start information was
      accepted or not. */
   virtual bool setWarmStart(const CoinWarmStart* warmstart);
   //@}
   
  //---------------------------------------------------------------------------
   /**@name Problem query methods
      
   Querying a problem that has no data associated with it will result in
   zeros for the number of rows and columns, and NULL pointers from
   the methods that return vectors.
   
   Const pointers returned from any data-query method are valid as
   long as the data is unchanged and the solver is not called.
   */
   //@{
   /// Get pointer to SYMPHONY environment (eventually we won't need this)
   sym_environment *getSymphonyEnvironment() const {return env_;}
   
   /// Get number of columns
   virtual int getNumCols() const;
   
   /// Get number of rows
   virtual int getNumRows() const;
  
   /// Get number of nonzero elements
   virtual int getNumElements() const;
  
   /// Get pointer to array[getNumCols()] of column lower bounds
   virtual const double * getColLower() const;
  
   /// Get pointer to array[getNumCols()] of column upper bounds
   virtual const double * getColUpper() const;
  
      /** Get pointer to array[getNumRows()] of row constraint senses.
  	<ul>
  	<li>'L': <= constraint
  	<li>'E': =  constraint
  	<li>'G': >= constraint
  	<li>'R': ranged constraint
  	<li>'N': free constraint
  	</ul>
      */
   virtual const char * getRowSense() const;
  
      /** Get pointer to array[getNumRows()] of row right-hand sides
  	<ul>
  	  <li> if getRowSense()[i] == 'L' then
	       getRightHandSide()[i] == getRowUpper()[i]
  	  <li> if getRowSense()[i] == 'G' then
	       getRightHandSide()[i] == getRowLower()[i]
  	  <li> if getRowSense()[i] == 'R' then
	       getRightHandSide()[i] == getRowUpper()[i]
  	  <li> if getRowSense()[i] == 'N' then
	       getRightHandSide()[i] == 0.0
  	</ul>
      */
   virtual const double * getRightHandSide() const;
  
      /** Get pointer to array[getNumRows()] of row ranges.
  	<ul>
            <li> if getRowSense()[i] == 'R' then
                    getRowRange()[i] == getRowUpper()[i] - getRowLower()[i]
            <li> if getRowSense()[i] != 'R' then
                    getRowRange()[i] is 0.0
          </ul>
      */
   virtual const double * getRowRange() const;
  
      /// Get pointer to array[getNumRows()] of row lower bounds
   virtual const double * getRowLower() const;
  
      /// Get pointer to array[getNumRows()] of row upper bounds
   virtual const double * getRowUpper() const;
  
      /// Get pointer to array[getNumCols()] of objective function coefficients
   virtual const double * getObjCoefficients() const;
  
   /** Get pointer to array[getNumCols()] of second 
       objective function coefficients if loaded before.
   */
   virtual const double * getObj2Coefficients() const;

      /// Get objective function sense (1 for min (default), -1 for max)
   virtual double getObjSense() const; 
 
      /// Return true if variable is continuous
   virtual bool isContinuous(int colIndex) const;
  
      /// Return true if variable is binary
   virtual bool isBinary(int colIndex) const;  

      /** Return true if column is integer.
          Note: This function returns true if the the column
          is binary or a general integer.
      */
   virtual bool isInteger(int colIndex) const;
  
      /// Return true if variable is general integer
   virtual bool isIntegerNonBinary(int colIndex) const;
  
      /// Return true if variable is binary and not fixed at either bound
   virtual bool isFreeBinary(int colIndex) const; 
   
      /// Get pointer to row-wise copy of matrix
   virtual const CoinPackedMatrix * getMatrixByRow() const;
  
      /// Get pointer to column-wise copy of matrix
   virtual const CoinPackedMatrix * getMatrixByCol() const;  

      /// Get solver's value for infinity
      virtual double getInfinity() const;

    //@}
    
    /**@name Solution query methods */
    //@{
      /// Get pointer to array[getNumCols()] of primal variable values
   virtual const double * getColSolution() const;
  
      /// Get pointer to array[getNumRows()] of dual variable values
      virtual const double * getRowPrice() const;
  
      /// Get a pointer to array[getNumCols()] of reduced costs
      virtual const double * getReducedCost() const;
  
      /** Get pointer to array[getNumRows()] of row activity levels (constraint
  	matrix times the solution vector). */
   virtual const double * getRowActivity() const;
  
      /// Get objective function value
   virtual double getObjValue() const;

      /// Get the current upper/lower bound
   virtual double getPrimalBound() const;

      /** Get the number of iterations it took to solve the problem (whatever
	  ``iteration'' means to the solver). */
   virtual int getIterationCount() const;
  
      /** Get as many dual rays as the solver can provide. In case of proven
          primal infeasibility there should be at least one.
     
          \note
	  Implementors of solver interfaces note that
          the double pointers in the vector should point to arrays of length
          getNumRows() and they should be allocated via new[].
     
          \note
	  Clients of solver interfaces note that
          it is the client's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getDualRays(int maxNumRays,
					       bool fullRay = false) const{
       throw CoinError("Error: Function not implemented",
		       "getDualRays", "OsiSymSolverInterface");
    }
      /** Get as many primal rays as the solver can provide. (In case of proven
          dual infeasibility there should be at least one.)
     
          <strong>NOTE for implementers of solver interfaces:</strong> <br>
          The double pointers in the vector should point to arrays of length
          getNumCols() and they should be allocated via new[]. <br>
     
          <strong>NOTE for users of solver interfaces:</strong> <br>
          It is the user's responsibility to free the double pointers in the
          vector using delete[].
      */
      virtual std::vector<double*> getPrimalRays(int maxNumRays) const{
       throw CoinError("Error: Function not implemented",
		       "getPrimalRays", "OsiSymSolverInterface");
    }
  
    //@}

    //-------------------------------------------------------------------------
    /**@name Methods to modify the objective, bounds, and solution

       For functions which take a set of indices as parameters
       (\c setObjCoeffSet(), \c setColSetBounds(), \c setRowSetBounds(),
       \c setRowSetTypes()), the parameters follow the C++ STL iterator
       convention: \c indexFirst points to the first index in the
       set, and \c indexLast points to a position one past the last index
       in the set.
    
    */
    //@{
      /** Set an objective function coefficient */
   virtual void setObjCoeff( int elementIndex, double elementValue );

      /** Set an objective function coefficient for the second objective */
   virtual void setObj2Coeff( int elementIndex, double elementValue );

   using OsiSolverInterface::setColLower ;
      /** Set a single column lower bound.
    	  Use -getInfinity() for -infinity. */
   virtual void setColLower( int elementIndex, double elementValue );
      
   using OsiSolverInterface::setColUpper ;
      /** Set a single column upper bound.
    	  Use getInfinity() for infinity. */
   virtual void setColUpper( int elementIndex, double elementValue );      

       /** Set a single row lower bound.
    	  Use -getInfinity() for -infinity. */
   virtual void setRowLower( int elementIndex, double elementValue );
      
      /** Set a single row upper bound.
    	  Use getInfinity() for infinity. */
   virtual void setRowUpper( int elementIndex, double elementValue );
    
      /** Set the type of a single row */
      virtual void setRowType(int index, char sense, double rightHandSide,
    			      double range);
    
    /// Set the objective function sense.
    /// (1 for min (default), -1 for max)
   virtual void setObjSense(double s);
    
    /** Set the primal solution variable values
    
	colsol[getNumCols()] is an array of values for the primal variables.
	These values are copied to memory owned by the solver interface object
	or the solver.  They will be returned as the result of getColSolution()
	until changed by another call to setColSolution() or by a call to any
	solver routine.  Whether the solver makes use of the solution in any
	way is solver-dependent.
    */
   virtual void setColSolution(const double *colsol);
   
    /** Set the a priori upper/lower bound */

   virtual void setPrimalBound(const double bound);

    /** Set dual solution variable values

	rowprice[getNumRows()] is an array of values for the dual
	variables. These values are copied to memory owned by the solver
	interface object or the solver.  They will be returned as the result of
	getRowPrice() until changed by another call to setRowPrice() or by a
	call to any solver routine.  Whether the solver makes use of the
	solution in any way is solver-dependent.
    */

   virtual void setRowPrice(const double * rowprice);

    //@}

    //-------------------------------------------------------------------------
    /**@name Methods to set variable type */
    //@{

      using OsiSolverInterface::setContinuous ;
      /** Set the index-th variable to be a continuous variable */
   virtual void setContinuous(int index);

      using OsiSolverInterface::setInteger ;
      /** Set the index-th variable to be an integer variable */
   virtual void setInteger(int index);


      using OsiSolverInterface::setColName ;
   virtual void setColName(char **colname);

    //@}
    //-------------------------------------------------------------------------
    
    //-------------------------------------------------------------------------
    /**@name Methods to expand a problem.

       Note that new columns are added as continuous variables.

    */
    //@{

      using OsiSolverInterface::addCol ;
      /** Add a column (primal variable) to the problem. */
      virtual void addCol(const CoinPackedVectorBase& vec,
			  const double collb, const double colub,   
			  const double obj);

      /** Remove a set of columns (primal variables) from the problem.  */
   virtual void deleteCols(const int num, const int * colIndices);
    
      using OsiSolverInterface::addRow ;
      /** Add a row (constraint) to the problem. */
   virtual void addRow(const CoinPackedVectorBase& vec,
		       const double rowlb, const double rowub);
      /** */
   virtual void addRow(const CoinPackedVectorBase& vec,
		       const char rowsen, const double rowrhs,   
		       const double rowrng);
   
   /** Delete a set of rows (constraints) from the problem. */
   virtual void deleteRows(const int num, const int * rowIndices);
    
    //@}

  //---------------------------------------------------------------------------

  /**@name Methods to input a problem */
  //@{

    virtual void loadProblem();
   
    /** Load in an problem by copying the arguments (the constraints on the
        rows are given by lower and upper bounds). If a pointer is 0 then the
        following values are the default:
        <ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
          <li> <code>rowub</code>: all rows have upper bound infinity
          <li> <code>rowlb</code>: all rows have lower bound -infinity
	  <li> <code>obj</code>: all variables have 0 objective coefficient
        </ul>
    */
    virtual void loadProblem(const CoinPackedMatrix& matrix,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub);
			    
    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by lower and upper bounds).
	For default values see the previous method.

	\warning
	The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       double*& rowlb, double*& rowub);

    /** Load in an problem by copying the arguments (the constraints on the
	rows are given by sense/rhs/range triplets). If a pointer is 0 then the
	following values are the default:
	<ul>
          <li> <code>colub</code>: all columns have upper bound infinity
          <li> <code>collb</code>: all columns have lower bound 0 
	  <li> <code>obj</code>: all variables have 0 objective coefficient
          <li> <code>rowsen</code>: all rows are >=
          <li> <code>rowrhs</code>: all right hand sides are 0
          <li> <code>rowrng</code>: 0 for the ranged rows
        </ul>
    */
    virtual void loadProblem(const CoinPackedMatrix& matrix,
			     const double* collb, const double* colub,
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng);

    /** Load in an problem by assuming ownership of the arguments (the
        constraints on the rows are given by sense/rhs/range triplets). For
        default values see the previous method.

	\warning
	The arguments passed to this method will be
	freed using the C++ <code>delete</code> and <code>delete[]</code>
	functions. 
    */
    virtual void assignProblem(CoinPackedMatrix*& matrix,
			       double*& collb, double*& colub, double*& obj,
			       char*& rowsen, double*& rowrhs,
			       double*& rowrng);

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const CoinBigIndex * start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const double* rowlb, const double* rowub);

    /** Just like the other loadProblem() methods except that the matrix is
	given in a standard column major ordered format (without gaps). */
    virtual void loadProblem(const int numcols, const int numrows,
			     const CoinBigIndex * start, const int* index,
			     const double* value,
			     const double* collb, const double* colub,   
			     const double* obj,
			     const char* rowsen, const double* rowrhs,   
			     const double* rowrng);

    /** Write the problem in MPS format to the specified file.

      If objSense is non-zero, a value of -1.0 causes the problem to be
      written with a maximization objective; +1.0 forces a minimization
      objective. If objSense is zero, the choice is left to implementation.
    */
    virtual void writeMps(const char *filename,
			  const char *extension = "mps",
			  double objSense=0.0) const;

   void parseCommandLine(int argc, char **argv);

   using OsiSolverInterface::readMps ;
   virtual int readMps(const char * infile, const char *extension = "mps");

   virtual int readGMPL(const char * modelFile, const char * dataFile=NULL);

   void findInitialBounds();

   int createPermanentCutPools();

  //@}

  //---------------------------------------------------------------------------

   enum keepCachedFlag {
      /// discard all cached data (default)
      KEEPCACHED_NONE    = 0,
      /// column information: objective values, lower and upper bounds, variable types
      KEEPCACHED_COLUMN  = 1,
      /// row information: right hand sides, ranges and senses, lower and upper bounds for row
      KEEPCACHED_ROW     = 2,
      /// problem matrix: matrix ordered by column and by row
      KEEPCACHED_MATRIX  = 4,
      /// LP solution: primal and dual solution, reduced costs, row activities
      KEEPCACHED_RESULTS = 8,
      /// only discard cached LP solution
      KEEPCACHED_PROBLEM = KEEPCACHED_COLUMN | KEEPCACHED_ROW | KEEPCACHED_MATRIX,
      /// keep all cached data (similar to getMutableLpPtr())
      KEEPCACHED_ALL     = KEEPCACHED_PROBLEM | KEEPCACHED_RESULTS,
      /// free only cached column and LP solution information
      FREECACHED_COLUMN  = KEEPCACHED_PROBLEM & ~KEEPCACHED_COLUMN,
      /// free only cached row and LP solution information
      FREECACHED_ROW     = KEEPCACHED_PROBLEM & ~KEEPCACHED_ROW,
      /// free only cached matrix and LP solution information
      FREECACHED_MATRIX  = KEEPCACHED_PROBLEM & ~KEEPCACHED_MATRIX,
      /// free only cached LP solution information
      FREECACHED_RESULTS = KEEPCACHED_ALL & ~KEEPCACHED_RESULTS
   };
   
  ///@name Constructors and destructors
  //@{
    /// Default Constructor
    OsiSymSolverInterface(); 
    
    /** Clone

      The result of calling clone(false) is defined to be equivalent to
      calling the default constructor OsiSolverInterface().
    */
   virtual OsiSolverInterface * clone(bool copyData = true) const;
  
    /// Copy constructor 
    OsiSymSolverInterface(const OsiSymSolverInterface &);
  
    /// Assignment operator 
    OsiSymSolverInterface & operator=(const OsiSymSolverInterface& rhs);
  
    /// Destructor 
    virtual ~OsiSymSolverInterface ();

    /** Reset the solver interface.

    A call to reset() returns the solver interface to the same state as
    it would have if it had just been constructed by calling the default
    constructor OsiSolverInterface().
    */
    virtual void reset();
  //@}

  //---------------------------------------------------------------------------

protected:
  ///@name Protected methods
  //@{
    /** Apply a row cut (append to the constraint matrix). */
   virtual void applyRowCut( const OsiRowCut & rc );

    /** Apply a column cut (adjust the bounds of one or more variables). */
   virtual void applyColCut( const OsiColCut & cc );

    /** Set OsiSolverInterface object state for default constructor

      This routine establishes the initial values of data fields in the
      OsiSolverInterface object when the object is created using the
      default constructor.
    */

    void setInitialData();
  //@}

private:

  /// The real work of the constructor
  void gutsOfConstructor();
  
  /// The real work of the destructor
  void gutsOfDestructor();

  /// free cached column rim vectors
  void freeCachedColRim();

  /// free cached row rim vectors
  void freeCachedRowRim();

  /// free cached result vectors
  void freeCachedResults();
  
  /// free cached matrices
  void freeCachedMatrix();

  /// free all cached data (except specified entries, see getLpPtr())
  void freeCachedData( int keepCached = KEEPCACHED_NONE );

  /// free all allocated memory
  void freeAllMemory();

  /**@name Private member data */
  //@{
   /// The pointer to the SYMPHONY problem environment
   sym_environment *env_;
  //@}
   
   /// Pointer to objective vector
   mutable double  *obj_;

   /// Pointer to second objective vector to be used in bicriteria solver
   mutable double  *obj2_;
   
   /// Pointer to dense vector of variable lower bounds
   mutable double  *collower_;
   
   /// Pointer to dense vector of variable lower bounds
   mutable double  *colupper_;

   /// Pointer to dense vector of variable lower bounds
   mutable double  *colredcost_;

   /// Pointer to dense vector of row sense indicators
   mutable char    *rowsense_;
  
   /// Pointer to dense vector of row right-hand side values
   mutable double  *rhs_;
  
   /** Pointer to dense vector of slack upper bounds for range constraints 
       (undefined for non-range rows) 
   */
   mutable double  *rowrange_;
   
   /// Pointer to dense vector of row lower bounds
   mutable double  *rowlower_;
   
   /// Pointer to dense vector of row upper bounds
   mutable double  *rowupper_;
   
   /// Pointer to dense vector of row prices
   mutable double  *rowprice_;

   /// Pointer to primal solution vector
   mutable double  *colsol_;
   
   /// Pointer to row activity (slack) vector
   mutable double  *rowact_;
   
   /// Pointer to row-wise copy of problem matrix coefficients.
   mutable CoinPackedMatrix *matrixByRow_;  
   
   /// Pointer to row-wise copy of problem matrix coefficients.
   mutable CoinPackedMatrix *matrixByCol_;  

};

//#############################################################################
/** A function that tests the methods in the OsiSymSolverInterface class. */
void OsiSymSolverInterfaceUnitTest(const std::string & mpsDir, const std::string & netlibDir);

#endif
