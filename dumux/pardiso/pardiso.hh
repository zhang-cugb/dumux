// $Id$

#ifndef DUNE_PARDISO_HH
#define DUNE_PARDISO_HH

#include <dune/istl/preconditioners.hh>

/* Change this, if your Fortran compiler does not append underscores. */
/* e.g. the AIX compiler:  #define F77_FUNC(func) func                */

#ifdef AIX
#define F77_FUNC(func)  func
#else
#define F77_FUNC(func)  func ## _
#endif


#ifdef HAVE_PARDISO
/* PARDISO prototype. */
extern "C" int F77_FUNC(pardisoinit)
    (void *, int *, int *);

extern "C" int F77_FUNC(pardiso)
    (void *, int *, int *, int *, int *, int *,
     double *, int *, int *, int *, int *, int *,
     int *, double *, double *, int *);
#endif

namespace Dune {


  /*! \brief The sequential Pardiso preconditioner.

     Put the Pardiso direct solver into the preconditioner framework.
   */
  template<class M, class X, class Y>
  class SeqPardiso : public Preconditioner<X,Y> {
  public:
    //! \brief The matrix type the preconditioner is for.
    typedef M matrix_type;
    //! \brief The domain type of the preconditioner.
    typedef X domain_type;
    //! \brief The range type of the preconditioner.
    typedef Y range_type;
    //! \brief The field type of the preconditioner.
    typedef typename X::field_type field_type;

    typedef typename M::RowIterator RowIterator;
    typedef typename M::ColIterator ColIterator;

    // define the category
    enum {
      //! \brief The category the preconditioner is part of
      category=SolverCategory::sequential
    };


    SeqPardiso()
    {

#ifdef HAVE_PARDISO
            mtype_ = 11;
            nrhs_ = 1;
            num_procs_ = 1;
            maxfct_ = 1;
            mnum_   = 1;
            msglvl_ = 0;
            error_  = 0;

            //F77_FUNC(pardisoinit) (pt_,  &mtype_, iparm_);
#else
        DUNE_THROW(NotImplemented, "no Pardiso library available, reconfigure with correct --with-pardiso options");
#endif
        }

    void factorize (M& A)
    {
#ifdef HAVE_PARDISO

    RowIterator i0 = A.begin();
    ColIterator j0 = (*i0).begin();

    systemsize_ = (*j0).N();
        n_ = A.N()*systemsize_;
        int nnz = 0;
        RowIterator endi = A.end();
    int rows = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
        rows++;
            //if (A.rowdim(i.index()) != 1)
            //    DUNE_THROW(NotImplemented, "SeqPardiso: row blocksize != 1.");
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
                //if (A.coldim(j.index()) != 1)
                //    DUNE_THROW(NotImplemented, "SeqPardiso: column blocksize != 1.");
                nnz += systemsize_*systemsize_;
            }
        }
          //std::cout << "rows = " << rows;
        std::cout << "SeqPardiso: dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;

        a_ = new double[nnz];
        ia_ = new int[n_+1];
        ja_ = new int[nnz];

        int count = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
        for (int iComp = 0; iComp < systemsize_; iComp++) {
            ia_[i.index()*systemsize_ + iComp] = count+1;
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
            for (int jComp = 0; jComp < systemsize_; jComp++) {
                a_[count] = (*j)[iComp][jComp];
                ja_[count] = j.index()*systemsize_ + jComp + 1;

                count++;
        }
            }
        }
        }
        ia_[n_] = count+1;

    /*std::cout << "systemsize_ =" << systemsize_ << ", n_ = " << n_ << ", nnz_ = " << nnz << std::endl;
    for (int i = 0; i <= n_; i++)
        std::cout << ia_[i] << std::endl;
    */
    /*std::cout << "ja_:" << std::endl;
    for (int i = 0; i <= nnz; i++)
        std::cout << ja_[i] << std::endl;
    std::cout << "a_:" << std::endl;
    for (int i = 0; i <= nnz; i++)
        std::cout << a_[i] << std::endl;
    */

           F77_FUNC(pardisoinit) (pt_,  &mtype_, iparm_);


        phase_ = 11;
        int idum;
        double ddum;
        iparm_[2]  = num_procs_;

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                   &n_, a_, ia_, ja_, &idum, &nrhs_,
                   iparm_, &msglvl_, &ddum, &ddum, &error_);

        if (error_ != 0)
            DUNE_THROW(MathError, "Constructor SeqPardiso: Reordering failed. Error code " << error_);

        std::cout << "  Reordering completed. Number of nonzeros in factors  = " << iparm_[17] << std::endl;

        phase_ = 22;

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                   &n_, a_, ia_, ja_, &idum, &nrhs_,
                   iparm_, &msglvl_, &ddum, &ddum, &error_);

        if (error_ != 0)
            DUNE_THROW(MathError, "Constructor SeqPardiso: Factorization failed. Error code " << error_);

        std::cout << "  Factorization completed." << std::endl;




#endif
    }

    /*! \brief Constructor.

    Constructor gets all parameters to operate the prec.
    \param A The matrix to operate on.
    \param n The number of iterations to perform.
    \param w The relaxation factor.
    */
    SeqPardiso (M& A)

    {
#ifdef HAVE_PARDISO

        mtype_ = 11;
        nrhs_ = 1;
        num_procs_ = 1;
        maxfct_ = 1;
        mnum_   = 1;
        msglvl_ = 0;
        error_  = 0;

    RowIterator i0 = A.begin();
    ColIterator j0 = (*i0).begin();

    systemsize_ = (*j0).N();
        n_ = A.N()*systemsize_;



        int nnz = 0;
        RowIterator endi = A.end();
    int rows = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
        rows++;
            //if (A.rowdim(i.index()) != 1)
            //    DUNE_THROW(NotImplemented, "SeqPardiso: row blocksize != 1.");
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
                //if (A.coldim(j.index()) != 1)
                //    DUNE_THROW(NotImplemented, "SeqPardiso: column blocksize != 1.");
                nnz += systemsize_*systemsize_;
            }
        }
          //std::cout << "rows = " << rows;
        std::cout << "SeqPardiso: dimension = " << n_ << ", number of nonzeros = " << nnz << std::endl;

        a_ = new double[nnz];
        ia_ = new int[n_+1];
        ja_ = new int[nnz];

        int count = 0;
        for (RowIterator i = A.begin(); i != endi; ++i)
        {
        for (int iComp = 0; iComp < systemsize_; iComp++) {
            ia_[i.index()*systemsize_ + iComp] = count+1;
            ColIterator endj = (*i).end();
            for (ColIterator j = (*i).begin(); j != endj; ++j) {
            for (int jComp = 0; jComp < systemsize_; jComp++) {
                a_[count] = (*j)[iComp][jComp];
                ja_[count] = j.index()*systemsize_ + jComp + 1;

                count++;
        }
            }
        }
        }
        ia_[n_] = count+1;

    /*std::cout << "systemsize_ =" << systemsize_ << ", n_ = " << n_ << ", nnz_ = " << nnz << std::endl;
    for (int i = 0; i <= n_; i++)
        std::cout << ia_[i] << std::endl;
    */
    /*std::cout << "ja_:" << std::endl;
    for (int i = 0; i <= nnz; i++)
        std::cout << ja_[i] << std::endl;
    std::cout << "a_:" << std::endl;
    for (int i = 0; i <= nnz; i++)
        std::cout << a_[i] << std::endl;
    */

        F77_FUNC(pardisoinit) (pt_,  &mtype_, iparm_);

        phase_ = 11;
        int idum;
        double ddum;
        iparm_[2]  = num_procs_;

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                   &n_, a_, ia_, ja_, &idum, &nrhs_,
                   iparm_, &msglvl_, &ddum, &ddum, &error_);

        if (error_ != 0)
            DUNE_THROW(MathError, "Constructor SeqPardiso: Reordering failed. Error code " << error_);

        std::cout << "  Reordering completed. Number of nonzeros in factors  = " << iparm_[17] << std::endl;

        phase_ = 22;

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                   &n_, a_, ia_, ja_, &idum, &nrhs_,
                   iparm_, &msglvl_, &ddum, &ddum, &error_);

        if (error_ != 0)
            DUNE_THROW(MathError, "Constructor SeqPardiso: Factorization failed. Error code " << error_);

        std::cout << "  Factorization completed." << std::endl;

#else
        DUNE_THROW(NotImplemented, "no Pardiso library available, reconfigure with correct --with-pardiso options");
#endif
    }

    /*!
      \brief Prepare the preconditioner.

      \copydoc Preconditioner::pre(X&,Y&)
    */
    virtual void pre (X& x, Y& b) {}

    /*!
      \brief Apply the preconditioner.

      \copydoc Preconditioner::apply(X&,const Y&)
    */
    virtual void apply (X& v, const Y& d)
    {
#ifdef HAVE_PARDISO
        phase_ = 33;

        iparm_[7] = 1;       /* Max numbers of iterative refinement steps. */
        int idum;
        double x[2*n_];
        double b[2*n_];
        for (typename X::size_type i = 0; i < v.size(); i++) {
        for (int comp = 0; comp < systemsize_; comp++) {
                x[i*systemsize_ + comp] = v[i][comp];
                b[i*systemsize_ + comp] = d[i][comp];
        }
        }

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                   &n_, a_, ia_, ja_, &idum, &nrhs_,
                   //&n_, &ddum, &idum, &idum, &idum, &nrhs_,
                   iparm_, &msglvl_, b, x, &error_);

        if (error_ != 0)
            DUNE_THROW(MathError, "SeqPardiso.apply: Backsolve failed. Error code " << error_);

        for (typename X::size_type i = 0; i < v.size(); i++)
        for (int comp = 0; comp < systemsize_; comp++)
                v[i][comp] = x[i*systemsize_ + comp];


        //phase_ = -1;                 // Release internal memory.
        //int idum;
        //double ddum;

        //F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
        //           &n_, &ddum, ia_, ja_, &idum, &nrhs_,
        //           iparm_, &msglvl_, &ddum, &ddum, &error_);
    //delete a_;
    //delete ia_;
    //delete ja_;
        //std::cout << "SeqPardiso: Backsolve completed." << std::endl;
#endif
    }

    /*!
      \brief Clean up.

      \copydoc Preconditioner::post(X&)
    */
    virtual void post (X& x)
    {
#ifdef HAVE_PARDISO
       phase_ = -1;                 // Release internal memory.
        int idum;
        double ddum;

        F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                   &n_, &ddum, ia_, ja_, &idum, &nrhs_,
                   iparm_, &msglvl_, &ddum, &ddum, &error_);
    delete a_;
    delete ia_;
    delete ja_;
#endif
    }

    ~SeqPardiso()
    {
#ifdef HAVE_PARDISO
        if (phase_ != -1) {
            phase_ = -1;                 // Release internal memory.
            int idum;
            double ddum;

            F77_FUNC(pardiso) (pt_, &maxfct_, &mnum_, &mtype_, &phase_,
                    &n_, &ddum, ia_, ja_, &idum, &nrhs_,
                    iparm_, &msglvl_, &ddum, &ddum, &error_);
            delete a_;
            delete ia_;
            delete ja_;
        }
#endif
    }

  private:
    //M A_; //!< The matrix we operate on.
    int n_; //!< dimension of the system
    double *a_; //!< matrix values
    int *ia_; //!< indices to rows
    int *ja_; //!< column indices
    int mtype_; //!< matrix type, currently only 11 (real unsymmetric matrix) is supported
    int nrhs_; //!< number of right hand sides
    void *pt_[64]; //!< internal solver memory pointer
    int iparm_[64]; //!< Pardiso control parameters.
    int maxfct_;    //!< Maximum number of numerical factorizations.
    int mnum_;  //!<        Which factorization to use.
    int msglvl_;    //!< flag to print statistical information
    int error_;      //!< error flag
    int num_procs_; //!< number of processors.
    int systemsize_;
    int phase_;
  };

}






#endif

