// $Id$
#ifndef DUNE_IMPESSUBPROBS_HH
#define DUNE_IMPESSUBPROBS_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/upscaledsaturation/preprocess/fractionalflowsubprobs.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, last changed by Markus Wolff
 */

namespace Dune {
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * coupled diffusion/transport problems
 */

template<class G, class DiffusionSubProbs, class TransportSubProbs, class VC> class IMPESSubProbs :
    public FractionalFlowSubProbs<G, DiffusionSubProbs, TransportSubProbs, VC> {
    typedef typename DiffusionSubProbs::RepresentationType PressType;
    typedef typename DiffusionSubProbs::NumberType RT;

public:
    typedef typename TransportSubProbs::RepresentationType RepresentationType;

    virtual void totalVelocity(const RT t=0) {
        this->calcTotalVelocity(t);
//        std::cout<<"velocity= "<<this->variables().velocity<<std::endl;
        return;
    }

    VC variables() const{
        return this->transproblem.variables;
    }

    virtual void initial() {
        double t = 0;
        this->initialTransport();
        this->pressure(t);
        totalVelocity(t);

        return;
    }

    virtual int update(const RT t, RT& dt, RepresentationType& updateVec,
            RT cFLFactor = 1) {
        int pressSize = variables().pressure.size();
        PressType pressOldIter(variables().pressure);
        PressType pressHelp(pressSize);
        int satSize = variables().saturation.size();
        RepresentationType saturation(variables().saturation);
        RepresentationType satOldIter(variables().saturation);
        RepresentationType satHelp(satSize);
        RepresentationType satDiff(satSize);
        RepresentationType updateOldIter(satSize);
        RepresentationType updateHelp(satSize);
        RepresentationType updateDiff(satSize);

        bool converg = false;
        int iter = 0;
        int iterTot = 0;
        updateOldIter = 0;
        while (!converg) {
            iter++;
            iterTot++;
            // update pressure
            pressure(t);
            totalVelocity(t);

            TransportSubProbs::update(t, dt, updateVec,cFLFactor);
            if (iterFlag) { // only needed if iteration has to be done
                variables().pressure *= omega;
                pressHelp = pressOldIter;
                pressHelp *= (1-omega);
                variables().pressure += pressHelp;
                pressOldIter = variables().pressure;

                updateHelp = updateVec;
                saturation = variables().saturation;
                saturation += (updateHelp *= (dt*cFLFactor));
                saturation *= omega;
                satHelp = satOldIter;
                satHelp *= (1-omega);
                saturation += satHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                satOldIter = saturation;
                updateOldIter = updateVec;
            }
            // break criteria for iteration loop
            if (iterFlag==2&& dt*updateDiff.two_norm()/(saturation).two_norm() <= maxDefect )
                converg = true;
            else if (iterFlag==2&& iter > nIter ) {
                std::cout << "Nonlinear loop in IMPES.update exceeded nIter = "
                        << nIter << " iterations."<< std::endl;
                return 1;
            } else if (iterFlag==1&& iter > nIter )
                converg = true;
            else if (iterFlag==0)
                converg = true;
        }
        // outputs
        if (iterFlag==2)
//            std::cout << "Iteration steps: "<< iterTot << std::endl;
//        std::cout.setf(std::ios::scientific, std::ios::floatfield);

        return 0;
    }

    //! Construct an IMPES object.
    IMPESSubProbs(DiffusionSubProbs& diff, TransportSubProbs& trans, int flag = 2, int nIt = 30,
            double maxDef = 1e-5, double om = 1) :
        FractionalFlowSubProbs<G, DiffusionSubProbs, TransportSubProbs, VC>(diff, trans),
                iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om) {
    }

protected:
    const int iterFlag;
    const int nIter;
    const double maxDefect;
    const double omega;
};
}
#endif
