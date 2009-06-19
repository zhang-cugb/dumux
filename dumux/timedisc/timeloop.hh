// $Id:$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/

#ifndef DUNE_TIMELOOP_HH
#define DUNE_TIMELOOP_HH

#include "dumux/timedisc/timestep.hh"
#include "dumux/timedisc/rungekuttastep.hh"
#include "dumux/timedisc/impliciteulerstep.hh"


namespace Dune
{

/** \todo Please doc me! */

template<class G, class Model, bool useMultiWriter = false>
class TimeLoop
{
public:
    void execute(Model& model)
    {
        // generate one meta vtk-file holding the individual timesteps
        char multiFileName[128];
        char fileNameVTK[128];
        sprintf(multiFileName, "multi-%s.pvd", fileName);
        std::ofstream multiFile(multiFileName);
        multiFile << "<?xml version=\"1.0\"?>" << std::endl
                  << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" "
                  << std::endl << "compressor=\"vtkZLibDataCompressor\">"
                  << std::endl << " <Collection>" << std::endl;

        // initialize solution with initial values
        model.initial();

        int k = 0;
        model.vtkout(fileName, k);
        switch (G::dimension)
        {
        case 1:
            sprintf(fileNameVTK, "%s-%05d.vtp", fileName, k);
            break;
        default:
            sprintf(fileNameVTK, "%s-%05d.vtu", fileName, k);
            break;
        }
        multiFile << "   <DataSet timestep=\"" << k << "\" file=\""
                  << fileNameVTK << "\"/>" << std::endl;

        // now do the time steps
        double t = tStart;
        //          double dtOriginal;
        //          if (fixed)
        //              dtOriginal = dt;
        while (t < tEnd)
        {
            k++;
            double dtOld = dt;

            if (t == tStart)
                timeStep.execute(model, t, dt, firstDt, tEnd, cFLFactor);
            else
                timeStep.execute(model, t, dt, maxDt, tEnd, cFLFactor);

            if (fixed)
            {
                if (dt > dtOld)
                {
                    t += dtOld;
                    t = std::min(t, tEnd);
                    if (dt > tEnd - t)
                    {
                        std::cout << "\t" << k << "\t" << t << "\t" << dtOld
                                  << "\t # timestep number k, time t, timestep size dt"
                                  << std::endl;
                        //std::cout << ", timestep: " << k << "\t t=" << t << "\t dt=" << (tEnd-t) << std::endl;
                    }
                    else
                    {
                        std::cout << "\t" << k << "\t" << t << "\t" << dtOld
                                  << "\t # timestep number k, time t, timestep size dt"
                                  << std::endl;
                        //std::cout << ", timestep: " << k << "\t t=" << t << "\t dt=" << dt << std::endl;
                    }
                }

                else
                {
                    t += dt;
                    t = std::min(t, tEnd);
                    std::cout << "\t" << k << "\t" << t << "\t" << dt
                              << "\t # timestep number k, time t, timestep size dt"
                              << std::endl;
                    //std::cout << ", timestep: " << k << "\t t=" << t << "\t dt=" << dt << std::endl;
                }
            }
            else
            {
                t += dt;
                t = std::min(t, tEnd);
                std::cout << ", timestep: " << k << "\t t=" << t << "\t dt="
                          << dt << std::endl;

            }

            // generate output
            if (k % modulo == 0)
            {
                model.vtkout(fileName, k / modulo);
                switch (G::dimension)
                {
                case 1:
                    sprintf(fileNameVTK, "%s-%05d.vtp", fileName, k / modulo);
                    break;
                default:
                    sprintf(fileNameVTK, "%s-%05d.vtu", fileName, k / modulo);
                    break;
                }
                multiFile << "   <DataSet timestep=\"" << t << "\" file=\""
                          << fileNameVTK << "\"/>" << std::endl;
            }
            //            if (fixed)
            //                dt = dtOriginal;
        }
        // finalize output
        multiFile << " </Collection>" << std::endl << "</VTKFile>" << std::endl;
        multiFile.close();

        return;
    }

    TimeLoop(const double ts, const double te, const char* name = "timeloop",
             const int mod = 1, const double cfl = 1, const double mdt = 1e100,
             const double fdt = 1e100, TimeStep<G, Model>& tist =
             *(new RungeKuttaStep<G, Model> (1))) :
        tStart(ts), tEnd(te), maxDt(mdt), firstDt(fdt), cFLFactor(cfl), modulo(
                                                                               mod), timeStep(tist), fileName(name), fixed(false)
    {
    }

    TimeLoop(const double ts, const double te, const double dtime = 1e100,
             const char* name = "timeloop", const int mod = 1, const double mdt =
             1e100, const double fdt = 1e100, TimeStep<G, Model>& tist =
             *(new ImplicitEulerStep<G, Model> )) :
        tStart(ts), tEnd(te), dt(dtime), maxDt(mdt), firstDt(fdt),
        cFLFactor(1), modulo(mod), timeStep(tist), fileName(name),
        fixed(true)
    {
    }

private:
    const double tStart;
    const double tEnd;
    double dt;
    const double maxDt;
    const double firstDt;
    const double cFLFactor;
    const int modulo;
    TimeStep<G, Model>& timeStep;
    const char* fileName;
    const bool fixed;
};

/** \todo Please doc me! */

template<class G, class Model>
class TimeLoop<G, Model, true>
{
public:
    // HACK: this function does the same as the execute function
    // above, but it uses the new newton method and the
    // vtkMultiWriter. it is a temorary measure to ease the
    // migration...
    template<class MultiWriter>
    void executeMultiWriter(Model& model, MultiWriter& writer,
                            bool writeRestart = false, bool restart = false, int restartNum = 0)
    {
        int k = 0;
        int countVtk = 0;

        writer.beginTimestep(0, model.grid().leafView());

        // initialize solution with initial values
        model.setVtkMultiWriter(&writer);
        if (!restart)
        {
            model.initial();
        }

        // initialize solution from restart file with the number restartNum
        else
        {
            model.restart(restartNum);
        }
        model.addvtkfields(writer);
        std::cout << ">>> writing initial output file" << std::endl;
        writer.endTimestep(); // writes output file

        // now do the time steps
        double t = tStart;
        //          double dtOriginal;
        //          if (fixed)
        //              dtOriginal = dt;
        while (t < tEnd)
        {
            k++;
            if (t == tStart)
                timeStep.execute(model, t, dt, firstDt, tEnd, cFLFactor);
            else
                timeStep.execute(model, t, dt, maxDt, tEnd, cFLFactor);

            t += dt;
            t = std::min(t, tEnd);
            std::cout << ", timestep: " << k << "\t t=" << t << "\t dt=" << dt
                      << std::endl;

            // generate output
            if (k % modulo == 0)
            {
                countVtk++;
                writer.beginTimestep(t, model.grid().leafView());
                model.addvtkfields(writer);
                std::cout << ">>> writing output-file number " << countVtk
                          << " at time : " << t << std::endl;
                writer.endTimestep();
            }
            // write restart file
            if (writeRestart == true)
            {
                if (k % modulo == 0) // writes restart file each timestep
                {
                    ++restartNum;
                    std::cout << ">>> writing restart-file number "
                              << restartNum << " at time : " << t << std::endl;
                    model.writerestartfile(restartNum);
                }
            }
        }

        return;
    }

    TimeLoop(const double ts, const double te, const char* name = "timeloop",
             const int mod = 1, const double cfl = 1, const double mdt = 1e100,
             const double fdt = 1e100, TimeStep<G, Model>& tist =
             *(new RungeKuttaStep<G, Model> (1))) :
        tStart(ts), tEnd(te), maxDt(mdt), firstDt(fdt), cFLFactor(cfl), modulo(
                                                                               mod), timeStep(tist), fileName(name), fixed(false)
    {
    }

    TimeLoop(const double ts, const double te, const double dtime = 1e100,
             const char* name = "timeloop", const int mod = 1, const double fdt =
             1e100, TimeStep<G, Model>& tist = *(new ImplicitEulerStep<
                                                 G, Model> )) :
        tStart(ts), tEnd(te), dt(dtime), maxDt(1e100), firstDt(fdt), cFLFactor(
                                                                               1), modulo(mod), timeStep(tist), fileName(name), fixed(true)
    {
    }

private:
    const double tStart;
    const double tEnd;
    double dt;
    const double maxDt;
    const double firstDt;
    const double cFLFactor;
    const int modulo;
    TimeStep<G, Model>& timeStep;
    const char* fileName;
    const bool fixed;
};
}
#endif
