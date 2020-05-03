// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Particles
 * \brief A simple particle vtk writer
 */

#ifndef DUMUX_PARTICLES_PARTICLE_VTK_WRITER_HH
#define DUMUX_PARTICLES_PARTICLE_VTK_WRITER_HH

#include <memory>
#include <fstream>
#include <iomanip>
#include <vector>

namespace Dumux {

/*!
 * \ingroup Particles
 * \brief A simple particle vtk writer
 */
template <class ParticleCloud>
class ParticleVTKWriter
{
    static constexpr std::size_t numBeforeLineBreak = 5;
public:

    ParticleVTKWriter(std::shared_ptr<ParticleCloud> cloud)
    : cloud_(cloud)
    {}

    void write(const std::string& name) const
    {
        std::ofstream vtpFile(name + ".vtp");
        writeHeader_(vtpFile, cloud_->size());

        // points
        writeCoordinates_(vtpFile);

        // particle data
        vtpFile << "<PointData Scalars=\"id\">\n";
        // always write out particle id to enable tracking paraview
        writeParticleId_(vtpFile, "id");
        for (std::size_t i = 0; i < particleData_.size(); ++i)
            writeParticleData_(vtpFile, names_[i], *(particleData_[i]));
        vtpFile << "</PointData>\n";

        writeFooter_(vtpFile);
    }

    //! add a particle data vector
    void addParticleData(const std::vector<double>& particleData, const std::string& name)
    {
        names_.push_back(name);
        particleData_.push_back(&particleData);
    }

private:
    std::shared_ptr<const ParticleCloud> cloud_;
    std::vector<std::string> names_;
    std::vector<const std::vector<double>*> particleData_;

    /*!
     * \brief Writes the header to an output stream
     */
    void writeHeader_(std::ostream& file, std::size_t numParticles) const
    {
        std::string header = "<?xml version=\"1.0\"?>\n";
                    header += "<VTKFile type=\"PolyData\" version=\"0.1\" byte_order=\"LittleEndian\">\n";
                    header += "<PolyData>\n";
                    header += "<Piece NumberOfLines=\"" + std::to_string(0) + "\" NumberOfPoints=\"" + std::to_string(numParticles) + "\">\n";
        file << header;
    }

    /*!
     * \brief Writes the coordinates to the file
     * \param file The output file
     */
    void writeCoordinates_(std::ostream& file) const
    {
        // write the positions to the file
        file << "<Points>\n";
        file << "<DataArray type=\"Float32\" Name=\"Coordinates\" NumberOfComponents=\"3\" format=\"ascii\">\n";
        std::size_t counter = 0;
        for (const auto& particle : particles(*cloud_))
        {
            const auto& pos = particle.position();
            file << pos << " ";
            // add missing coordinates if dimension is < 3
            for (int dimIdx = pos.size(); dimIdx < 3; ++dimIdx)
                file << "0.0 ";

            // introduce a line break after a certain time
            if (++counter; counter > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
        file << "</Points>\n";
    }

    void writeParticleData_(std::ostream& file, const std::string& name, const std::vector<double>& data) const
    {
        file << "<DataArray type=\"Float32\" Name=\"" << name << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        std::size_t counter = 0;
        for (const auto& particle : particles(*cloud_))
        {
            file << data[particle.id()] << " ";

            // introduce a line break after a certain time
            if(++counter; counter > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    void writeParticleId_(std::ostream& file, const std::string& name) const
    {
        file << "<DataArray type=\"Int64\" Name=\"" << name << "\" NumberOfComponents=\"1\" format=\"ascii\">\n";
        std::size_t counter = 0;
        for (const auto& particle : particles(*cloud_))
        {
            file << particle.id() << " ";

            // introduce a line break after a certain time
            if(++counter; counter > numBeforeLineBreak)
            {
                file << '\n';
                counter = 0;
            }
        }
        file << "\n</DataArray>\n";
    }

    /*!
     * \brief Writes the footer to the file
     */
    void writeFooter_(std::ostream& file) const
    {
        file << "</Piece>\n";
        file << "</PolyData>\n";
        file << "</VTKFile>";
    }
};

/**
 * \file
 * \ingroup Particles
 * \brief A pvd writer (vtk time series) for particles
 */
template <class ParticleCloud>
class ParticlePVDWriter
{
public:
    ParticlePVDWriter(std::shared_ptr<ParticleCloud> cloud, const std::string& name)
    : vtkWriter_(cloud)
    , name_(name)
    {}

    //! write a time point
    void write(double time)
    {
        const auto count = timeSteps_.size();
        timeSteps_.push_back(time);

        // write vtp file
        vtkWriter_.write(seqName_(count));

        // write pvd file
        std::ofstream pvdFile;
        pvdFile.exceptions(std::ios_base::badbit | std::ios_base::failbit |
                           std::ios_base::eofbit);
        std::string pvdname = name_ + ".pvd";
        pvdFile.open(pvdname.c_str());
        pvdFile << "<?xml version=\"1.0\"?> \n"
                << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\"> \n"
                << "<Collection> \n";

        for (std::size_t i = 0; i <= count; i++)
        {
            pvdFile << "<DataSet timestep=\"" << timeSteps_[i]
                    << "\" group=\"\" part=\"0\" name=\"\" file=\""
                    << seqName_(i) << ".vtp" << "\"/> \n";
        }

        pvdFile << "</Collection> \n"
                << "</VTKFile> \n" << std::flush;
    }

    //! add particle data vector
    void addParticleData(const std::vector<double>& particleData, const std::string& name)
    { vtkWriter_.addParticleData(particleData, name); }

private:
    std::string seqName_(std::size_t count) const
    {
        std::stringstream n;
        n << name_ << "-" << std::setfill('0') << std::setw(5) << count;
        return n.str();
    }

    ParticleVTKWriter<ParticleCloud> vtkWriter_;
    const std::string name_;
    std::vector<double> timeSteps_;
};

} // end namespace Dumux

#endif
