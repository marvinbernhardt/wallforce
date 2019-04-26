/*
 * This file is part of the GROMACS molecular simulation package.
 *
 * Copyright (c) 2011,2012,2013,2014,2015, by the GROMACS development team, led by
 * Mark Abraham, David van der Spoel, Berk Hess, and Erik Lindahl,
 * and including many others, as listed in the AUTHORS file in the
 * top-level source directory and at http://www.gromacs.org.
 *
 * GROMACS is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1
 * of the License, or (at your option) any later version.
 *
 * GROMACS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with GROMACS; if not, see
 * http://www.gnu.org/licenses, or write to the Free Software Foundation,
 * Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
 *
 * If you want to redistribute modifications to GROMACS, please
 * consider that scientific software is very special. Version
 * control is crucial - bugs must be traceable. We will be happy to
 * consider code for inclusion in the official distribution, but
 * derived work must not be called official GROMACS. Details are found
 * in the README & COPYING files - if they are missing, get the
 * official version at http://www.gromacs.org.
 *
 * To help us fund GROMACS development, we humbly ask that you cite
 * the research papers on the package. Check out http://www.gromacs.org.
 */
#include <string>
#include <vector>

#include <gromacs/trajectoryanalysis.h>

using namespace gmx;


enum AxisType
{
    eAxisType_x,
    eAxisType_y,
    eAxisType_z
};
//! Strings corresponding to AxisType
const char *const c_axisTypes[] = { "x", "y", "z"};


class WallForceCalculator : public TrajectoryAnalysisModule
{
    public:
        WallForceCalculator();

        virtual void initOptions(IOptionsContainer          *options,
                TrajectoryAnalysisSettings *settings);
        virtual void initAnalysis(const TrajectoryAnalysisSettings &settings,
                const TopologyInformation        &top);

        virtual void analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
                TrajectoryAnalysisModuleData *pdata);

        virtual void finishAnalysis(int nframes);
        virtual void writeOutput();

    private:
        class ModuleData;

        std::string                      fnDist_;
        Selection                        sel_;
        real                             wall_pos_;
        AxisType                         wall_axis_;
        real                             wall_k_;

        AnalysisData                     data_;
        AnalysisDataAverageModulePointer avem_;
};


WallForceCalculator::WallForceCalculator()
    : wall_pos_(0.0), wall_axis_(eAxisType_z), wall_k_(0.0)
{
    registerAnalysisDataset(&data_, "avedist");
}


    void
WallForceCalculator::initOptions(IOptionsContainer          *options,
        TrajectoryAnalysisSettings *settings)
{
    static const char *const desc[] = {
        "Calculate the forces of wall restraints.[PAR]",
        "You can use selections."
    };

    settings->setHelpText(desc);

    options->addOption(FileNameOption("o")
            .filetype(eftPlot).outputFile()
            .store(&fnDist_).defaultBasename("avedist")
            .description("Average distances from reference group"));

    options->addOption(SelectionOption("select")
            .store(&sel_).required()
            .description("Group to calculate wall force from"));
    options->addOption(RealOption("wallr").store(&wall_pos_)
                           .description("Wall position as distance from the origin"));

    options->addOption(EnumOption<AxisType>("axis").store(&wall_axis_)
                           .enumValue(c_axisTypes)
                           .description("Axis normal to the wall"));

    options->addOption(RealOption("wallk").store(&wall_k_)
                           .description("Wall force constant"));

    settings->setFlag(TrajectoryAnalysisSettings::efRequireTop);
}


    void
WallForceCalculator::initAnalysis(const TrajectoryAnalysisSettings &settings,
        const TopologyInformation         & /*top*/)
{

    data_.setColumnCount(0, 1);

    avem_.reset(new AnalysisDataAverageModule());
    data_.addModule(avem_);

    if (!fnDist_.empty())
    {
        AnalysisDataPlotModulePointer plotm(
                new AnalysisDataPlotModule(settings.plotSettings()));
        plotm->setFileName(fnDist_);
        plotm->setTitle("Average forece");
        plotm->setXAxisIsTime();
        plotm->setYLabel("Force (kJ/mol/nm)");
        data_.addModule(plotm);
    }
}


    void
WallForceCalculator::analyzeFrame(int frnr, const t_trxframe &fr, t_pbc *pbc,
        TrajectoryAnalysisModuleData *pdata)
{
    AnalysisDataHandle         dh     = pdata->dataHandle(data_);
    const Selection           &sel = pdata->parallelSelection(sel_);

    dh.startFrame(frnr, fr.time);
    int nr = sel.posCount();
    real ftot = 0.0;
    real k_wall = 0.1;

    int xyz;
    if (wall_axis_ == eAxisType_x) xyz = 0;
    else if (wall_axis_ == eAxisType_y) xyz = 1;
    else if (wall_axis_ == eAxisType_z) xyz = 2;
    else throw std::invalid_argument("unknown axis");

    for (int i = 0; i < nr; ++i)
    {
        SelectionPosition p = sel.position(i);
        if (p.x()[xyz] > wall_pos_)
        {
            ftot += -k_wall * (p.x()[xyz] - wall_pos_);
        }
    }
    dh.setPoint(0, ftot);
    dh.finishFrame();
}


    void
WallForceCalculator::finishAnalysis(int nframes)
{
}


    void
WallForceCalculator::writeOutput()
{
    // We print out the average of the mean distances for each group.
    fprintf(stderr, "Average force for '%s': %.3f kJ/mol/nm\n",
            sel_.name(), avem_->average(0, 0));
}

/*! \brief
 * The main function for the analysis template.
 */
    int
main(int argc, char *argv[])
{
    return gmx::TrajectoryAnalysisCommandLineRunner::runAsMain<WallForceCalculator>(argc, argv);
}
