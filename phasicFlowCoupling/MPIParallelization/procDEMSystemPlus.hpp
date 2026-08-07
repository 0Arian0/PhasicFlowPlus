/*------------------------------- phasicFlow ---------------------------------
      O        C enter of
     O O       E ngineering and
    O   O      M ultiscale modeling of
   OOOOOOO     F luid flow
------------------------------------------------------------------------------
  Copyright (C): www.cemf.ir
  email: hamid.r.norouzi AT gmail.com
------------------------------------------------------------------------------
Licence:
  This file is part of phasicFlow code. It is a free software for simulating
  granular and multiphase flows. You can redistribute it and/or modify it under
  the terms of GNU General Public License v3 or any other later versions.

  phasicFlow is distributed to help others in their research in the field of
  granular and multiphase flows, but WITHOUT ANY WARRANTY; without even the
  implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

-----------------------------------------------------------------------------*/

#ifndef pFlow_procDEMSystemPlus_hpp
#define pFlow_procDEMSystemPlus_hpp

#include "DEMSystem.hpp"
#include "processorPlus.hpp"
#include "procVectorPlus.hpp"
#include <vector>
#include <string>

namespace pFlow
{
namespace Plus
{

/**
 * @class procDEMSystem
 * @brief Manages top-level DEM execution and parallel MPI memory routing.
 *
 * Serves as the principal interface between OpenFOAM processors and the
 * underlying Kokkos GPU/CPU particle system data streams.
 */
class procDEMSystem
{
protected:

    //- protected members

        /// Managed pointer to the underlying DEM system (valid only on master).
        uniquePtr<DEMSystem>    demSystem_ = nullptr;

        /// Simulation start time synchronised across all MPI ranks.
        real                    startTime_ = 0;

public:

    //- constructors

        procDEMSystem(
            word    demSystemName,
            int     argc,
            char*   argv[],
            bool    requireRVel = false);

        virtual ~procDEMSystem() = default;

    //- public methods

        // ================================================================= //
        // Global Time & Control
        // ================================================================= //

        inline
        real startTime() const
        {
            return startTime_;
        }

        inline
        bool getDataFromDEM()
        {
            if (demSystem_) return demSystem_->beforeIteration();
            return true;
        }

        inline
        bool updateParticleDistribution(
            real                    extentFraction,
            const procVector<box>&  domains)
        {
            if (demSystem_)
            {
                return demSystem_->updateParticleDistribution(
                    extentFraction, 
                    domains);
            }
            return true;
        }

        // ================================================================= //
        // Mechanical Field Accessors
        // ================================================================= //

        inline
        span<const int32> parIndexInDomain(int32 di) const
        {
            if (demSystem_) return demSystem_->parIndexInDomain(di);
            return span<const int32>();
        }

        inline
        span<uint32> particleIdAllMaster() const
        {
            if (demSystem_) return demSystem_->particleId();
            return span<uint32>();
        }

        inline
        span<realx3> particlesCenterMassAllMaster()
        {
            if (demSystem_) return demSystem_->position();
            return span<realx3>();
        }

        inline
        span<realx3> particlesVelocityAllMaster()
        {
            if (demSystem_) return demSystem_->velocity();
            return span<realx3>();
        }

        inline
        span<realx3> particlesRVelocityAllMaster()
        {
            if (demSystem_) return demSystem_->rVelocity();
            return span<realx3>();
        }

        inline
        span<realx3> particlesFluidForceAllMaster()
        {
            if (demSystem_) return demSystem_->parFluidForce();
            return span<realx3>();
        }

        inline
        span<realx3> particlesAccelerationAllMaster()
        {
            if (demSystem_) return demSystem_->acceleration();
            return span<realx3>();
        }

        inline
        span<realx3> particlesFluidTorqueAllMaster()
        {
            if (demSystem_) return demSystem_->parFluidTorque();
            return span<realx3>();
        }

        inline
        std::vector<real> shapeDiametersAllMaster() const
        {
            if (demSystem_) return demSystem_->shapeDiameters();
            return std::vector<real>{};
        }

        inline
        span<real> particlesDiameterAllMaster()
        {
            if (demSystem_) return demSystem_->diameter();
            return span<real>();
        }

        inline
        span<real> particlesCourseGrainFactorMasterAllMaster()
        {
            if (demSystem_) return demSystem_->courseGrainFactor();
            return span<real>();
        }

        inline
        procVector<int32> numParInDomainMaster() const
        {
            if (demSystem_) return demSystem_->numParInDomains();
            return procVector<int32>(true);
        }

        inline
        procVector<span<const int32>> parIndexInDomainsMaster() const
        {
            procVector<span<const int32>> parIndex(true);
            if (demSystem_)
            {
                for (size_t i = 0; i < parIndex.size(); ++i)
                {
                    parIndex[i] = demSystem_->parIndexInDomain(i);
                }
            }
            return parIndex;
        }

        // ================================================================= //
        // Thermal Coupling Extensions
        // ================================================================= //

        inline
        span<real> temperature()
        {
            if (demSystem_) return demSystem_->temperature();
            return span<real>();
        }

        inline
        span<real> emissivity()
        {
            if (demSystem_) return demSystem_->emissivity();
            return span<real>();
        }

        inline
        span<real> radSumTemp()
        {
            if (demSystem_) return demSystem_->radSumTemp();
            return span<real>();
        }

        inline
        span<uint32> radNumPrt()
        {
            if (demSystem_) return demSystem_->radNumPrt();
            return span<uint32>();
        }

        inline
        span<real> parFluidHeatSourceConv()
        {
            if (demSystem_) return demSystem_->parFluidHeatSourceConv();
            return span<real>();
        }

        inline
        span<real> parFluidHeatSourceRad()
        {
            if (demSystem_) return demSystem_->parFluidHeatSourceRad();
            return span<real>();
        }

        inline
        span<real> parFluidKappa()
        {
            if (demSystem_) return demSystem_->parFluidKappa();
            return span<real>();
        }

        inline
        span<real> parFluidAlpha()
        {
            if (demSystem_) return demSystem_->parFluidAlpha();
            return span<real>();
        }

        // ================================================================= //
        // Chemical Reaction Extensions
        // ================================================================= //

        inline
        span<real> solidMassFractions()
        {
            if (demSystem_) return demSystem_->solidMassFractions();
            return span<real>();
        }

        inline
        span<real> gasMassSource()
        {
            if (demSystem_) return demSystem_->gasMassSource();
            return span<real>();
        }

        inline
        span<real> gasMassSourceSp()
        {
            if (demSystem_) return demSystem_->gasMassSourceSp();
            return span<real>();
        }

        inline
        span<real> gasConcentrations()
        {
            if (demSystem_) return demSystem_->gasConcentrations();
            return span<real>();
        }

        /// @brief Per-particle solid-side reaction heat: (1-n)*Q_rxn [W].
        inline
        span<real> reactionHeat()
        {
            if (demSystem_) return demSystem_->reactionHeat();
            return span<real>();
        }

        /// @brief Per-particle fluid-side reaction heat: n*Q_rxn [W].
        ///        Zero-filled when n = 0 (default for surface reactions).
        inline
        span<real> reactionHeatFluid()
        {
            if (demSystem_) return demSystem_->reactionHeatFluid();
            return span<real>();
        }

        /**
         * @brief Returns gas species names from the DEM kinetics order.
         *        Used by buildSpeciesMapping() cross-check at CFD startup.
         */
        inline
        std::vector<std::string> gasSpeciesNames() const
        {
            if (demSystem_) return demSystem_->gasSpeciesNames();
            return std::vector<std::string>();
        }

        /**
         * @brief Returns gas species molar masses [kg/mol] from the DEM
         *        kinetics, in the same order as gasSpeciesNames().
         *        Used by buildSpeciesMapping() cross-check at CFD startup to
         *        catch a stale/mismatched transportProperties/gasMw entry.
         */
        inline
        std::vector<real> gasMolarMasses() const
        {
            if (demSystem_) return demSystem_->gasMolarMasses();
            return std::vector<real>();
        }

        // NOTE (R9): sendReactionDataToDEM() wrapper removed.
        //   The base DEMSystem::sendReactionDataToDEM() was never called from
        //   the CFD coupling code and returned false (indicating failure).
        //   Removing dead code prevents accidental calls that silently
        //   indicate failure without any visible effect.
        //   Reaction data flow:
        //     DEM->CFD: reactionDataHostUpdatedSync() in getDataFromDEM()
        //     CFD->DEM: sendGasConcentrationsToDEM() (below)

        // ================================================================= //
        // Sync Dispatchers — Host -> DEM Device
        // ================================================================= //

        inline
        bool sendFluidForceToDEM()
        {
            if (demSystem_) return demSystem_->sendFluidForceToDEM();
            return true;
        }

        inline
        bool sendFluidTorqueToDEM()
        {
            if (demSystem_) return demSystem_->sendFluidTorqueToDEM();
            return true;
        }

        inline
        bool sendFluidHeatSourcesToDEM()
        {
            if (demSystem_) return demSystem_->sendFluidHeatSourcesToDEM();
            return true;
        }

        inline
        bool sendFluidPropertiesToDEM()
        {
            if (demSystem_) return demSystem_->sendFluidPropertiesToDEM();
            return true;
        }

        inline
        bool sendGasConcentrationsToDEM()
        {
            if (demSystem_) return demSystem_->sendGasConcentrationsToDEM();
            return true;
        }

        // ================================================================= //
        // Timestep Execution
        // ================================================================= //

        inline
        bool iterate(real upToTime, bool writeTime, const word& timeName)
        {
            if (demSystem_)
            {
                if (writeTime)
                {
                    return demSystem_->iterate(upToTime, upToTime, timeName);
                }
                else
                {
                    return demSystem_->iterate(upToTime);
                }
            }
            return true;
        }

        inline
        bool iterate(real upToTime)
        {
            if (demSystem_) return demSystem_->iterate(upToTime);
            return true;
        }

        inline
        Timers* getTimers()
        {
            if (demSystem_) return &demSystem_->Control().timers();
            return nullptr;
        }

        // Legacy backward-compatibility alias
        inline
        span<real> particlesTemperatureAllMaster()
        {
            return temperature();
        }

}; // procDEMSystem

} // Plus
} // pFlow

#endif // pFlow_procDEMSystemPlus_hpp


