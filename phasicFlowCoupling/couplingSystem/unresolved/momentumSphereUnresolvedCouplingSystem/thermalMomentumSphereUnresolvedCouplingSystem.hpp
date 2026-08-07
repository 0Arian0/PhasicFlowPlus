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

#ifndef pFlow_thermalMomentumSphereUnresolvedCouplingSystem_hpp
#define pFlow_thermalMomentumSphereUnresolvedCouplingSystem_hpp

#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "momentumInteraction.hpp"
#include "heatInteraction.hpp"

namespace pFlow
{
namespace coupling
{

/**
 * @class thermalMomentumSphereUnresolvedCouplingSystem
 * @brief Two-way coupling manager for momentum and heat transfer.
 *
 * @details
 * Handles momentum (drag/lift) and heat transfer (convection, radiation, PFP)
 * coupling between the Eulerian CFD mesh and Lagrangian DEM particles.
 *
 * The protected porosityCoupling() accessor exists specifically to allow
 * derived reactive classes to obtain a valid const porosity& reference
 * during their member-initialiser lists, avoiding the undefined-behaviour
 * static_cast that would otherwise be required on the volScalarField
 * returned by alpha().
 */
class thermalMomentumSphereUnresolvedCouplingSystem
:
    public unresolvedCouplingSystem
{
public:

    //- Type info

        TypeInfo(
            "thermalSphereUnresolvedCouplingSystem<thermalMomentum>");

private:

    //- private members

        // --- Section 1: Interaction Models & Managers ---

        /// @brief Porosity calculator (fluid volume fraction alpha).
        uniquePtr<porosity>         porosity_ = nullptr;

        /// @brief Drag + lift + virtual-mass interaction manager.
        momentumInteraction         momentumInteraction_;

        /// @brief Convective + radiative heat interaction manager.
        heatInteraction             heatInteraction_;

        /// @brief Performance timer for porosity computation.
        Timer                       porosityTimer_;

        /// @brief True when Gaussian / diffusion distribution is required.
        bool                        requiresDistribution_ = false;

        // --- Section 2: MPI Communication Fields (ProcCMFields) ---

        ///< Tp  [K]
        Plus::realProcCMField       particleTemperature_;
        
        ///< Q_conv [W]
        Plus::realProcCMField       fluidHeatSourceConv_;
        
        ///< Q_rad  [W]
        Plus::realProcCMField       fluidHeatSourceRad_;
        
        ///< epsilon  [-]
        Plus::realProcCMField       emissivity_;
        
        ///< Sigma(T) neighbours [K]
        Plus::realProcCMField       radSumTemp_;
        
        ///< N neighbours [-]
        Plus::uint32ProcCMField     radNumPrt_;

        // PFP (Particle-Fluid-Particle) pipeline fields
        Plus::realProcCMField       fluidKappa_;
        
        Plus::realProcCMField       fluidAlpha_;

        // --- Section 2b: Cumulative diagnostic counters ---
        //
        // These accumulate for the ENTIRE run rather than being reset every
        // timestep, so a chronic issue that stays just under the per-step
        // reporting threshold (and would therefore never print a warning on
        // its own) is still visible once its cumulative effect becomes
        // significant. See collectFluidProperties() and
        // sendFluidPropertiesToDEM() for where each is incremented.

        /// @brief Total particle-timesteps, across the whole run, for which
        /// kappa had to be floored to a positive value after MPI collection.
        mutable uint64              cumulativeBadKappaCount_ = 0;

        /// @brief Total particle-timesteps, across the whole run, for which
        /// alpha had to be clamped into [0,1] after MPI collection.
        mutable uint64              cumulativeBadAlphaCount_ = 0;

        /// @brief Total particle-timesteps, across the whole run, for which a
        /// particle's mapped fluid cell index was invalid when sampling
        /// fluidKappa/fluidAlpha for the DEM side PFP model.
        mutable uint64              cumulativeInvalidCellCount_ = 0;

        /// @brief Milestone (in cumulative event count) at which the next
        /// "still occurring" reminder for the bad kappa/alpha counters is
        /// printed, so repeated occurrences are reported periodically rather
        /// than either flooding the log every step or never appearing again
        /// after the very first report.
        mutable uint64              nextKappaAlphaReportMilestone_ = 1;

        /// @brief Same milestone mechanism as above, for the invalid-cell
        /// counter in sendFluidPropertiesToDEM().
        mutable uint64              nextInvalidCellReportMilestone_ = 1;

    //- private methods

        // --- Section 4: Private Helper Methods ---

        bool collectFluidHeatSource();
        
        bool collectFluidProperties();

        void sendFluidHeatSourceToDEM();
        
        void sendFluidPropertiesToDEM();

protected:

    //- protected methods

        // --- Section 3: Protected Interface for Derived Classes ---

        /**
         * @brief Returns a const reference to the internally-owned 
         * porosity object.
         * @return const porosity& - the unique-ptr-owned porosity instance.
         */
        inline
        const porosity& porosityCoupling() const
        {
            return porosity_();
        }

        /**
         * @brief Scatters DEM particle fields from the MPI master to 
         * all workers.
         */
        bool distributeParticleFields() override;

public:

    //- constructors

        // --- Section 5: Constructors / Destructor ---

        thermalMomentumSphereUnresolvedCouplingSystem(
            word            shapeTypeName,
            word            couplingSystemType,
            Foam::fvMesh&   mesh,
            int             argc,
            char*           argv[]);

        thermalMomentumSphereUnresolvedCouplingSystem(
            const thermalMomentumSphereUnresolvedCouplingSystem&) = delete;

        thermalMomentumSphereUnresolvedCouplingSystem& operator=(
            const thermalMomentumSphereUnresolvedCouplingSystem&) = delete;

        ~thermalMomentumSphereUnresolvedCouplingSystem() 
            override = default;

    //- public methods

        add_vCtor(
            unresolvedCouplingSystem,
            thermalMomentumSphereUnresolvedCouplingSystem,
            word
        );

        // --- Section 6: Physical Coupling Calculations ---

        void calculatePorosity() override;
        
        void calculateMomentumCoupling() override;
        
        void calculateHeatCoupling() override;
        
        void calculateMassCoupling() override;

        // --- Section 7: Source-term Accessors ---

        inline
        const Foam::volScalarField& alpha() const override
        {
            return porosity_().alpha();
        }

        Foam::tmp<Foam::volScalarField> Sp() const override;
        
        Foam::tmp<Foam::volVectorField> Su() const override;
        
        Foam::tmp<Foam::volScalarField> heatSource() const override;

        /// @brief Implicit coefficient for the fluid energy equation 
        /// matrix [W/(m^3.K)].
        inline
        Foam::tmp<Foam::volScalarField> heatSp() const
        {
            return Foam::tmp<Foam::volScalarField>(
                heatInteraction_.heatSp());
        }

        /// @brief Explicit source term for the fluid energy equation 
        /// matrix [W/m^3].
        inline
        Foam::tmp<Foam::volScalarField> heatSu() const
        {
            return Foam::tmp<Foam::volScalarField>(
                heatInteraction_.heatSu());
        }

        /// @brief Solid velocity field mapped to the Eulerian mesh.
        Foam::tmp<Foam::volVectorField> Us() const override;

        /// @brief Fallback mass source method required by base class.
        inline
        Foam::tmp<Foam::volScalarField> massSource(
            const word& /*specieName*/) const override
        {
            notImplementedFunction;
            return Foam::tmp<Foam::volScalarField>(nullptr);
        }

        // --- Section 8: Identity ---

        inline
        word shapeTypeName() const override
        {
            return "sphere";
        }
        
        inline
        word couplingSystemType() const override
        {
            return "thermalMomentum";
        }
        
        inline
        bool requireCellDistribution() const override
        {
            return requiresDistribution_;
        }

        // --- Section 9: Data Synchronisation ---

        bool sendDataToDEM(real t, real dt) override;

        inline
        Plus::realProcCMField& particleTemperature()
        {
            return particleTemperature_;
        }
        
        inline
        Plus::realProcCMField& fluidHeatSourceConv()
        {
            return fluidHeatSourceConv_;
        }
        
        inline
        Plus::realProcCMField& fluidHeatSourceRad()
        {
            return fluidHeatSourceRad_;
        }
        
        inline
        Plus::realProcCMField& fluidKappa()
        {
            return fluidKappa_;
        }
        
        inline
        Plus::realProcCMField& fluidAlpha()
        {
            return fluidAlpha_;
        }

        // --- Section 10: Diagnostic Accessors ---

        /// @brief Total particle-timesteps for which kappa or alpha needed
        /// clamping after MPI collection, accumulated over the whole run.
        inline
        uint64 cumulativeBadKappaCount() const
        {
            return cumulativeBadKappaCount_;
        }
        
        inline
        uint64 cumulativeBadAlphaCount() const
        {
            return cumulativeBadAlphaCount_;
        }

        /// @brief Total particle-timesteps for which the mapped fluid cell was
        /// invalid while sampling fluidKappa/fluidAlpha for the DEM side,
        /// accumulated over the whole run.
        inline
        uint64 cumulativeInvalidCellCount() const
        {
            return cumulativeInvalidCellCount_;
        }

}; // thermalMomentumSphereUnresolvedCouplingSystem

} // coupling
} // pFlow

#endif // pFlow_thermalMomentumSphereUnresolvedCouplingSystem_hpp



