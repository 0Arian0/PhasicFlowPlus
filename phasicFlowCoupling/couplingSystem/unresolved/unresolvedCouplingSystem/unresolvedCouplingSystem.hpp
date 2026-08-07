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

#ifndef pFlow_unresolvedCouplingSystem_hpp
#define pFlow_unresolvedCouplingSystem_hpp

#include "couplingSystem.hpp"
#include "distributionBase.hpp"
#include "virtualConstructor.hpp"

namespace pFlow
{
namespace coupling
{

/**
 * @class unresolvedCouplingSystem
 * @ingroup couplingSystem
 * @brief Base class for Eulerian-Lagrangian coupling systems in unresolved 
 *        particle-fluid flows.
 *
 * The `unresolvedCouplingSystem` class serves as the abstract base for 
 * implementing unresolved coupling between particles and fluid in 
 * multiphase flows. 
 *
 * ## Key Responsibilities
 *
 * - **Particle Distribution**: Manages how particle-scale properties (velocity,
 *   forces, etc.) are mapped to fluid cells using configurable distribution 
 *   methods (e.g., Gaussian, diffusion, PCM, sub-division methods).
 * - **Mesh Management**: Handles the coupling mesh that bridges particles and 
 *   fluid cells.
 * - **Abstract Interface**: Defines the interface for derived classes to 
 *   implement specific coupling physics (momentum, heat, mass transfer).
 *
 * ## Architecture
 *
 * This class is part of a three-level hierarchy:
 * - **unresolvedCouplingSystem** (base): Defines distribution and mesh 
 *   management.
 * - **Derived Classes**: (e.g., momentumSphereUnresolvedCouplingSystem).
 *   Implement specific physical coupling mechanisms.
 *
 * ## Distribution Methods
 *
 * The distribution method determines how particle data is mapped to cells:
 * - **PCM**: Direct assignment to containing cell.
 * - **Gaussian**: Gaussian kernel-based distribution.
 * - **GaussianIntegral**: Volume-weighted Gaussian integration.
 * - **adaptiveGaussian**: Adaptive kernel based on cell-to-particle size ratio.
 * - **diffusion**: Laplacian diffusion-based smoothing.
 * - **subDivision9/29**: Geometric sub-division methods.
 *
 * @see couplingSystem, momentumSphereUnresolvedCouplingSystem
 * @see distributionBase, porosity, momentumInteraction
 */
class unresolvedCouplingSystem
:
    public couplingSystem
{
public:

    //- Type info

        /// Virtual constructor macro for registering derived classes.
        create_vCtor(
            unresolvedCouplingSystem,
            word,
            (
                word            shapeTypeName,
                word            couplingSystemType, 
                Foam::fvMesh&   mesh,
                int             argc, 
                char*           argv[]
            ),
            (shapeTypeName, couplingSystemType, mesh, argc, argv)
        );

private:

    //- private members

        /// Distribution method for mapping particle properties to fluid cells.
        uniquePtr<distributionBase> distribution_;

public:

    //- constructors

        /// Constructor initializing the unresolved coupling system.
        unresolvedCouplingSystem(
            word            shapeTypeName,
            word            couplingSystemType, 
            Foam::fvMesh&   mesh,
            int             argc, 
            char*           argv[]);

        /// Deleted copy constructor to prevent unintended copying.
        unresolvedCouplingSystem(
            const unresolvedCouplingSystem&) = delete;

        /// Deleted copy assignment operator to prevent unintended copying.
        unresolvedCouplingSystem& operator=(
            const unresolvedCouplingSystem&) = delete;

        /// Deleted move constructor to prevent unintended moving.
        unresolvedCouplingSystem(
            unresolvedCouplingSystem&&) = delete;

        /// Deleted move assignment operator to prevent unintended moving.
        unresolvedCouplingSystem& operator=(
            unresolvedCouplingSystem&&) = delete;

        /// Virtual destructor for proper cleanup of derived classes.
        ~unresolvedCouplingSystem() override = default;

    //- public methods

        /// Access the `unresolved` coupling configuration subdictionary.
        const Foam::dictionary& unresolvedDict() const;

        /// Const reference to the distribution method.
        inline
        const distributionBase& distribution() const
        {
            return *distribution_;
        }

        /// Get the name of the distribution method.
        inline
        const Foam::word distributionMethodName() const
        {
            return distribution_->distributionMethodName();
        }

        /// Const reference to the cell indices where particles are located.
        inline
        const Plus::procCMField<Foam::label>& parCellIndex() const
        {
            return cMesh().parCellIndex();
        }

        /// Update distribution weights based on current particle positions.
        inline
        void updateDistributionWeights()
        {
            if (distribution_)
            {
                distribution_->updateWeights(this->particleDiameter());
            }
        }

        /// Calculate local fluid volume fraction (porosity) in cells.
        virtual void calculatePorosity() = 0;

        /// Calculate momentum coupling between particles and fluid.
        virtual void calculateMomentumCoupling() = 0;

        /// Calculate heat coupling between particles and fluid.
        virtual void calculateHeatCoupling() = 0;

        /// Calculate mass coupling between particles and fluid.
        virtual void calculateMassCoupling() = 0;

        /// Returns the fluid volume fraction field (porosity).
        virtual const Foam::volScalarField& alpha() const = 0;

        /// Returns the implicit coefficient of momentum source term (Sp).
        virtual Foam::tmp<Foam::volScalarField> Sp() const = 0;

        /// Returns the explicit part of momentum source term (Su).
        virtual Foam::tmp<Foam::volVectorField> Su() const = 0;

        /// Returns the heat source term for the energy equation.
        virtual Foam::tmp<Foam::volScalarField> heatSource() const = 0;

        /// Returns the mass source term for the species transport equation.
        virtual Foam::tmp<Foam::volScalarField> massSource(
            const word& specieName) const = 0;

        /// Returns the particle shape type name (e.g., sphere, grain).
        virtual word shapeTypeName() const = 0;

        /// Returns the coupling system type (e.g., momentum, heatMomentum).
        virtual word couplingSystemType() const = 0;

        /// Returns the velocity of the injected mass source.
        virtual Foam::tmp<Foam::volVectorField> Us() const;
        
        /// Returns explicit reaction heat source (Fallback: Zero)
        virtual Foam::tmp<Foam::volScalarField> rxnHeatSu() const
        {
            return Foam::tmp<Foam::volScalarField>::New(
                Foam::IOobject(
                    "dummyRxnHeat", 
                    this->cMesh().mesh().time().timeName(), 
                    this->cMesh().mesh(), 
                    Foam::IOobject::NO_READ, 
                    Foam::IOobject::NO_WRITE, 
                    false),
                this->cMesh().mesh(),
                Foam::dimensionedScalar(
                    "zero", 
                    Foam::dimEnergy/Foam::dimVolume/Foam::dimTime, 
                    0.0));
        }

        /// Returns species mass sources (Fallback: Empty List)
        virtual const Foam::PtrList<Foam::volScalarField>& speciesSu() const
        {
            static Foam::PtrList<Foam::volScalarField> dummyList;
            return dummyList;
        }

        /// Returns species implicit sink (Fallback: Empty List)
        virtual const Foam::PtrList<Foam::volScalarField>& speciesSp() const
        {
            static Foam::PtrList<Foam::volScalarField> dummyList;
            return dummyList;
        }

        /// Indicates if the coupling requires cell-level distribution.
        virtual bool requireCellDistribution() const = 0;

        /// Static factory method to create derived coupling system instances.
        static uniquePtr<unresolvedCouplingSystem> create(
            word            shapeTypeName,
            word            couplingSystemType, 
            Foam::fvMesh&   mesh,
            int             argc, 
            char*           argv[]);

}; 

} // coupling
} // pFlow

#endif // pFlow_unresolvedCouplingSystem_hpp


