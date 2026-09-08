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

#ifndef pFlow_heatInteraction_hpp
#define pFlow_heatInteraction_hpp

#include "OFCompatibleHeader.hpp"
#include "virtualConstructor.hpp"
#include "Timer.hpp"
#include "procCMFields.hpp"
#include "heatTransfer.hpp"
#include "PCM.hpp"

namespace pFlow
{
namespace coupling
{

// Forward declarations
class unresolvedCouplingSystem;
class porosity;
class fluidAveraging;
class solidAveraging;

/**
 * @class heatInteraction
 * @brief Manages the high-level mapping and execution of fluid-particle 
 *        heat transfer.
 *
 * @details
 * Acts as the interface between the coupling system and the underlying 
 * thermodynamic closure models. It handles the spatial distribution mapping 
 * between the discrete Lagrangian particles and the continuous Eulerian mesh.
 *
 * Mapping Strategies Supported:
 * - **Cell (PCM)**: Particle-in-Cell method. The entire heat source of a 
 *   particle is assigned exclusively to the single fluid cell containing 
 *   its center.
 * - **Distribution**: Advanced volumetric mapping (e.g., Gaussian) where a 
 *   particle's heat source is smoothly distributed across neighboring cells.
 */
class heatInteraction
{
public:

    //- Type info

        TypeInfo("heatInteraction");

private:
    
    //- private members

        const porosity&             porosity_;
        
        bool                        heatExchangeDistribute_;
        
        bool                        requireCellDistribution_ = false;

        uniquePtr<heatTransfer>     heatTransfer_;
        
        uniquePtr<PCM>              noDistribution_ = nullptr;
        
        Timer                       heatInteractionTimer_;

public:
    
    //- constructors

        heatInteraction(
            const unresolvedCouplingSystem& uCS, 
            const porosity&                 prsty);
        
        virtual ~heatInteraction() = default;

    //- public methods

        /// @brief Implicit source coefficient (Sp) for the fluid energy 
        /// equation matrix.
        inline
        const Foam::volScalarField& heatSp() const 
        { 
            return std::as_const<const heatTransfer&>(*heatTransfer_).Sp(); 
        }

        /// @brief Explicit source vector (Su) for the fluid energy equation 
        /// matrix.
        inline
        const Foam::volScalarField& heatSu() const 
        { 
            return std::as_const<const heatTransfer&>(*heatTransfer_).Su(); 
        }

        /// @brief Alias for heatSu, returning the total volumetric heat source.
        inline
        const Foam::volScalarField& Sh() const 
        { 
            return heatSu(); 
        }

        inline
        const porosity& Porosity() const
        {
            return porosity_;
        }

        const unresolvedCouplingSystem& uCS() const;

        inline
        bool requireCellDistribution() const 
        { 
            return requireCellDistribution_; 
        }

        const Foam::dictionary& dict() const;

        static const Foam::dictionary& getDict(
            const unresolvedCouplingSystem& uCS);

        /**
         * @brief Calculates convective and radiative heat exchanges.
         * @details Gathers required physical fields and delegates the actual 
         * calculation to the instantiated `heatTransfer` model, applying the 
         * user-selected distribution mapping.
         *
         * @param fluidVelocity Interpolated fluid velocity at particle centers.
         * @param parVelocity   Discrete particle velocity.
         * @param dp            Particle diameter.
         * @param Tp            Particle temperature.
         * @param Qp            [OUT] Computed convective heat source.
         * @param emissivity    Particle emissivity.
         * @param radSumTemp    Sum of neighboring temperatures for radiation.
         * @param radNumPrt     Number of neighbors considered in radiation.
         * @param QpRad         [OUT] Computed radiative heat source.
         */
        virtual void calculateCoupling(
            const fluidAveraging&           fluidVelocity,
            const solidAveraging&           parVelocity,
            const Plus::realProcCMField&    dp,
            const Plus::realProcCMField&    Tp,
            Plus::realProcCMField&          Qp,
            const Plus::realProcCMField&    emissivity,
            const Plus::realProcCMField&    radSumTemp,
            const Plus::uint32ProcCMField&  radNumPrt,
            Plus::realProcCMField&          QpRad);

}; // heatInteraction

} // coupling
} // pFlow

#endif // pFlow_heatInteraction_hpp
