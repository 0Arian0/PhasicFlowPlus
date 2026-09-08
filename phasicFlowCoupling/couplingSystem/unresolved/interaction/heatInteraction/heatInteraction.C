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

#include "heatInteraction.hpp"
#include "unresolvedCouplingSystem.hpp"
#include "porosity.hpp"
#include "fluidAveraging.hpp"
#include "solidAveraging.hpp"

namespace pFlow
{
namespace coupling
{

//----------------------------- constructors ----------------------------------

heatInteraction::heatInteraction(
    const unresolvedCouplingSystem& uCS, 
    const porosity&                 prsty)
:
    porosity_(prsty),
    heatInteractionTimer_(
        "heatInteraction", 
        &uCS.couplingTimers())
{
    // Extract the user's chosen mapping method from the simulation dictionary
    auto heatExch = dict().get<Foam::word>("heatExchange");

    if (heatExch == "distribution")
    {
        // Advanced mapping: Heat source is spread across multiple fluid cells
        heatExchangeDistribute_ = true;
    }
    else if (heatExch == "cell")
    {
        // Point-Center mapping: Heat source is dumped into a single fluid cell
        heatExchangeDistribute_ = false;
    }
    else
    {
        Foam::Info << "Unknown heat exchange method: " << heatExch 
                   << " in " << dict().name() << Foam::endl;
        Plus::processor::abort(0);
    }

    // Instantiate the physical calculation model (delegation)
    heatTransfer_ = heatTransfer::create(uCS, porosity_);
    requireCellDistribution_ = heatExchangeDistribute_;

    // heatExchangeDistribute_ == true (dictionary value "distribution")
    // is the branch that constructs and uses a PCM mapper below.
    if (heatExchangeDistribute_)
    {
        noDistribution_ = makeUnique<PCM>(uCS.cMesh(), uCS.centerMass());
    }
}

//---------------------------- public methods ---------------------------------

const unresolvedCouplingSystem& heatInteraction::uCS() const 
{ 
    return porosity_.uCS(); 
}

const Foam::dictionary& heatInteraction::dict() const 
{ 
    return heatInteraction::getDict(uCS()); 
}

const Foam::dictionary& heatInteraction::getDict(
    const unresolvedCouplingSystem& uCS)
{
    return uCS.unresolvedDict().subDict("heatInteraction");
}

void heatInteraction::calculateCoupling(
    const fluidAveraging&           fluidVelocity,
    const solidAveraging&           parVelocity,
    const Plus::realProcCMField&    dp,
    const Plus::realProcCMField&    Tp,
    Plus::realProcCMField&          Qp,
    const Plus::realProcCMField&    emissivity,
    const Plus::realProcCMField&    radSumTemp,
    const Plus::uint32ProcCMField&  radNumPrt,
    Plus::realProcCMField&          QpRad)
{
    heatInteractionTimer_.start();

    // Delegate the actual Nusselt number and Stefan-Boltzmann calculations
    // to the heatTransfer_ object, passing the appropriate distribution mapper.
    heatTransfer_->calculateHeatTransfer(
        fluidVelocity,
        parVelocity,
        dp, 
        Tp, 
        heatExchangeDistribute_ ? noDistribution_() : uCS().distribution(),
        Qp,
        emissivity,
        radSumTemp,
        radNumPrt,
        QpRad);

    heatInteractionTimer_.end();

    // Print profiling information to the standard output
    Foam::Info << Blue_Text("Heat interaction time: ") 
               << Yellow_Text(heatInteractionTimer_.lastTime())
               << Yellow_Text(" s") << Foam::endl;
}

//+ + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + +

} // coupling
} // pFlow


