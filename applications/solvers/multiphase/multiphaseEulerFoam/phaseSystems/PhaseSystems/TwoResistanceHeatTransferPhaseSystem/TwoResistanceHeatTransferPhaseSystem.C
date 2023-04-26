/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "TwoResistanceHeatTransferPhaseSystem.H"
#include "BlendedInterfacialModel.H"
#include "heatTransferModel.H"
#include "fvmSup.H"
#include "rhoReactionThermo.H"

// * * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefs
(
    const phaseSystem::dmdtfTable& dmdtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const latentHeatScheme scheme,
    const latentHeatTransfer transfer,
    phaseSystem::heatTransferTable& eqns
) const
{
    HeatTransferPhaseSystem<BasePhaseSystem>::addDmdtHefsWithoutL
    (
        dmdtfs,
        Tfs,
        scheme,
        eqns
    );

    // Loop the pairs
    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs, dmdtfIter)
    {
        const phasePairKey& key = dmdtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const volScalarField dmdtf(Pair<word>::compare(pair, key)**dmdtfIter());

        const volScalarField& Tf = *Tfs[key];

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();

        // Transfer coefficients
        const volScalarField H1(heatTransferModels_[key].first()->K());
        const volScalarField H2(heatTransferModels_[key].second()->K());
        const volScalarField H1Fac(H1/(H1 + H2));
        const volScalarField HEff(H1Fac*H2);

        // Latent heat contribution
        switch (transfer)
        {
            case latentHeatTransfer::heat:
            {
                *eqns[phase1.name()] +=
                  - HEff*(thermo2.T() - thermo1.T()) + H1*(Tf - thermo1.T());

                *eqns[phase2.name()] +=
                  - HEff*(thermo1.T() - thermo2.T()) + H2*(Tf - thermo2.T());

                break;
            }
            case latentHeatTransfer::mass:
            {
                const volScalarField L(this->L(pair, dmdtf, Tf, scheme));

                *eqns[phase1.name()] += H1Fac*dmdtf*L;
                *eqns[phase2.name()] += (1 - H1Fac)*dmdtf*L;

                break;
            }
        }
    }
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefs
(
    const phaseSystem::dmidtfTable& dmidtfs,
    const phaseSystem::dmdtfTable& Tfs,
    const latentHeatScheme scheme,
    const latentHeatTransfer transfer,
    phaseSystem::heatTransferTable& eqns
) const
{
    HeatTransferPhaseSystem<BasePhaseSystem>::addDmidtHefsWithoutL
    (
        dmidtfs,
        Tfs,
        scheme,
        eqns
    );

    // Loop the pairs
    forAllConstIter(phaseSystem::dmidtfTable, dmidtfs, dmidtfIter)
    {
        const phasePairKey& key = dmidtfIter.key();
        const phasePair& pair(this->phasePairs_[key]);

        const volScalarField& Tf = *Tfs[key];

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        const rhoThermo& thermo1 = phase1.thermo();
        const rhoThermo& thermo2 = phase2.thermo();

        // Transfer coefficients
        const volScalarField H1(heatTransferModels_[key].first()->K());
        const volScalarField H2(heatTransferModels_[key].second()->K());
        const volScalarField H1Fac(H1/(H1 + H2));
        const volScalarField HEff(H1Fac*H2);

        // Loop the species
        forAllConstIter(HashPtrTable<volScalarField>, *dmidtfIter(), dmidtfJter)
        {
            const word& specie = dmidtfJter.key();

            // Mass transfer rates
            const volScalarField dmidtf
            (
                Pair<word>::compare(pair, key)**dmidtfJter()
            );

            // Latent heat contribution
            switch (transfer)
            {
                case latentHeatTransfer::heat:
                {
                    // Do nothing. This term is handled outside the specie loop.

                    break;
                }
                case latentHeatTransfer::mass:
                {
                    const volScalarField Li
                    (
                        this->Li(pair, specie, dmidtf, Tf, scheme)
                    );

                    *eqns[phase1.name()] += H1Fac*dmidtf*Li;
                    *eqns[phase2.name()] += (1 - H1Fac)*dmidtf*Li;

                    break;
                }
            }
        }

        // Latent heat contribution
        switch (transfer)
        {
            case latentHeatTransfer::heat:
            {
                *eqns[phase1.name()] +=
                  - HEff*(thermo2.T() - thermo1.T()) + H1*(Tf - thermo1.T());

                *eqns[phase2.name()] +=
                  - HEff*(thermo1.T() - thermo2.T()) + H2*(Tf - thermo2.T());

                break;
            }
            case latentHeatTransfer::mass:
            {
                // Do nothing. This term is handled inside the specie loop.

                break;
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
TwoResistanceHeatTransferPhaseSystem
(
    const fvMesh& mesh
)
:
    HeatTransferPhaseSystem<BasePhaseSystem>(mesh)
{
    Info << "[TwoResistance] TwoResistanceHeatTransferPhaseSystem Init." << endl;
    this->generatePairsAndSubModels
    (
        "heatTransfer",
        heatTransferModels_,
        false
    );
    Info << "[TwoResistance] size of heatTransferModels_" << heatTransferModels_.size() << endl;
    // Check that models have been specified on both sides of the interfaces
    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];
        Info << "[TwoResistance] " << pair << " model: " << heatTransferModels_[pair].first()->name() << " and " << heatTransferModels_[pair].second()->name() << endl;

        if (!heatTransferModels_[pair].first().valid())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase1().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
        if (!heatTransferModels_[pair].second().valid())
        {
            FatalErrorInFunction
                << "A heat transfer model for the " << pair.phase2().name()
                << " side of the " << pair << " pair is not specified"
                << exit(FatalError);
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
~TwoResistanceHeatTransferPhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    Info << "[TwoResistance] Enter heatTransfer()" << endl;
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr
    (
        new phaseSystem::heatTransferTable()
    );

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    forAll(this->phaseModels_, phasei)
    {
        const phaseModel& phase = this->phaseModels_[phasei];

        eqns.insert
        (
            phase.name(),
            new fvScalarMatrix(phase.thermo().he(), dimEnergy/dimTime)
        );
    }

    forAllConstIter
    (
        heatTransferModelTable,
        heatTransferModels_,
        heatTransferModelIter
    )
    {
        Info << "[TwoResistance] Calculate heatTransfer Iter" << endl;
        const phasePair& pair = this->phasePairs_[heatTransferModelIter.key()];

        const phaseModel& phase1 = pair.phase1();
        const phaseModel& phase2 = pair.phase2();
        Info << "[TwoResistance] Calculate he" << endl;
        const volScalarField& he1 = phase1.thermo().he();
        const volScalarField& he2 = phase2.thermo().he();
        Info << "[TwoResistance] Calculate Cpv" << endl;
        const volScalarField Cpv1(phase1.thermo().Cpv());
        const volScalarField Cpv2(phase2.thermo().Cpv());

        Info << "[TwoResistance] Calculate K" << endl;
        const volScalarField H1(heatTransferModelIter().first()->K());
        const volScalarField H2(heatTransferModelIter().second()->K());
        Info << "[TwoResistance] Calculate HEff" << endl;
        Info << "[TwoResistance]min(H1) = " << min(H1) << endl;
        Info << "[TwoResistance]min(H2) = " << min(H2) << endl;
        Info << "[TwoResistance]min(H1+H2) = " << min(H1 + H2) << endl;
        const volScalarField HEff(H1*H2/(H1 + H2));
        Info<< "[TwoResistance] he1." << pair.name()
            << ": min = " <<      min(he1.primitiveField())
            << ", mean = " << average(he1.primitiveField())
            << ", max = " <<      max(he1.primitiveField())
            << endl;
        Info<< "[TwoResistance] he2." << pair.name()
            << ": min = " <<      min(he2.primitiveField())
            << ", mean = " << average(he2.primitiveField())
            << ", max = " <<      max(he2.primitiveField())
            << endl;
        Info<< "[TwoResistance] T1." << pair.name()
            << ": min = " <<      min(pair.phase1().thermo().T().primitiveField())
            << ", mean = " << average(pair.phase1().thermo().T().primitiveField())
            << ", max = " <<      max(pair.phase1().thermo().T().primitiveField())
            << endl;
        Info<< "[TwoResistance] T2." << pair.name()
            << ": min = " <<      min(pair.phase2().thermo().T().primitiveField())
            << ", mean = " << average(pair.phase2().thermo().T().primitiveField())
            << ", max = " <<      max(pair.phase2().thermo().T().primitiveField())
            << endl;
        Info<< "[TwoResistance] HEff." << pair.name()
            << ": min = " <<      min(HEff.primitiveField())
            << ", mean = " << average(HEff.primitiveField())
            << ", max = " <<      max(HEff.primitiveField())
            << endl;

        Info << "[TwoResistance] Adding heatTrasfer source term to equations. " << endl;
        Info << "[TwoResistance] Adding source term to " << phase1.name() << endl;
        *eqns[phase1.name()] +=
            HEff*(phase2.thermo().T() - phase1.thermo().T())
           + H1/Cpv1*he1 - fvm::Sp(H1/Cpv1, he1);

        Info << "[TwoResistance] Adding source term to " << phase2.name() << endl;
        *eqns[phase2.name()] +=
            HEff*(phase1.thermo().T() - phase2.thermo().T())
           + H2/Cpv2*he2 - fvm::Sp(H2/Cpv2, he2);
        Info << "[TwoResistance] Complete adding heatTrasfer source term to equations. " << endl;
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::
correctEnergyTransport()
{
    BasePhaseSystem::correctEnergyTransport();

    Info << "[TwoResistance] Before entering correctInterfaceThermo()" << endl;
    correctInterfaceThermo();
    Info << "[TwoResistance] After entering correctInterfaceThermo()" << endl;
}


template<class BasePhaseSystem>
bool Foam::TwoResistanceHeatTransferPhaseSystem<BasePhaseSystem>::read()
{
    if (BasePhaseSystem::read())
    {
        bool readOK = true;

        // Models ...

        return readOK;
    }
    else
    {
        return false;
    }
}


// ************************************************************************* //
