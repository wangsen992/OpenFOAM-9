/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2015-2021 OpenFOAM Foundation
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

#include "InterfaceCompositionPhaseChangePhaseSystem.H"
#include "interfaceCompositionModel.H"
#include "heatTransferModel.H"
#include "diffusiveMassTransferModel.H"
#include "fvmSup.H"

// * * * * * * * * * * * * Private Member Functions * * * * * * * * * * * * //

template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctDmdtfs()
{
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        *dmdtfs_[pair] = Zero;

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;
                Info << "[InterfaceCompositionPhaseChangePhaseSystem|correctDmdtfs] forAllConstIter on species, pair, pairPhase, index: " 
                  << specie << ", " << pair.name() << ", " << phase.name() 
                  << ", " << pairIter.index()
                  << endl;

                *dmdtfs_[pair] +=
                    (pairIter.index() == 0 ? +1 : -1)
                   *(
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie)
                    );
            }
        }
    }
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::dmidtfTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
totalDmidtfs() const
{
    autoPtr<phaseSystem::dmidtfTable> totalDmidtfsPtr
    (
        new phaseSystem::dmidtfTable
    );
    phaseSystem::dmidtfTable& totalDmidtfs = totalDmidtfsPtr();

    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        if (!totalDmidtfs.found(pair))
        {
            totalDmidtfs.insert(pair, new HashPtrTable<volScalarField>());
        }

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;


                tmp<volScalarField> dmidtf
                (
                    (pairIter.index() == 0 ? +1 : -1)
                   *(
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie)
                    )
                );

                if (totalDmidtfs[pair]->found(specie))
                {
                    *(*totalDmidtfs[pair])[specie] += dmidtf;
                }
                else
                {
                    totalDmidtfs[pair]->insert(specie, dmidtf.ptr());
                }
            }
        }
    }

    return totalDmidtfsPtr;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
InterfaceCompositionPhaseChangePhaseSystem
(
    const fvMesh& mesh
)
:
    BasePhaseSystem(mesh),
    nInterfaceCorrectors_
    (
        this->template lookupOrDefault<label>("nInterfaceCorrectors", 1)
    )
{
    this->generatePairsAndSubModels
    (
        "interfaceComposition",
        interfaceCompositionModels_
    );

    this->generatePairsAndSubModels
    (
        "diffusiveMassTransfer",
        diffusiveMassTransferModels_,
        false
    );

    // Check that models have been specified in the correct combinations
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        this->template validateMassTransfer<interfaceCompositionModel>(pair);

        if (!this->diffusiveMassTransferModels_.found(pair))
        {
            FatalErrorInFunction
                << "A diffusive mass transfer model the " << pair
                << " pair is not specified. This is required by the "
                << "corresponding interface composition model."
                << exit(FatalError);
        }

        forAllConstIter(phasePair, pair, pairIter)
        {
            if
            (
                interfaceCompositionModelIter()[pairIter.index()].valid()
             && !diffusiveMassTransferModels_[pair][pairIter.index()].valid()
            )
            {
                FatalErrorInFunction
                    << "A mass transfer model for the " << (*pairIter).name()
                    << " side of the " << pair << " pair is not "
                    << "specified. This is required by the corresponding "
                    << "interface composition model."
                    << exit(FatalError);
            }
        }

        if
        (
            !this->heatTransferModels_.found(pair)
         || !this->heatTransferModels_[pair].first().valid()
         || !this->heatTransferModels_[pair].second().valid()
        )
        {
             FatalErrorInFunction
                 << "A heat transfer model for both sides of the " << pair
                 << "pair is not specified. This is required by the "
                 << "corresponding interface composition model"
                 << exit(FatalError);
        }
    }

    // Generate mass transfer fields, initially assumed to be zero
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        dmdtfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "interfaceCompositionPhaseChange:dmdtf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, Zero)
            )
        );

        dmidtfSus_.insert(pair, new HashPtrTable<volScalarField>());

        dmidtfSps_.insert(pair, new HashPtrTable<volScalarField>());

        Tfs_.insert
        (
            pair,
            new volScalarField
            (
                IOobject
                (
                    IOobject::groupName
                    (
                        "interfaceCompositionPhaseChange:Tf",
                        pair.name()
                    ),
                    this->mesh().time().timeName(),
                    this->mesh(),
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                (pair.phase1().thermo().T() + pair.phase2().thermo().T())/2
            )
        );

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                dmidtfSus_[pair]->insert
                (
                    specie,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "interfaceCompositionPhaseChange:dmidtfSu",
                                    specie
                                ),
                                pair.name()
                            ),
                            this->mesh().time().timeName(),
                            this->mesh()
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );

                dmidtfSps_[pair]->insert
                (
                    specie,
                    new volScalarField
                    (
                        IOobject
                        (
                            IOobject::groupName
                            (
                                IOobject::groupName
                                (
                                    "interfaceCompositionPhaseChange:dmidtfSp",
                                    specie
                                ),
                                pair.name()
                            ),
                            this->mesh().time().timeName(),
                            this->mesh()
                        ),
                        this->mesh(),
                        dimensionedScalar(dimDensity/dimTime, 0)
                    )
                );
            }
        }
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
~InterfaceCompositionPhaseChangePhaseSystem()
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

template<class BasePhaseSystem>
Foam::tmp<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdtf
(
    const phasePairKey& key
) const
{
    tmp<volScalarField> tDmdtf = BasePhaseSystem::dmdtf(key);

    Info << "[InterfaceCompositionPhaseChangePhaseSystem|dmdtf] key = " << key << ";";
    if (dmdtfs_.found(key))
    {
        const label dmdtSign(Pair<word>::compare(this->phasePairs_[key], key));
        Info << "found key, sign = " << dmdtSign;

        tDmdtf.ref() += dmdtSign**dmdtfs_[key];
    }
    Info << endl;

    return tDmdtf;
}


template<class BasePhaseSystem>
Foam::PtrList<Foam::volScalarField>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::dmdts() const
{
    PtrList<volScalarField> dmdts(BasePhaseSystem::dmdts());

    forAllConstIter(phaseSystem::dmdtfTable, dmdtfs_, dmdtfsIter)
    {
        const phasePair& pair = this->phasePairs_[dmdtfsIter.key()];
        const phaseModel& phase = pair.phase1();
        const phaseModel& otherPhase = pair.phase2();

        addField(phase, "dmdt", *dmdtfsIter(), dmdts);
        addField(otherPhase, "dmdt", - *dmdtfsIter(), dmdts);
    }

    return dmdts;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransfer()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransfer();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::momentumTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
momentumTransferf()
{
    autoPtr<phaseSystem::momentumTransferTable> eqnsPtr =
        BasePhaseSystem::momentumTransferf();

    phaseSystem::momentumTransferTable& eqns = eqnsPtr();

    this->addDmdtUfs(dmdtfs_, eqns);

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::heatTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
heatTransfer() const
{
    Info << "[InterfaceComposition] Enter heatTransfer()" << endl;
    autoPtr<phaseSystem::heatTransferTable> eqnsPtr =
        BasePhaseSystem::heatTransfer();

    phaseSystem::heatTransferTable& eqns = eqnsPtr();

    Info << "[InterfaceComposition] Apply addDmidtHefs()" << endl;
    this->addDmidtHefs
    (
        totalDmidtfs(),
        Tfs_,
        latentHeatScheme::symmetric,
        latentHeatTransfer::mass,
        eqns
    );

    return eqnsPtr;
}


template<class BasePhaseSystem>
Foam::autoPtr<Foam::phaseSystem::specieTransferTable>
Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
specieTransfer() const
{
    Info << "[InterfaceCompositionPhaseChangePhaseSystem] specieTransfer() " << endl;
    autoPtr<phaseSystem::specieTransferTable> eqnsPtr =
        BasePhaseSystem::specieTransfer();

    phaseSystem::specieTransferTable& eqns = eqnsPtr();

    // Explicit
    /*
    this->addDmidtYf(totalDmidtfs(), eqns);
    */

    // Semi-implicit
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];
        Info << "[InterfaceCompositionPhaseChangePhaseSystem] forAllConstIter on interfaceCompositionModels_: " << pair.name() << endl;

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;
            const phaseModel& otherPhase = pairIter.otherPhase();
            Info << "[InterfaceCompositionPhaseChangePhaseSystem] forAllConstIter on phasePair: "
                 << phase.name() << " and " << otherPhase.name() << endl;

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;
                Info << "[InterfaceCompositionPhaseChangePhaseSystem] forAllConstIter on species: "
                     << specie << endl;

                volScalarField dmidtf
                (
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie)
                );
                Info<< "[InterfaceCompositionPhaseChangePhaseSystem] phase.Y(specie) of " 
                    << ": gMin = " <<      gMin(phase.Y(specie).primitiveField())
                    << ", mean = " << gAverage(phase.Y(specie).primitiveField())
                    << ", max = " <<      gMax(phase.Y(specie).primitiveField())
                    << endl;
                Info<< "[InterfaceCompositionPhaseChangePhaseSystem] dmidtf (explicit) " 
                    << ": min = " <<      gMin(dmidtf.primitiveField())
                    << ", mean = " << gAverage(dmidtf.primitiveField())
                    << ", max = " <<      gMax(dmidtf.primitiveField())
                    << endl;

                // Implicit transport through this phase
                // [Test] using explicit transport for species equation here
                *eqns[phase.Y(specie).name()] +=
                    *(*dmidtfSus_[pair])[specie]
                  + fvm::Sp(*(*dmidtfSps_[pair])[specie], phase.Y(specie));

                // Explicit transport out of the other phase
                if (eqns.found(IOobject::groupName(specie, otherPhase.name())))
                {
                    *eqns[otherPhase.Y(specie).name()] -=
                        *(*dmidtfSus_[pair])[specie]
                      + *(*dmidtfSps_[pair])[specie]*phase.Y(specie);
                }
            }
        }
    }

    return eqnsPtr;
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correct()
{
    BasePhaseSystem::correct();

    // Sum up the contribution from each interface composition model
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        const volScalarField& Tf(*this->Tfs_[pair]);

        forAllConstIter(phasePair, pair, pairIter)
        {
            const autoPtr<interfaceCompositionModel>& compositionModelPtr =
                interfaceCompositionModelIter()[pairIter.index()];

            if (!compositionModelPtr.valid()) continue;

            const interfaceCompositionModel& compositionModel =
                compositionModelPtr();

            const phaseModel& phase = *pairIter;

            Info<< "[InterfaceCompositionPhaseChangePhaseSystem] diffusiveMassTransfer " 
                << diffusiveMassTransferModels_[pair][pairIter.index()]->name()
                << endl;
            const volScalarField K
            (
                diffusiveMassTransferModels_[pair][pairIter.index()]->K()
            );

            forAllConstIter
            (
                hashedWordList,
                compositionModel.species(),
                specieIter
            )
            {
                const word& specie = *specieIter;

                const volScalarField KD(K*compositionModel.D(specie));
                const volScalarField Yf(compositionModel.Yf(specie, Tf));
                Info<< "[InterfaceCompositionPhaseChangePhaseSystem] K for  "  << specie
                    << ": min = " <<      min(K.primitiveField())
                    << ", mean = " << average(K.primitiveField())
                    << ", max = " <<      max(K.primitiveField())
                    << endl;
                Info<< "[InterfaceCompositionPhaseChangePhaseSystem] D for  "  << specie
                    << ": min = " <<      min(compositionModel.D(specie)->primitiveField())
                    << ", mean = " << average(compositionModel.D(specie)->primitiveField())
                    << ", max = " <<      max(compositionModel.D(specie)->primitiveField())
                    << endl;
                Info<< "[InterfaceCompositionPhaseChangePhaseSystem] Yf for  "  << specie
                    << ": min = " <<      min(Yf.primitiveField())
                    << ", mean = " << average(Yf.primitiveField())
                    << ", max = " <<      max(Yf.primitiveField())
                    << endl;

                *(*dmidtfSus_[pair])[specie] = phase.rho()*KD*Yf;
                *(*dmidtfSps_[pair])[specie] = - phase.rho()*KD;
            }
        }
    }

    correctDmdtfs();
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctSpecies()
{
    BasePhaseSystem::correctSpecies();

    correctDmdtfs();
}


template<class BasePhaseSystem>
void Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::
correctInterfaceThermo()
{
    // This loop solves for the interface temperatures, Tf, and updates the
    // interface composition models.
    //
    // In the presence of thermally-coupled mass transfer, the relation between
    // heat transfers across the interface between phases 1 and 2 is:
    //
    //                         Q1 + Q2 = mDot*L
    //     H1*(Tf - T1) + H2*(Tf - T1) = K*rho*(Yfi - Yi)*Li
    //
    // Where Q1 and Q2 are the net transfer into phases 1 and 2 respectively,
    // H1 and H2 are the heat transfer coefficients on either side, Tf is the
    // temperature at the interface, mDot is the mass transfer rate from phase
    // 2 to phase 1, and L is the latent heat of phase 2 minus phase 1, K is
    // the diffusive mass transfer coefficient, Yfi - Yi is the concentration
    // difference of a transferring specie between the interface and the bulk
    // driving the transfer, Li is the latent heat change of the specie, and
    // rho is the density in the phase in which the diffusive mass transfer is
    // being represented.
    //
    // Yfi is likely to be a strong non-linear (typically exponential) function
    // of Tf, so the solution for the temperature is newton-accelerated.

    Info << "[InterfaceComposition] Entering correctInterfaceThermo() " << endl;
    BasePhaseSystem::correctInterfaceThermo();

    // First loop 
    forAllConstIter
    (
        interfaceCompositionModelTable,
        interfaceCompositionModels_,
        interfaceCompositionModelIter
    )
    {
        const phasePair& pair =
            this->phasePairs_[interfaceCompositionModelIter.key()];

        Info << "[InterfaceComposition] First loop, Pair " << pair.name() << endl;
        Info << "[InterfaceComposition] heatTransferModels_[pair].first() " << pair.name() << endl;
        const volScalarField H1(this->heatTransferModels_[pair].first()->K());
        Info << "[InterfaceComposition] heatTransferModels_[pair].second() " << pair.name() << endl;
        const volScalarField H2(this->heatTransferModels_[pair].second()->K());
        const dimensionedScalar HSmall("small", heatTransferModel::dimK, small);

        const typename diffusiveMassTransferModelTable::value_type&
            diffusiveMassTransfers = this->diffusiveMassTransferModels_[pair];

        const typename interfaceCompositionModelTable::value_type &
            interfaceCompositions = this->interfaceCompositionModels_[pair];

        volScalarField& Tf = *this->Tfs_[pair];

        for (label i = 0; i < nInterfaceCorrectors_; ++ i)
        {
            Info << "[interfaceComposition] interfaceCorrector: " << i << endl;
            tmp<volScalarField> dmdtLf =
                volScalarField::New
                (
                    IOobject::groupName("dmdtLf", pair.name()),
                    this->mesh(),
                    dimensionedScalar(dimEnergy/dimVolume/dimTime, 0)
                );
            tmp<volScalarField> dmdtLfPrime =
                volScalarField::New
                (
                    IOobject::groupName("dmdtLfPrime", pair.name()),
                    this->mesh(),
                    dimensionedScalar(dmdtLf().dimensions()/dimTemperature, 0)
                );
            tmp<volScalarField> Li =
                volScalarField::New
                (
                    IOobject::groupName("Li", pair.name()),
                    this->mesh(),
                    dimensionedScalar(dimEnergy/dimMass, 0)
                );

            // Add latent heats from forward and backward models
            forAllConstIter(phasePair, pair, pairIter)
            {
                Info << "[InterfaceComposition] Second loop, Pair " << pairIter.index() << " of " << pair.name() << endl;
                if (interfaceCompositions[pairIter.index()].valid())
                {
                    const BlendedInterfacialModel<diffusiveMassTransferModel>&
                        diffusiveMassTransfer =
                        diffusiveMassTransfers[pairIter.index()];

                    const interfaceCompositionModel&
                        interfaceComposition =
                        interfaceCompositions[pairIter.index()];

                    const label sign = pairIter.index() == 0 ? 1 : -1;

                    forAllConstIter
                    (
                        hashedWordList,
                        interfaceComposition.species(),
                        specieIter
                    )
                    {
                        const word& specie = *specieIter;

                        const volScalarField dY
                        (
                            interfaceComposition.dY(specie, Tf)
                        );

                        const volScalarField dYfPrime
                        (
                            interfaceComposition.dYfPrime(specie, Tf)
                        );

                        Info<< "[InterfaceComposition] Tf." << pair.name()
                            << ": min = " << min(Tf.primitiveField())
                            << ", mean = " << average(Tf.primitiveField())
                            << ", max = " << max(Tf.primitiveField())
                            << endl;
                        Li.ref() = this->Li
                                    (
                                        pair,
                                        specie,
                                        dY,
                                        Tf,
                                        latentHeatScheme::symmetric
                                    );
                        Info<< "[InterfaceComposition] Li." << pair.name()
                            << ": min = " <<      min(Li->primitiveField())
                            << ", mean = " << average(Li->primitiveField())
                            << ", max = " <<      max(Li->primitiveField())
                            << endl;
                        Info<< "[InterfaceComposition] dY." << pair.name()
                            << ": gMin = " <<      gMin(dY.primitiveField())
                            << ", mean = " << gAverage(dY.primitiveField())
                            << ", max = " <<      gMax(dY.primitiveField())
                            << endl;
                        Info<< "[InterfaceComposition] dYfPrime." << pair.name()
                            << ": min = " <<      gMin(dYfPrime.primitiveField())
                            << ", mean = " << gAverage(dYfPrime.primitiveField())
                            << ", max = " <<      gMax(dYfPrime.primitiveField())
                            << endl;
                        const volScalarField rhoKDL
                        (
                            pairIter().thermo().rho()
                           *diffusiveMassTransfer.K()
                           *interfaceComposition.D(specie)
                           *this->Li
                            (
                                pair,
                                specie,
                                dY,
                                Tf,
                                latentHeatScheme::symmetric
                            )
                        );

                        dmdtLf.ref() += sign*rhoKDL*dY;
                        dmdtLfPrime.ref() += sign*rhoKDL*dYfPrime;
                    }
                }
            }

            // Update the interface temperature by applying one step of newton's
            // method to the interface relation
            Info<< "[InterfaceComposition] Tf." << pair.name()
                << ": min = " << min(Tf.primitiveField())
                << ", mean = " << average(Tf.primitiveField())
                << ", max = " << max(Tf.primitiveField())
                << endl;
            Info<< "[InterfaceComposition] T1." << pair.name()
                << ": min = " <<      min(pair.phase1().thermo().T().primitiveField())
                << ", mean = " << average(pair.phase1().thermo().T().primitiveField())
                << ", max = " <<      max(pair.phase1().thermo().T().primitiveField())
                << endl;
            Info<< "[InterfaceComposition] T2." << pair.name()
                << ": min = " <<      min(pair.phase2().thermo().T().primitiveField())
                << ", mean = " << average(pair.phase2().thermo().T().primitiveField())
                << ", max = " <<      max(pair.phase2().thermo().T().primitiveField())
                << endl;
            Info<< "[InterfaceComposition] H1." << pair.name()
                << ": min = " <<      min(H1.primitiveField())
                << ", mean = " << average(H1.primitiveField())
                << ", max = " <<      max(H1.primitiveField())
                << endl;
            Info<< "[InterfaceComposition] H2." << pair.name()
                << ": min = " <<      min(H2.primitiveField())
                << ", mean = " << average(H2.primitiveField())
                << ", max = " <<      max(H2.primitiveField())
                << endl;
            Info<< "[InterfaceComposition] dmdtLf." << pair.name()
                << ": min = " <<      min(dmdtLf->primitiveField())
                << ", mean = " << average(dmdtLf->primitiveField())
                << ", max = " <<      max(dmdtLf->primitiveField())
                << endl;
            Info<< "[InterfaceComposition] dmdtLfPrime." << pair.name()
                << ": min = " <<      min(dmdtLfPrime->primitiveField())
                << ", mean = " << average(dmdtLfPrime->primitiveField())
                << ", max = " <<      max(dmdtLfPrime->primitiveField())
                << endl;
            Info << "[InterfaceComposition] HSmall = " << HSmall << endl;
            // tmp<volScalarField> TfDenom =
            //     volScalarField::New
            //     (
            //         IOobject::groupName("TfDenom", pair.name()),
            //         this->mesh(),
            //         dimensionedScalar(H1.dimensions(), 0)
            //     );
            // TfDenom.ref() = H1 + H2 - dmdtLfPrime;
            // Info<< "[InterfaceComposition] H1 + H2 - dmdtLfPrime" << pair.name()
            //     << ": min = " <<      min(TfDenom->primitiveField())
            //     << ", mean = " << average(TfDenom->primitiveField())
            //     << ", max = " <<      max(TfDenom->primitiveField())
            //     << endl;
            Tf -=
                (
                    H1*(Tf - pair.phase1().thermo().T())
                  + H2*(Tf - pair.phase2().thermo().T())
                  - dmdtLf
                )
               /(
                    max(H1 + H2 - dmdtLfPrime, HSmall)
                );

            Tf.correctBoundaryConditions();

            Info<< "[InterfaceComposition] Tf." << pair.name()
                << ": min = " << min(Tf.primitiveField())
                << ", mean = " << average(Tf.primitiveField())
                << ", max = " << max(Tf.primitiveField())
                << endl;

            // Update the interface compositions
            if (this->interfaceCompositionModels_[pair].first().valid())
            {
                this->interfaceCompositionModels_[pair].first()->update(Tf);
            }
            if (this->interfaceCompositionModels_[pair].second().valid())
            {
                this->interfaceCompositionModels_[pair].second()->update(Tf);
            }
        }
    }
}


template<class BasePhaseSystem>
bool Foam::InterfaceCompositionPhaseChangePhaseSystem<BasePhaseSystem>::read()
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
