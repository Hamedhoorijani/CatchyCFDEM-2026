/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#include "segregationExtent.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(segregationExtent, 0);
    addToRunTimeSelectionTable(functionObject, segregationExtent, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::segregationExtent::initialise
(
    const dictionary& dict
)
{}


void Foam::functionObjects::segregationExtent::writeFileHeader
(
    const label i
)
{
    writeCommented(file(), "Time");

    file() << tab << "s";

    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::segregationExtent::segregationExtent
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, runTime, dict),
    logFiles(obr_, name),
    direction_(readLabel(dict.lookup("direction"))),
    phases_(dict.lookup("solidPhases"))
{
    read(dict);
    resetName(name);
    createFiles();
}


Foam::functionObjects::segregationExtent::segregationExtent
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict
)
:
    fvMeshFunctionObject(name, obr, dict),
    logFiles(obr_, name),
    direction_(readLabel(dict.lookup("direction"))),
    phases_(dict.lookup("solidPhases"))
{
    read(dict);
    resetName(name);
    createFiles();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjects::segregationExtent::~segregationExtent()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::segregationExtent::read
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::functionObjects::segregationExtent::execute()
{
    return true;
}


bool Foam::functionObjects::segregationExtent::write()
{
    Log << type() << "  " << name() << " write " << nl;

    if (Pstream::master())
    {
        writeTime(file());
    }

    const volScalarField& alpha_small = mesh_.lookupObject<volScalarField>(IOobject::groupName("alpha", phases_[0]));
    const volScalarField& alpha_large = mesh_.lookupObject<volScalarField>(IOobject::groupName("alpha", phases_[1]));
    scalar vol_small = 0.0;
    scalar vol_large = 0.0;
    scalar h_small = 0.0;
    scalar h_large = 0.0;

/*  volScalarField vol_smalli = alpha_small*mesh_.V();
  volScalarField vol_largei = alpha_large*mesh_.V();
  volScalarfield h_smalli = alpha_small*mesh_.V()*mesh_.C()[direction_];
  volScalarField h_largei = alpha_large*mesh_.V()*mesh_.C()[direction_];
  
	scalar vol_small = gSum(vol_smalli);
 	scalar vol_large = gSum(vol_largei);
 	scalar h_small = gSum(h_smalli);
 	scalar h_large = gSum(h_largei); */

    forAll(alpha_small, celli)
    {
        vol_small += alpha_small[celli]*mesh_.V()[celli];
        vol_large += alpha_large[celli]*mesh_.V()[celli];
        h_small += alpha_small[celli]*mesh_.V()[celli]*mesh_.C()[celli][direction_];
        h_large += alpha_large[celli]*mesh_.V()[celli]*mesh_.C()[celli][direction_];
    }

    h_small /= vol_small;
    h_large /= vol_large;
    scalar x_small = vol_small/(vol_small + vol_large);

    if (Pstream::master())
    {
        file()<< tab << (h_small/h_large - 1.0)/((2.0-x_small)/(1.0-x_small)-1.0);
    }

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
