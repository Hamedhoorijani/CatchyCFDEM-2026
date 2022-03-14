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

#include "segregationExtentRad.H"
#include "fvMesh.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjectsCFDEM
{
    defineTypeNameAndDebug(segregationExtentRad, 0);
    addToRunTimeSelectionTable(functionObjectCFDEM, segregationExtentRad, dictionary);
}
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjectsCFDEM::segregationExtentRad::initialise
(
    const dictionary& dict
)
{}


void Foam::functionObjectsCFDEM::segregationExtentRad::writeFileHeader
(
    const label i
)
{
    writeCommented(file(), "Time");

    file() << tab << "s" << tab << "rav_small" << tab << "rav_large" << tab << "x_small" << tab << "x_large";

    file() << endl;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjectsCFDEM::segregationExtentRad::segregationExtentRad
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    cfdemCloud& sm
)
:
    fvMeshFunctionObjectCFDEM(name, runTime, dict, sm),
    logFiles(obr_, name),
    directionExclude_(readLabel(dict.lookup("directionExclude"))),
    particleCloud_(sm)
{
    read(dict);
    resetName(name);
    createFiles();
}


Foam::functionObjectsCFDEM::segregationExtentRad::segregationExtentRad
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    cfdemCloud& sm
)
:
    fvMeshFunctionObjectCFDEM(name, obr, dict, sm),
    logFiles(obr_, name),
    directionExclude_(readLabel(dict.lookup("directionExclude"))),
    particleCloud_(sm)
{
    read(dict);
    resetName(name);
    createFiles();
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::functionObjectsCFDEM::segregationExtentRad::~segregationExtentRad()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjectsCFDEM::segregationExtentRad::read
(
    const dictionary& dict
)
{
    return true;
}


bool Foam::functionObjectsCFDEM::segregationExtentRad::execute()
{
    return true;
}


bool Foam::functionObjectsCFDEM::segregationExtentRad::write()
{
    Log << type() << "  " << name() << " write " << nl;

    if (Pstream::master())
    {
        writeTime(file());
    }

    scalarField directions_(2,0.0);
    if(directionExclude_ == 0)
    {
        directions_[0] = 1;
        directions_[1] = 2;
    }

    if(directionExclude_ == 1)
    {
        directions_[0] = 0;
        directions_[1] = 2;
    }

    if(directionExclude_ == 2)
    {
        directions_[0] = 0;
        directions_[1] = 1;
    }
    //maximum radial position of GSVU geometry
    dimensionedScalar rmax = dimensionedScalar("rmax",dimensionSet(0,1,0,0,0,0,0),0.04);

    //scalarField r_(particleCloud_.numberOfParticles(),0.0);  
    //scalarField density_(particleCloud_.numberOfParticles(),0.0); 
    double **r_(NULL);
    double **density_(NULL); 

    particleCloud_.dataExchangeM().getData("radius", "scalar-atom", r_);
    particleCloud_.dataExchangeM().getData("density", "scalar-atom", density_);

    scalar rmin_=1e20;
    scalar rmax_=0;

    for(int i = 0; i < particleCloud_.numberOfParticles(); i++)
    {
        if(r_[i][0] < rmin_) rmin_ = r_[i][0];
	if(r_[i][0] > rmax_) rmax_ = r_[i][0];
    }

    scalar rcut_ = (rmax_+rmin_)/2;

    scalar rav_small(0);
    scalar rav_large(0);
    scalar n_small(0);
    scalar rho_small(0);
    scalar n_large(0);
    scalar rho_large(0);

    vector position_(0,0,0);

    for(int i = 0; i < particleCloud_.numberOfParticles(); i++)
    {
        //get radial positions of small and large particle, store in array + calculate average
        position_ = particleCloud_.position(i);

        if(r_[i][0] < rcut_)
        {
            rav_small += sqrt(pow(position_[directions_[0]],2) + pow(position_[directions_[1]],2));
            n_small++; 
            rho_small += density_[i][0];
        }
        else
        {
            rav_large += sqrt(pow(position_[directions_[0]],2) + pow(position_[directions_[1]],2));
            n_large++; 
            rho_large += density_[i][0];
        }

        rav_small /= n_small;
        rav_large /= n_large;
        rho_small /= n_small;
        rho_large /= n_large;
    }

    //determine mass fraction small and large partices
    scalar x_small = n_small*rho_small/(n_small*rho_small + n_large*rho_large);

    if (Pstream::master())
    {
        file()<< tab << (rav_small/rav_large - 1.0)/((2.0-x_small)/(1.0-x_small)-1.0) << tab << rav_small << tab << rav_large << tab << x_small << tab << 1-x_small;
    }

    if (Pstream::master())
    {
        file()<< endl;
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
