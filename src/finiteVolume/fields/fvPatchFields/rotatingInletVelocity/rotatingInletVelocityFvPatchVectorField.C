/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2019 OpenFOAM Foundation
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

#include "rotatingInletVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::rotatingInletVelocityFvPatchVectorField::
rotatingInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF),
    origin_(),
    axis_(Zero),
    omega_(0),
    volumetricFlowRate_(0),
    massFlowRate_(0),
    rhoName_("rho")
{}


Foam::rotatingInletVelocityFvPatchVectorField::
rotatingInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict, false),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(Function1<scalar>::New("omega", dict)),
    volumetricFlowRate_
    (
	    dict.found("volumetricFlowRate")
      ? Function1<scalar>::New("volumetricFlowRate", dict)
      : autoPtr<Function1<scalar>>(nullptr)
    ),
    massFlowRate_
    (
		dict.found("massFlowRate")
      ? Function1<scalar>::New("massFlowRate", dict)
      : autoPtr<Function1<scalar>>(nullptr)
    ),
	rhoName_(dict.lookupOrDefault<word>("rho", "rho"))
{
    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        // Evaluate the wall velocity
        updateCoeffs();
    }
}


Foam::rotatingInletVelocityFvPatchVectorField::
rotatingInletVelocityFvPatchVectorField
(
    const rotatingInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_, false),
    volumetricFlowRate_(ptf.volumetricFlowRate_, false),
    massFlowRate_(ptf.massFlowRate_, false),
    rhoName_(ptf.rhoName_)    
{}


Foam::rotatingInletVelocityFvPatchVectorField::
rotatingInletVelocityFvPatchVectorField
(
    const rotatingInletVelocityFvPatchVectorField& rivpvf
)
:
    fixedValueFvPatchField<vector>(rivpvf),
    origin_(rivpvf.origin_),
    axis_(rivpvf.axis_),
    omega_(rivpvf.omega_, false),
	volumetricFlowRate_(rivpvf.volumetricFlowRate_, false),
    massFlowRate_(rivpvf.massFlowRate_, false),
    rhoName_(rivpvf.rhoName_)
{}


Foam::rotatingInletVelocityFvPatchVectorField::
rotatingInletVelocityFvPatchVectorField
(
    const rotatingInletVelocityFvPatchVectorField& rivpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(rivpvf, iF),
    origin_(rivpvf.origin_),
    axis_(rivpvf.axis_),
    omega_(rivpvf.omega_, false),
	volumetricFlowRate_(rivpvf.volumetricFlowRate_, false),
    massFlowRate_(rivpvf.massFlowRate_, false),
    rhoName_(rivpvf.rhoName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::rotatingInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar t = this->db().time().timeOutputValue();
    scalar om = omega_->value(t); // solid-body rotation like rotatingWallVelocity

    // Calculate the rotating wall velocity from the specification of the motion
    const vectorField Up
    (
        (-om)*((patch().Cf() - origin_) ^ (axis_/mag(axis_)))
    );

    // Remove the component of Up normal to the wall
    // just in case it is not exactly circular
    const vectorField n(patch().nf());
    const vectorField Utan = Up - n *(n&Up); // vectorField Velcity of rotational wall
    
   // -- set uniform normal component from flow rate specification
    const scalarField& magSf = patch().magSf();
    const scalar Atot = gSum(magSf);

    scalar vn = 0.0; // uniform normal speed (negative => into domain; n points outward)

     if (massFlowRate_.valid())
    {
        const scalar mDot = massFlowRate_->value(t);

        // area-weighted rho on the patch
        scalar rhoAvg = 1.0;
        if (this->db().foundObject<volScalarField>(rhoName_))
        {
            const volScalarField& rho =
                this->db().lookupObject<volScalarField>(rhoName_);
            const scalarField rhoP = rho.boundaryField()[patch().index()];
            rhoAvg = gSum(rhoP*magSf)/max(SMALL, Atot);
        }
        else
        {
            WarningInFunction
                << "density field '" << rhoName_
                << "' not found; assuming rho=1 for massFlowRate." << nl;
        }
       // Convention: positive mDot => inflow into the domain
        vn = -(mDot/max(SMALL, rhoAvg))/max(SMALL, Atot);
    }
    else if (volumetricFlowRate_.valid())
    {
        const scalar Q = volumetricFlowRate_->value(t);

        // Convention: positive Q => inflow into the domain
        vn = -(Q/max(SMALL, Atot));
    }
    else
    {
        // no flow-rate specified: just tangential swirl
        vn = 0.0;
    }

    // final inlet velocity
    vectorField::operator=(Utan + n*vn);
    fixedValueFvPatchVectorField::updateCoeffs();
}


void Foam::rotatingInletVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "axis", axis_);
    writeEntry(os, omega_());
    if (massFlowRate_.valid())
    {
        writeEntry(os, massFlowRate_());
        writeEntry(os, "rho", rhoName_);
    }
    else if (volumetricFlowRate_.valid())
    {
        writeEntry(os, volumetricFlowRate_());
    }
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchVectorField,
        rotatingInletVelocityFvPatchVectorField
    );
}

// ************************************************************************* //
