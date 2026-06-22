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

#include "rotatingWallPartialSlipVelocityFvPatchVectorField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

rotatingWallPartialSlipVelocityFvPatchVectorField::
rotatingWallPartialSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    origin_(Zero),
    axis_(Zero),
    omega_(nullptr),
    alpha_(p.size(), scalar(1))
{}


rotatingWallPartialSlipVelocityFvPatchVectorField::
rotatingWallPartialSlipVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF, dict, false),
    origin_(dict.lookup("origin")),
    axis_(dict.lookup("axis")),
    omega_(Function1<scalar>::New("omega", dict)),
    alpha_(p.size(), scalar(1))
{
    if (dict.found("alpha"))
    {
        alpha_ = scalarField("alpha", dict, p.size());
    }

    if (dict.found("value"))
    {
        fvPatchField<vector>::operator=
        (
            vectorField("value", dict, p.size())
        );
    }
    else
    {
        updateCoeffs();
    }
}


rotatingWallPartialSlipVelocityFvPatchVectorField::
rotatingWallPartialSlipVelocityFvPatchVectorField
(
    const rotatingWallPartialSlipVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    origin_(ptf.origin_),
    axis_(ptf.axis_),
    omega_(ptf.omega_, false),
    alpha_(ptf.alpha_)
{}


rotatingWallPartialSlipVelocityFvPatchVectorField::
rotatingWallPartialSlipVelocityFvPatchVectorField
(
    const rotatingWallPartialSlipVelocityFvPatchVectorField& rwvpvf
)
:
    fixedValueFvPatchVectorField(rwvpvf),
    origin_(rwvpvf.origin_),
    axis_(rwvpvf.axis_),
    omega_(rwvpvf.omega_, false),
    alpha_(rwvpvf.alpha_)
{}


rotatingWallPartialSlipVelocityFvPatchVectorField::
rotatingWallPartialSlipVelocityFvPatchVectorField
(
    const rotatingWallPartialSlipVelocityFvPatchVectorField& rwvpvf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(rwvpvf, iF),
    origin_(rwvpvf.origin_),
    axis_(rwvpvf.axis_),
    omega_(rwvpvf.omega_, false),
    alpha_(rwvpvf.alpha_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void rotatingWallPartialSlipVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (!omega_.valid())
    {
        FatalErrorInFunction
            << "omega Function1 is not set"
            << exit(FatalError);
    }

    const scalar axisMag = mag(axis_);
    if (axisMag < VSMALL)
    {
        FatalErrorInFunction
            << "Rotation axis has zero magnitude: axis = " << axis_
            << exit(FatalError);
    }

    const scalar t = this->db().time().timeOutputValue();
    const scalar om = omega_->value(t);

    const vector eAxis = axis_/axisMag;

    const vectorField Cf(patch().Cf());
    const vectorField n(patch().nf());

    // Rigid-body rotating wall velocity
    vectorField Uwall((-om)*((Cf - origin_) ^ eAxis));

    // Tangential part of wall velocity (impermeable)
    Uwall -= n*(n & Uwall);

    // Tangential part of internal field
    const vectorField Ui(this->patchInternalField());
    vectorField Uint(Ui - n*(n & Ui));

    // Blend: alpha=1 -> no-slip rotating wall, alpha=0 -> full slip
    vectorField::operator=(alpha_*Uwall + (scalar(1) - alpha_)*Uint);

    fixedValueFvPatchVectorField::updateCoeffs();
}


void rotatingWallPartialSlipVelocityFvPatchVectorField::write(Ostream& os) const
{
    fvPatchVectorField::write(os);
    writeEntry(os, "origin", origin_);
    writeEntry(os, "axis", axis_);
    writeEntry(os, omega_());
    writeEntry(os, "alpha", alpha_);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchVectorField,
    rotatingWallPartialSlipVelocityFvPatchVectorField
);

} // End namespace Foam

// ************************************************************************* //
