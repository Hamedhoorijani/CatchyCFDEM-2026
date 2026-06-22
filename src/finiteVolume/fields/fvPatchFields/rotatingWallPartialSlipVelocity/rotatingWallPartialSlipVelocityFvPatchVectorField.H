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

Class
    Foam::rotatingWallPartialSlipVelocityFvPatchVectorField

Description
    Rotating-wall velocity boundary condition with controllable slip.

    The wall tangential velocity is computed from rigid-body rotation:
        U_wall = omega(t) x (x - origin)

    The normal component is removed to enforce impermeability.

    The applied patch velocity is a blend between the rotating-wall tangential
    velocity and the internal-field tangential velocity:
        U_patch = alpha * U_wall,t + (1 - alpha) * U_int,t

    where alpha is in [0,1].
        alpha = 1 -> rotating no-slip
        alpha = 0 -> full slip (rotation does not drive tangential velocity)

Usage
    \table
        Property | Description | Required | Default
        origin   | Origin point on rotation axis | yes | -
        axis     | Rotation axis direction (non-zero) | yes | -
        omega    | Angular speed Function1 [rad/s] | yes | -
        alpha    | Slip blending factor [0..1] (scalarField) | no | 1
    \endtable

    Example:
    \verbatim
    <patchName>
    {
        type        rotatingWallPartialSlipVelocity;
        origin      (0 0 0);
        axis        (0 0 1);

        omega       constant 100;      // rad/s (Function1)

        alpha       uniform 0.2;       // 0 slip, 1 no-slip

        value       uniform (0 0 0);   // placeholder
    }
    \endverbatim

    The \c omega entry is a Function1 of time, see Foam::Function1s.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s

SourceFiles
    rotatingWallPartialSlipVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingWallPartialSlipVelocityFvPatchVectorField_H
#define rotatingWallPartialSlipVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
   Class rotatingWallPartialSlipVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class rotatingWallPartialSlipVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private Data

        //- Origin point on the rotation axis
        vector origin_;

        //- Axis of rotation (direction, need not be unit length)
        vector axis_;

        //- Angular speed (Function1 of time), [rad/s]
        autoPtr<Function1<scalar>> omega_;

        //- Slip blending factor [0..1] (0=slip, 1=no-slip), per-face
        scalarField alpha_;


public:

    //- Runtime type information
    TypeName("rotatingWallPartialSlipVelocity");


    // Constructors

        //- Construct from patch and internal field
        rotatingWallPartialSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rotatingWallPartialSlipVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping onto a new patch
        rotatingWallPartialSlipVelocityFvPatchVectorField
        (
            const rotatingWallPartialSlipVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        rotatingWallPartialSlipVelocityFvPatchVectorField
        (
            const rotatingWallPartialSlipVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new rotatingWallPartialSlipVelocityFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        rotatingWallPartialSlipVelocityFvPatchVectorField
        (
            const rotatingWallPartialSlipVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new rotatingWallPartialSlipVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member Functions

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};

} // End namespace Foam

#endif

// ************************************************************************* //
