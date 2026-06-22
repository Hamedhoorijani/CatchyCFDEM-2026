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
    Foam::rotatingWallVelocityFvPatchVectorField

Description
    This boundary condition provides a rotational velocity condition.

Usage
    \table
        Property     | Description             | Required    | Default value
        origin       | origin of rotation in Cartesian co-ordinates | yes|
        axis         | axis of rotation        | yes         |
        omega        | angular velocity of the frame [rad/s] | yes    |
        value	     | inlet velocity [m/s]
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            rotatingInletVelocity;
        origin          (0 0 0);
        axis            (0 0 1);
        omega           100;
        value 		uniforma (0 0 0);
    }
    \endverbatim

    The \c omega entry is a Function1 of time, see Foam::Function1s.

See also
    Foam::fixedValueFvPatchField
    Foam::Function1s

SourceFiles
    rotatingWallVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef rotatingInletVelocityFvPatchVectorField_H
#define rotatingInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
           Class rotatingWallVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class rotatingInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchField<vector>
{
    // Private Data

        //- Origin of the rotation
        vector origin_;

        //- Axis of the rotation
        vector axis_;

        //- Rotational speed
        autoPtr<Function1<scalar>> omega_;

	// flow control
	autoPtr<Function1<scalar>> volumetricFlowRate_;	//volumetric flow rate [m3/s]
	autoPtr<Function1<scalar>> massFlowRate_;		//mass flow rate [kg/s]
	word rhoName_;						//density of the gas phase
public:

    //- Runtime type information
    TypeName("rotatingInletVelocity");


    // Constructors

        //- Construct from patch and internal field
        rotatingInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        rotatingInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given rotatingWallVelocityFvPatchVectorField
        //  onto a new patch
        rotatingInletVelocityFvPatchVectorField
        (
            const rotatingInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Copy constructor
        rotatingInletVelocityFvPatchVectorField
        (
            const rotatingInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new rotatingInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Copy constructor setting internal field reference
        rotatingInletVelocityFvPatchVectorField
        (
            const rotatingInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Update the coefficients associated with the patch field
        virtual void updateCoeffs();

        //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
