/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2019 OpenFOAM Foundation
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
    Foam::functionObjects::fvMeshFunctionObject

Description
    Specialization of Foam::functionObject for an Foam::fvMesh, providing a
    reference to the Foam::fvMesh.

    If the selected region is not an Foam::fvMesh a Foam::FatalError will be
    generated.

See also
    Foam::regionFunctionObject
    Foam::functionObject

SourceFiles
    fvMeshFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_fvMeshFunctionObjectCFDEM_H
#define functionObjects_fvMeshFunctionObjectCFDEM_H

#include "regionFunctionObject.H"
#include "cfdemCloud.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes
class fvMesh;

namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                     Class fvMeshFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class fvMeshFunctionObjectCFDEM
:
    public regionFunctionObject
{

protected:

    // Protected member data

        //- Reference to the fvMesh
        const fvMesh& mesh_;
        const cfdemCloud& particleCloud_;


public:

    //- Runtime type information
    TypeName("fvMeshFunctionObjectCFDEM");


    // Constructors

        //- Construct from Time and dictionary
        fvMeshFunctionObjectCFDEM
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );

        //- Construct from the region objectRegistry and dictionary
        fvMeshFunctionObjectCFDEM
        (
            const word& name,
            const objectRegistry& obr,
            const dictionary& dict
        );

        //- Disallow default bitwise copy construction
        fvMeshFunctionObjectCFDEM(const fvMeshFunctionObjectCFDEM&) = delete;


    //- Destructor
    virtual ~fvMeshFunctionObjectCFDEM();


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const fvMeshFunctionObjectCFDEM&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
