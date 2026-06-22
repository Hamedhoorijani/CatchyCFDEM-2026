/*---------------------------------------------------------------------------*\
    CFDEMcoupling - Open Source CFD-DEM coupling

    CFDEMcoupling is part of the CFDEMproject
    www.cfdem.com
                                Christoph Goniva, christoph.goniva@cfdem.com
                                Copyright 2009-2012 JKU Linz
                                Copyright 2012-     DCS Computing GmbH, Linz
-------------------------------------------------------------------------------
License
    This file is part of CFDEMcoupling.

    CFDEMcoupling is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 3 of the License, or (at your
    option) any later version.

    CFDEMcoupling is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with CFDEMcoupling; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

    Copyright (C) Behrad Esgandari 2022
    Department of Particulate Flow Modeling
    Johannes Kepler University, Linz, Austria

    Description
    Drag model of Tang et al. for flow past fixed assemblies of monodisperse 
    spheres with direct numerical simulations with the Immersed boundary method.
    Range of applicability: 0 < alphaParticles < 0.6, Re <= 1000

    Reference:
    Tang, Y., Peters, E.A.J.F., Kuipers, J.A.M., Kriebitzsch, S.H.L. and 
    van der Hoef, M.A., 2015. A new drag correlation from fully resolved 
    simulations of flow past monodisperse static arrays of spheres. 
    AIChE journal, 61(2), pp.688-698. 


Class
    TangDrag

SourceFiles
    TangDrag.C

\*---------------------------------------------------------------------------*/

#ifndef TangDrag_H
#define TangDrag_H

#include "forceModel.H"
#include "interpolationCellPoint.H"
//#include "interpolationCellPoint.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class TangDrag Declaration
\*---------------------------------------------------------------------------*/

class TangDrag
:
    public forceModel
{
private:
    dictionary propsDict_;

    word velFieldName_;

    const volVectorField& U_;

    word voidfractionFieldName_;

    const volScalarField& voidfraction_;

    const scalar minVoidfraction_;

    word UsFieldName_;

    const volVectorField& UsField_;  // the average particle velocity field (for implicit/expliti force split)

    mutable scalar scaleDia_;

    mutable scalar scaleDrag_;

    double F(double, double) const;

public:

    //- Runtime type information
    TypeName("TangDrag");


    // Constructors

        //- Construct from components
        TangDrag
        (
            const dictionary& dict,
            cfdemCloud& sm
        );

    // Destructor

        ~TangDrag();


    // Member Functions
        void setForce() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
