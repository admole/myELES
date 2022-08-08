/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2013-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::functionObjects::omegaLES

Group
    grpFieldFunctionObjects

Description
    Calculates and outputs the second invariant of the velocity gradient tensor
    [1/s^2].

    \f[
        omegaLES = 0.5(sqr(tr(\nabla U)) - tr(((\nabla U) \cdot (\nabla U))))
    \f]

    where
    \vartable
        U           | velocity [m/s]
    \endvartable

Usage
    Example of function object specification to calculate omegaLES:
    \verbatim
    omegaLES1
    {
        type        omegaLES;
        libs        ("libfieldFunctionObjects.so");
        ...
    }
    \endverbatim

    Where the entries comprise:
    \table
        Property     | Description             | Required    | Default value
        type         | type name: omegaLES            | yes         |
        U            | Name of velocity field  | no          | U
        result       | Name of omegaLES field         | no          | \<function name\>
        log          | Log to standard output  | no          | yes
    \endtable

See also
    Foam::functionObjects::fieldExpression
    Foam::functionObjects::fvMeshFunctionObject

SourceFiles
    omegaLES.C

\*---------------------------------------------------------------------------*/

#ifndef functionObjects_omegaLES_H
#define functionObjects_omegaLES_H

#include "fieldExpression.H"
#include "turbulenceModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{

/*---------------------------------------------------------------------------*\
                          Class omegaLES Declaration
\*---------------------------------------------------------------------------*/

class omegaLES
:
    public fieldExpression
{
    // Private Member Functions

        //- Calculate the omegaLES field and return true if successful
        virtual bool calc();


public:

    //- Runtime type information
    TypeName("omegaLES");


    // Constructors

        //- Construct from Time and dictionary
        omegaLES
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );


    //- Destructor
    virtual ~omegaLES();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
