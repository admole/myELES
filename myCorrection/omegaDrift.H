/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2015 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2015 OpenCFD Ltd.
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
    Foam::fv::omegaDrift

Group
    grpFvOptionsConstraints

Description
    Constraint for velocity to dampen velocity fluctuations in
    steady simulations

    In case of velocity exceeding limit applies
    a source sink to remove the excess velocity by adding a term
    that drives the velocity to 0 (note that we cannot use the old
    term velocity since it most likely is already excessive). This
    additional constraint is scaled with (U-UMax)/dt
    where dt is a local time scale (= velocity/cellsize) so it switches off
    if the velocity is below UMax.

    Constraint described by:

        omegaDriftCoeffs
        {
            UMax        100;

            // Optional: name of velocity field (default: U)
            //U         U;
            //UNames    (U);
        }


SourceFiles
    omegaDrift.C

\*---------------------------------------------------------------------------*/

#ifndef omegaDrift_H
#define omegaDrift_H

#include "cellSetOption.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace fv
{

/*---------------------------------------------------------------------------*\
                  Class omegaDrift Declaration
\*---------------------------------------------------------------------------*/

class omegaDrift
:
    public cellSetOption
{

protected:

    // Protected data

        //- Maximum velocity
        vector b_max;
        vector b_min;
        scalar rampStart;
        scalar rampEnd;
        scalar C1_;


    // Protected Member Functions

        void drift(fvMatrix<scalar>& eqn);


private:

    // Private Member Functions

        //- No copy construct
        omegaDrift(const omegaDrift&) = delete;

        //- No copy assignment
        void operator=(const omegaDrift&) = delete;


public:

    //- Runtime type information
    TypeName("omegaDrift");


    // Constructors

        //- Construct from components
        omegaDrift
        (
            const word& name,
            const word& modelType,
            const dictionary& dict,
            const fvMesh& mesh
        );


    //- Destructor
    virtual ~omegaDrift() = default;


    // Member Functions

        // Set values directly

            //- add Sup scalar matrix
            virtual void addSup(fvMatrix<scalar>& eqn, const label fieldi);

            // correct BC
            virtual void correct(volScalarField& omega);


    // I-O

            //- Write data
            virtual void writeData(Ostream& os) const;

            //- Read dictionary
            virtual bool read(const dictionary& dict);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace fv
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
