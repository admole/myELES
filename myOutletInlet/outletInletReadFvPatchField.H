/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::outletInletReadFvPatchField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a generic inflow condition, with
    specified outflow for the case of reverse flow.

Usage
    \table
        Property     | Description             | Required    | Default value
        phi          | Flux field name         | no          | phi
        outletValue  | Outlet value for reverse flow | yes   |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            outletInletRead;
        phi             phi;
        outletValue     uniform 0;
        value           uniform 0;
    }
    \endverbatim

    The mode of operation is determined by the sign of the flux across the
    patch faces.

Note
    Sign conventions:
    - Positive flux (out of domain): apply the "outletValue" fixed-value
    - Negative flux (into of domain): apply zero-gradient condition

See also
    Foam::mixedFvPatchField
    Foam::zeroGradientFvPatchField
    Foam::inletOutletFvPatchField

SourceFiles
    outletInletReadFvPatchField.C

\*---------------------------------------------------------------------------*/

#ifndef outletInletReadFvPatchField_H
#define outletInletReadFvPatchField_H

#include "fvPatchFields.H"
#include "mixedFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                     Class outletInletReadFvPatch Declaration
\*---------------------------------------------------------------------------*/

class outletInletReadFvPatchScalarField
:
    public mixedFvPatchScalarField
{

protected:

    // Protected data

        //- Name of flux field
        word phiName_;
        word FieldName_;
        word GradFieldName_;
        Field<scalar> valueField_;
        Field<scalar> gradField_;
        bool inlet_outlet_;


public:

    //- Runtime type information
    TypeName("outletInletRead");


    // Constructors

        //- Construct from patch and internal field
        outletInletReadFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        outletInletReadFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given outletInletReadFvPatchField onto a new patch
        outletInletReadFvPatchScalarField
        (
            const outletInletReadFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        outletInletReadFvPatchScalarField
        (
            const outletInletReadFvPatchScalarField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchScalarField> clone() const
        {
            return tmp<fvPatchScalarField>
            (
                new outletInletReadFvPatchScalarField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        outletInletReadFvPatchScalarField
        (
            const outletInletReadFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new outletInletReadFvPatchScalarField(*this, iF)
            );
        }


    // Member functions

        const scalarField& valueField() const
        {
            return valueField_;
        }

        const scalarField& gradField() const
        {
            return gradField_;
        }

        scalarField& valueField()
        {
            return valueField_;
        }

        scalarField& gradField()
        {
            return gradField_;
        }

        // Evaluation functions

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
