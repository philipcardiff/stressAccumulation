/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    stressAccumulationFunctionObject

Description
    Function that calculates the accumulated shear stress field, defined as

        stressAccum = integral-in-time sigma d(time)

    where sigma is a scaled-version of the equivalent (von Mises) stress
    tensor.

    This approach is related to the description in Apel et al. (2001),
    Assessment of Hemolysis Related Quantities in a Microaxial Blood Pump by
    Computational Fluid Dynamics, Artificial Organs, 25(5):341–347.

    The residence time and a damage variable are also calculated, where the
    damage is defined in Alemu and Bluestein (2007), Flow-induced Platelet
    Activation and Damage Accumulation in a Mechanical Heart Valve: Numerical
    Studies, Artificial Organs, 31(9):677–688.

Author
    Philip Cardiff, UCD. All rights reserved.

SourceFiles
    stressAccumulationFunctionObject.C

\*---------------------------------------------------------------------------*/

#ifndef stressAccumulationFunctionObject_H
#define stressAccumulationFunctionObject_H

#include "functionObject.H"
#include "dictionary.H"
#include "fvMesh.H"
#include "OFstream.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                        Class stressAccumulationFunctionObject Declaration
\*---------------------------------------------------------------------------*/

class stressAccumulationFunctionObject
:
    public functionObject
{
    // Private data

        //- Name
        const word name_;

        //- Reference to main object registry
        const Time& time_;

        //- Reference to mesh
        const fvMesh& mesh_;

        //- Outlet patch ID
        label outletPatchID_;

        //- History file ptr for writing the outlet integrals
        autoPtr<OFstream> historyFilePtr_;

        //- Instantaneous equivalent stress field
        volScalarField sigma_;

        //- Accmulated stress field
        volScalarField sigmaAccum_;

        //- Damage field
        volScalarField damage_;

        //- Create the residence time field
        volScalarField residenceTime_;

        //- Damage parameters
        const dimensionedScalar sigma0_;
        const dimensionedScalar r_;
        const dimensionedScalar k_;

        //- Number of non-orthogonal correctors
        const label nNonOrthoCorr_;


    // Private Member Functions

        //- Write data
        bool writeData();

        //- Disallow default bitwise copy construct
        stressAccumulationFunctionObject
        (
            const stressAccumulationFunctionObject&
        );

        //- Disallow default bitwise assignment
        void operator=(const stressAccumulationFunctionObject&);


public:

    //- Runtime type information
    TypeName("stressAccumulation");


    // Constructors

        //- Construct from components
        stressAccumulationFunctionObject
        (
            const word& name,
            const Time&,
            const dictionary&
        );


    // Member Functions

        //- start is called at the start of the time-loop
        virtual bool start();

        //- execute is called at each ++ or += of the time-loop
        virtual bool execute();

        //- Called when time was set at the end of the Time::operator++
        virtual bool timeSet()
        {
            return true;
        }

        //- Read and set the function object if its data has changed
        virtual bool read(const dictionary& dict);

        //- Write
        virtual bool write();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
