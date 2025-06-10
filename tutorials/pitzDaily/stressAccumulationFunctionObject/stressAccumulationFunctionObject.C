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

\*----------------------------------------------------------------------------*/

#include "stressAccumulationFunctionObject.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "pointFields.H"
#include "OSspecific.H"
#include "turbulentTransportModel.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(stressAccumulationFunctionObject, 0);

    addToRunTimeSelectionTable
    (
        functionObject,
        stressAccumulationFunctionObject,
        dictionary
    );
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

bool Foam::stressAccumulationFunctionObject::writeData()
{
    if (mesh_.foundObject<volVectorField>("U"))
    {
        // Lookup the velocity field
        const volVectorField& U =
            mesh_.lookupObject<volVectorField>("U");

        // Lookup the flux field
        const surfaceScalarField& phi =
            mesh_.lookupObject<surfaceScalarField>("phi");

        // Lookup the turbulence model
        const incompressible::turbulenceModel& turbulence =
            mesh_.lookupObject<incompressible::turbulenceModel>
            (
                "turbulenceProperties"
            );

        // Lookup the transport properties
        const IOdictionary& transportProps =
            mesh_.lookupObject<IOdictionary>
            (
                "transportProperties"
            );

        // Lookup the density
        const dimensionedScalar rho("rho", transportProps);

        // Calculate the deviatoric stress
        const volSymmTensorField tau(rho*turbulence.devReff(U));

        // Calculate scalar equivalent stress - scaled version of the von Mises
        // stress (2/3 rather than 3/2)
        sigma_ = sqrt((2.0/3.0)*magSqr(tau));

        // Time step
        const dimensionedScalar deltaT(time_.deltaT());

        // Advect the old sigmaAccum field
        for (int i = 0; i < nNonOrthoCorr_; ++i)
        {
            solve
            (
                fvm::ddt(sigmaAccum_)
              + fvm::div(phi, sigmaAccum_)
            );
        }

        // Add sigmaAccum creation in this time step
        sigmaAccum_ += sigma_*deltaT;

        // Lookup the solution control

        // Advect the old damage field
        for (int i = 0; i < nNonOrthoCorr_; ++i)
        {
            solve
            (
                fvm::ddt(damage_)
              + fvm::div(phi, damage_)
            );
        }

        // Increment damage
        damage_ +=
            Foam::pow(sigma_/sigma0_, r_)
           *deltaT.value()/Foam::pow(1.0 - damage_, k_);

        // Bound damage
        damage_ = max(min(damage_, 1.0), 0.0);

        // Advect the old residenceTime field
        for (int i = 0; i < nNonOrthoCorr_; ++i)
        {
            solve
            (
                fvm::ddt(residenceTime_)
              + fvm::div(phi, residenceTime_)
            );
        }

        // Increment damage
        residenceTime_ += deltaT;

        if (outletPatchID_ != -1)
        {
            // Calculate the integral of stressAccum through the outlet
            integralSigmaAccum_ +=
                gSum
                (
                    phi.boundaryField()[outletPatchID_]
                   *sigmaAccum_.boundaryField()[outletPatchID_]
                );

            // Calculate the integral of damage through the outlet
            integralDamage_ +=
                gSum
                (
                    phi.boundaryField()[outletPatchID_]
                   *damage_.boundaryField()[outletPatchID_]
                );

            // Write outlet integrals to file
            if (Pstream::master())
            {
                historyFilePtr_()
                    << time_.time().value()
                    << " " << integralSigmaAccum_
                    << " " << integralDamage_
                    << endl;
            }
        }
    }
    else
    {
        InfoIn(this->name() + " function object constructor")
            << "U not found" << endl;
    }

    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::stressAccumulationFunctionObject::stressAccumulationFunctionObject
(
    const word& name,
    const Time& t,
    const dictionary& dict
)
:
    functionObject(name),
    name_(name),
    time_(t),
    mesh_(time_.lookupObject<fvMesh>("region0")),
    outletPatchID_(-1),
    historyFilePtr_(),
    sigma_
    (
        IOobject
        (
            "sigma",
            time_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar(dimPressure, 0.0)
    ),
    sigmaAccum_
    (
        IOobject
        (
            "sigmaAccum",
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    damage_
    (
        IOobject
        (
            "damage",
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    residenceTime_
    (
        IOobject
        (
            "residenceTime",
            time_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    sigma0_("sigma0", dict),
    r_("r", dict),
    k_("k", dict),
    nNonOrthoCorr_(readInt(dict.lookup("nNonOrthoCorr"))),
    integralSigmaAccum_(0.0),
    integralDamage_(0.0)
{
    Info<< "Creating " << this->name() << " function object" << endl;

    word outletPatchName("notSpecified");

    if (dict.found("outletPatch"))
    {
        dict.lookup("outletPatch") >> outletPatchName;
    }
    else
    {
        WarningIn(this->name() + " function object constructor")
            << "stressAccumulationFunctionObject: outletPatch not specified" << endl;
    }

    outletPatchID_ = mesh_.boundaryMesh().findPatchID(outletPatchName);

    if (outletPatchID_ == -1)
    {
        WarningIn(this->name() + " function object constructor")
            << "outlet patch " << outletPatchName << " not found"
            << endl;
    }

    // Create outlet file if not already created
    if (historyFilePtr_.empty() && outletPatchID_ != -1)
    {
        // File update
        if (Pstream::master())
        {
            fileName historyDir;

            word startTimeName =
                time_.timeName(time_.startTime().value());

            if (Pstream::parRun())
            {
                // Put in undecomposed case (Note: gives problems for
                // distributed data running)
                historyDir = time_.path()/".."/"postProcessing"/startTimeName;
            }
            else
            {
                historyDir = time_.path()/"postProcessing"/startTimeName;
            }

            // Create directory if does not exist.
            mkDir(historyDir);

            // Open new file at start up
            historyFilePtr_.reset
            (
                new OFstream
                (
                    historyDir/"stressAccumulationFunctionObject"
                   + outletPatchName + ".dat"
                )
            );

            // Add headers to output data
            if (historyFilePtr_.valid())
            {
                historyFilePtr_()
                    << "# Time"
                    << " " << "sigmaAccum"
                    << " " << "damage"
                    << endl;
            }
        }
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::stressAccumulationFunctionObject::start()
{
    return writeData();
}


bool Foam::stressAccumulationFunctionObject::execute()
{
    return writeData();
}


bool Foam::stressAccumulationFunctionObject::read(const dictionary& dict)
{
    return true;
}


bool Foam::stressAccumulationFunctionObject::write()
{
    return false;
}

// ************************************************************************* //
