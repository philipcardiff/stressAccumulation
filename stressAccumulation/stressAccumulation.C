/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016-2023 OpenCFD Ltd.
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

Application
    stressAccumulation

Description
    Utility that reads in the velocity field for all the time steps and
    calculates the accumulated shear stress field, defined as

        stressAccum = \int_time sigmaEq dtime

    where sigmaEq is the equivalent (von Mises) stress tensor.

Author
    Philip Cardiff, UCD.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();

    #include "setRootCase.H"
    #include "createTime.H"

    // Lookup the time directories
    instantList timeDirs = timeSelector::select0(runTime, args);

    #include "createMesh.H"

    // Create the stress accumulation field
    volScalarField sigmaAccum
    (
        IOobject
        (
            "sigmaAccum",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Read the transport properties dict
    IOdictionary transportProps
    (
        IOobject
        (
            "transportProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Lookup the kinematic viscosity
    const dimensionedScalar nu
    (
        "nu", dimPressure*dimTime/dimDensity,
        transportProps
    );

    // Lookup the density
    const dimensionedScalar rho("rho", dimDensity, transportProps);

    // Calculate the dynamic viscosity
    const dimensionedScalar mu(nu*rho);

    // Create the simple control
    simpleControl simple(mesh);

    // Loop over the time directories
    forAll(timeDirs, timei)
    {
        // Skip the 0 time
        if (timei == 0)
        {
            continue;
        }

        runTime.setTime(timeDirs[timei], timei);

        Info<< nl << "Time: " << runTime.timeName() << endl;

        // Check for a mesh update
        mesh.readUpdate();

        // Read the velocity field
        Info<< nl << "Reading U" << endl;
        volVectorField U
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        // Calculate the velocity gradient
        const volTensorField gradU(fvc::grad(U));

        // Calculate the deviatoric stress tensor
        const volSymmTensorField tau("tau", mu*dev(twoSymm(gradU)));

        // Calculate equivalent (von Mises) stress
        const volScalarField sigmaEq
        (
            "sigmaEq", sqrt((3.0/2.0)*magSqr(tau))
        );

        // Read the volume flux field
        const surfaceScalarField phi
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            ),
            mesh
        );

        // Advect the old sigmaAccum field
        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(sigmaAccum)
              + fvm::div(phi, sigmaAccum)
            );
        }

        // Add sigmaAccum creation in this time step
        sigmaAccum += sigmaEq*runTime.deltaT();

        // Write the stress fields
        Info<< "Writing tau, sigmaEq and sigmaAccum" << endl;
        tau.write();
        sigmaEq.write();
        sigmaAccum.write();
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
