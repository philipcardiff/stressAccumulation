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

        stressAccum = integral-in-time sigma d(time)

    where sigma is a scaled-version of the equivalent (von Mises) stress
    tensor.

    This approach is related to the description in Apel et al. (2001),
    Assessment of Hemolysis Related Quantities in a Microaxial Blood Pump by
    Computational Fluid Dynamics, Artificial Organs, 25(5):341–347.

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

    // Create the damage field
    volScalarField damage
    (
        IOobject
        (
            "damage",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    // Create the residence time field
    volScalarField residenceTime
    (
        IOobject
        (
            "residenceTime",
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

    // Lookup the damage parameters
    const dimensionedScalar sigma0("sigma0", transportProps);
    const dimensionedScalar r("r", transportProps);
    const dimensionedScalar k("k", transportProps);

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

        // Note: this method is only correct is all time steps are writtne to
        // disk!
        const dimensionedScalar deltaT(runTime.deltaT());

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

        // Calculate scalar equivalent stress - scaled version of the von Mises
        // stress (2/3 rather than 3/2)
        const volScalarField sigma
        (
            "sigma", sqrt((2.0/3.0)*magSqr(tau))
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
        sigmaAccum += sigma*deltaT;

        // Advect the old damage field
        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(damage)
              + fvm::div(phi, damage)
            );
        }

        // Increment damage
        damage +=
            Foam::pow(sigma/sigma0, r)
           *deltaT.value()/Foam::pow(1.0 - damage, k);

        // Advect the old residenceTime field
        while (simple.correctNonOrthogonal())
        {
            solve
            (
                fvm::ddt(residenceTime)
              + fvm::div(phi, residenceTime)
            );
        }

        // Increment damage
        residenceTime += deltaT;

        // Write the stress fields
        Info<< "Writing tau, sigma, sigmaAccum, damage and residenceTime to "
            << "Time = " << runTime.timeName() << endl;
        tau.write();
        sigma.write();
        sigmaAccum.write();
        damage.write();
        residenceTime.write();
    }

    Info<< nl << "End" << nl << endl;

    return 0;
}


// ************************************************************************* //
