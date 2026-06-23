//==============================================================================
//  mmS.H — build Qe_mms, Ql_mms at the current time (uses mmsParams.H).
//==============================================================================
{
    scalar tNow  = runTime.value();
    scalar decay = Foam::exp(-tNow/tau_mms);

    forAll(mesh.C(), cellI)
    {
        scalar z   = mesh.C()[cellI].z() - z0_mms;
        scalar shp = Foam::cos(kz_mms*z);
        Qe_mms[cellI] = coefQe * shp * decay;
        Ql_mms[cellI] = coefQl * shp * decay;
    }
}

