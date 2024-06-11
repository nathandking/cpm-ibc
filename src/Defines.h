#pragma once
namespace cpm
{
    // #define X_ONE_OVER_DIAGONAL
    // #define SINGLE_STEP
    // #define SHOW_PROGRESS
    // #define SHOW_ITERATION

    // #define USE_IDENTITY_PRECONDITIONER
    // #define USE_DIAGONAL_PRECONDITIONER
    // #define USE_JACOBI_PRECONDITIONER

    enum Preconditioner
    {
        identity,
        diagonal,
        jacobi
    };

// #define TRACK_WHERE

}
