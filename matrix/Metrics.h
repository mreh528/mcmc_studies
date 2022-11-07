#ifndef _METRICS_H_
#define _METRICS_H_

#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TVectorD.h"
#include "TMatrixD.h"
#include "TMatrixDSym.h"
#include "TMatrixDEigen.h"
#include "TMatrixDSymEigen.h"

#include <iostream>

// ---- Linear algebra functions needed for matrix similarity metrics ---- //

// Calculates the matrix square root of a square matrix via eigen decomposition
TMatrixD SqrtMat(TMatrixD* mat) {
    if (mat->GetNcols() != mat->GetNrows()) {
        std::cout << "ERROR: matrix must be square for sqrt function." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Get the eigenvalues and eigenvectors
    TVectorD eigenvals(mat->GetNcols());
    TMatrixD eigenvecs(mat->EigenVectors(eigenvals));
    TMatrixD eigenvecs_inverted(eigenvecs);
    eigenvecs_inverted.Invert();

    // Put the eigenvalues into diagonal matrix form and sqrt the elements
    TMatrixD diag_mat(mat->GetNcols(), mat->GetNcols());
    diag_mat.Zero();
    for (int i = 0; i < mat->GetNcols(); ++i) {
        diag_mat(i,i) = TMath::Sqrt(eigenvals(i));
    }

    // Return the eigen-decomposed square root: (V*sqrt(D)*V^-1)
    eigenvecs *= diag_mat;
    eigenvecs *= eigenvecs_inverted;
    return eigenvecs;
}


// Calculates the log of a square matrix via eigen decomposition
TMatrixD LogMat(TMatrixD* mat) {
    if (mat->GetNcols() != mat->GetNrows()) {
        std::cout << "ERROR: matrix must be square for log function." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Get the eigenvalues and eigenvectors
    TVectorD eigenvals(mat->GetNcols());
    TMatrixD eigenvecs(mat->EigenVectors(eigenvals));
    TMatrixD eigenvecs_inverted(eigenvecs);
    eigenvecs_inverted.Invert();

    // Put the eigenvalues into diagonal matrix form and take the log of the elements
    TMatrixD diag_mat(mat->GetNcols(), mat->GetNcols());
    diag_mat.Zero();
    for (int i = 0; i < mat->GetNcols(); ++i) {
        diag_mat(i,i) = TMath::Log(eigenvals(i));
    }

    // Return the eigen-decomposed matrix log: (V*log(D)*V^-1)
    eigenvecs *= diag_mat;
    eigenvecs *= eigenvecs_inverted;
    return eigenvecs;
}


// Calculates the L2 Norm of a vector
double L2Norm(TVectorD* vec) {
    double total = 0.;
    for (int i = 0; i < vec->GetNrows(); ++i) {
        total += (*vec)(i) * (*vec)(i);
    }

    return TMath::Sqrt(total);
}


// Calculates the Frobenius Norm of a matrix (similar to L2 vector norm, but for matrices)
double FrobeniusNorm(TMatrixD* mat) {
    double total = 0.;
    for (int i = 0; i < mat->GetNrows(); ++i) {
        for (int j = 0; j < mat->GetNcols(); ++j) {
            total += (*mat)(i,j) * (*mat)(i,j);
        }
    }

    return TMath::Sqrt(total);
}


// Calculates the trace of a square matrix
double TraceMat(TMatrixD* mat) {
    if (mat->GetNrows() != mat->GetNcols()) {
        std::cout << "ERROR: Will only take the trace of a square matrix." << std::endl;
        exit(EXIT_FAILURE);
    }

    // Sum the diagonal elements
    double total = 0.;
    for (int i = 0; i < mat->GetNrows(); ++i) {
        total += (*mat)(i,i);
    }

    return total;
}


// ---- Metrics for calculating similiarity between two distributions ---- //

// Calculates the AIRM distance between two matrices.
// Matrices must be of the same square shape.
double AffineInvariantRiemannianMetric(TMatrixD* mat1, TMatrixD* mat2) {
    if (mat1->GetNcols() != mat2->GetNcols() ||
        mat1->GetNrows() != mat2->GetNrows() ||
        mat1->GetNcols() != mat1->GetNrows()) {
        std::cout << "ERROR: Input matrices for calculating AIRM must be of the same square shape." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Calculating AIRM between " << mat1->GetName() << " and "
              << mat2->GetName() << "..." << std::endl;

    TMatrixD* mat1_sqrt = new TMatrixD(SqrtMat(mat1));
    mat1_sqrt->Invert();
    TMatrixD* tmp = new TMatrixD((*mat1_sqrt)*((*mat2)*(*mat1_sqrt)));
    TMatrixD* logtmp = new TMatrixD(LogMat(tmp));

    return FrobeniusNorm(logtmp);
}


// Calculates the LERM distance between two matrices.
// Matrices must be of the same square shape.
double LogEuclideanRiemannMetric(TMatrixD* mat1, TMatrixD* mat2) {
    if (mat1->GetNcols() != mat2->GetNcols() ||
        mat1->GetNrows() != mat2->GetNrows() ||
        mat1->GetNcols() != mat1->GetNrows()) {
        std::cout << "ERROR: Input matrices for calculating LERM must be of the same square shape." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Calculating LERM between " << mat1->GetName() << " and "
              << mat2->GetName() << "..." << std::endl;

    TMatrixD* mat1_log = new TMatrixD(LogMat(mat1));
    TMatrixD* mat2_log = new TMatrixD(LogMat(mat2));
    (*mat1_log) -= (*mat2_log);

    return FrobeniusNorm(mat1_log);
}


// Calculates the WM distance between two multivariate normal distributions.
// Distributions must be the same dimentionality.
double WassersteinMetric(TMatrixD* cov1, TVectorD* mean1, TMatrixD* cov2, TVectorD* mean2) {
    if (cov1->GetNrows() != cov1->GetNcols() ||
        cov1->GetNrows() != cov2->GetNrows() ||
        cov1->GetNrows() != cov2->GetNcols() ||
        cov1->GetNrows() != mean1->GetNrows() ||
        cov1->GetNrows() != mean2->GetNrows()) {
        std::cout << "ERROR: Input matrices and vectors for calculating Wasserstein must have the same dimensionality." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Calculating Wasserstein Metric between " << cov1->GetName() << " and "
              << cov2->GetName() << "..." << std::endl;

    // Prepare matrices and diff vector -- sorry for the mess :/
    TMatrixD* tmp1 = new TMatrixD(SqrtMat(cov2));
    TMatrixD* tmp2 = new TMatrixD((*tmp1) * ((*cov1) * (*tmp1)));
    TMatrixD* tmp3 = new TMatrixD(2.*SqrtMat(tmp2));
    TMatrixD* tmp4 = new TMatrixD((*cov1) + (*cov2) - (*tmp3));
    TVectorD* mean_diff = new TVectorD((*mean1) - (*mean2));

    // Calculate Wasserstein metric
    double total = L2Norm(mean_diff);
    total *= total;
    total += TraceMat(tmp4);

    return TMath::Sqrt(TMath::Abs(total));
}


// Calculates the Roberts-Rosenthal suboptimality factor between two matrices.
// Matrices must be of the same square shape.
double RRSuboptimalityFactor(TMatrixD* mat1, TMatrixD* mat2) {
    if (mat1->GetNcols() != mat2->GetNcols() ||
        mat1->GetNrows() != mat2->GetNrows() ||
        mat1->GetNcols() != mat1->GetNrows()) {
        std::cout << "ERROR: Input matrices for calculating RR Suboptimality must be of the same square shape." << std::endl;
        exit(EXIT_FAILURE);
    }

    std::cout << "Calculating Suboptimality Factor between " << mat1->GetName() << " and "
              << mat2->GetName() << "..." << std::endl;

    // Construct the matrix tmp = mat1^1/2 * mat2^-1/2
    TMatrixD* mat1_sqrt = new TMatrixD(SqrtMat(mat1));
    TMatrixD* mat2_sqrt = new TMatrixD(SqrtMat(mat2));
    mat2_sqrt->Invert();
    TMatrixD* tmp = new TMatrixD((*mat1_sqrt) * (*mat2_sqrt));

    // Get the eigenvalues (vectors not needed)
    int N = mat1->GetNcols();
    TVectorD* eigenvals = new TVectorD(N);
    tmp->EigenVectors(*eigenvals);

    // Evaluate the Roberts-Rosenthal sum
    double numerator = 0.;
    double denominator = 0.;
    for (int i = 0; i < N; ++i) {
        numerator += 1. / ((*eigenvals)(i) * (*eigenvals)(i));
        denominator += 1. / (*eigenvals)(i);
    }
    denominator *= denominator;

    return (double)N * numerator / denominator;
}

#endif

