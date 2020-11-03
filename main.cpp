#include <iostream>
#include <math.h>
#include <fstream>
#include <string.h>
#include <vector> 
#include <bits/stdc++.h> 

using namespace std;

int main() {

    int num_nodes, num_cells, num_faces, num_intFaces;
    int num_bc = 3, faceSum = 1, num_totBC = 0, sol_time = 0;		// initializing Number of BC exists in Mesh
    int num_bcFaces[num_bc];
    int **faces;
    double **vertices;
    ifstream myfile("flowOverCyl.msh");								// Name of Mesh File Name
	#include "readMesh.h"											// Reading Mesh
    cout << " Mesh Reading Finished Successfully " << endl;
	#include "preProcess.h"											// Pre Processing the Mesh
    cout << " PreProcessing Finished Successfully " << endl;

	// Initilaing Problem Parameters
    double Re, dt, nu, rho, d_length, U_freeStream;
    Re 	= 100; 		// Reynolds number
    dt 	= 1e-02; 	// time step
    rho = 1.0; 		// density
    nu 	= 1 / Re;	// kinematic Viscosity
    U_freeStream = 1.0;

    vector<double> u_curr(num_cells + 1 + num_totBC, U_freeStream);		// Current Velocity
    vector<double> v_curr(num_cells + 1 + num_totBC, 0);
    vector<double> u_pred(num_cells + 1 + num_totBC, U_freeStream);		// Predicted Velocity
    vector<double> v_pred(num_cells + 1 + num_totBC, 0);
    vector<double> u_Next(num_cells + 1 + num_totBC, U_freeStream);		// Corrected Velocity for Next Time Step
    vector<double> v_Next(num_cells + 1 + num_totBC, 0);
    vector<double> p_Next(num_cells + 1 + num_totBC, 0);				// Pressure
    
    vector<double> u_gradX(num_cells + 1 + num_totBC, 0);				// Gradients in X and Y for Velocity
    vector<double> u_gradY(num_cells + 1 + num_totBC, 0);
    vector<double> v_gradX(num_cells + 1 + num_totBC, 0);
    vector<double> v_gradY(num_cells + 1 + num_totBC, 0);

    vector<double> uc_pred(num_nodes + 1, 0);							// Nodal Values of Velocity and Pressure for Export
    vector<double> vc_pred(num_nodes + 1, 0);
    vector<double> uc_Next(num_nodes + 1, 0);
    vector<double> vc_Next(num_nodes + 1, 0);
    vector<double> pc_Next(num_nodes + 1, 0);

    vector<vector<double> > flux(num_cells + 1, vector<double> (4, 0));		// flux or flow rate
    vector<vector<double> > pdiff(num_cells + 1, vector<double> (4, 0));

    vector<double> udiff(4, 0);
    vector<double> vdiff(4, 0);
    vector<double> PSx(4, 0);
    vector<double> PSy(4, 0);

    vector<double> flux_noPressureCorrection(4, 0);
    vector<double> uu(4, 0);
    vector<double> vv(4, 0);

    int iter = 0;
    double error = 1;

    cout << " Starting time loop \n\n " << endl;

    while (error > 1e-07) {							// Time Loop for Steady state or Transient time 

        iter = iter + 1;
        cout << " Time = " << iter * dt << endl;	// Printing Solution Time

        // Predicted Velocity
        
        int iteration = 0;
        double residual_error = 1;
        double uPred_residual, vPred_residual, uPred_error, vPred_error;

        while (residual_error > 1e-09) {
            uPred_residual = 0;
            vPred_residual = 0;

            iteration = iteration + 1;

            for (int i = 1; i <= num_cells; i++) {		// Calculating Cell Gradients for Linear Upwind

                u_gradX[i] = 0.0;
                u_gradY[i] = 0.0;
                v_gradX[i] = 0.0;
                v_gradY[i] = 0.0;

                for (int k = 0; k < 4; k++) {

                    u_gradX[i] = u_gradX[i] + (u_pred[i] * faceWt[i][k] + u_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k] / Vp[i];
                    u_gradY[i] = u_gradY[i] + (u_pred[i] * faceWt[i][k] + u_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k] / Vp[i];
                    v_gradX[i] = v_gradX[i] + (v_pred[i] * faceWt[i][k] + v_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k] / Vp[i];
                    v_gradY[i] = v_gradY[i] + (v_pred[i] * faceWt[i][k] + v_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k] / Vp[i];

                }
            }

            for (int i = 1; i <= num_cells; i++) {

                double convection_CentralCoeff = 0;
                double uPredConvection = 0, vPredConvection = 0;
                double uPredDiff_CentralCoeff = 0, uTdiff = 0;
                double vPredDiff_CentralCoeff = 0, vTdiff = 0;

				// Convection Term
                for (int k = 0; k < 4; k++) {
                    // Linear Upwind for Convection
                    if (flux[i][k] >= 0) {
                        uu[k] = u_pred[i] + (u_gradX[i])*(fcX[i][k] - xc[i]) + (u_gradY[i])*(fcY[i][k] - yc[i]);
                        vv[k] = v_pred[i] + (v_gradX[i])*(fcX[i][k] - xc[i]) + (v_gradY[i])*(fcY[i][k] - yc[i]);
                        convection_CentralCoeff = convection_CentralCoeff + flux[i][k];
                    } else {
                        uu[k] = u_pred[fout[i][k]] + (u_gradX[fout[i][k]])*(fcX[i][k] - xc[fout[i][k]]) + (u_gradY[fout[i][k]])*(fcY[i][k] - yc[fout[i][k]]);
                        vv[k] = v_pred[fout[i][k]] + (v_gradX[fout[i][k]])*(fcX[i][k] - xc[fout[i][k]]) + (v_gradY[fout[i][k]])*(fcY[i][k] - yc[fout[i][k]]);
                    }
                    // Total Convection for each cell for U and V
                    uPredConvection = uPredConvection + flux[i][k] * uu[k];
                    vPredConvection = vPredConvection + flux[i][k] * vv[k];

                    // Volume Interpolation for Nodal Values of U and V
                    uc_pred[lnc[i][k]] = 0;
                    vc_pred[lnc[i][k]] = 0;
                    for (int z = 0; z < max_cellShare; z++) {
                        if (lcn[lnc[i][k]][z] > 0) {
                            uc_pred[lnc[i][k]] = uc_pred[lnc[i][k]] + pointWt[lnc[i][k]][z] * u_curr[lcn[lnc[i][k]][z]];
                            vc_pred[lnc[i][k]] = vc_pred[lnc[i][k]] + pointWt[lnc[i][k]][z] * v_curr[lcn[lnc[i][k]][z]];
                        }
                    }

                }
                
				// Diffusion Term
                for (int k = 0; k < 4; k++) {
                    udiff[k] = normCoff[i][k] * (u_pred[fout[i][k]] - u_pred[i])
                            + crossCoff[i][k] * snSign[i][k]* (uc_pred[faces[lcf[i][k]][1]] - uc_pred[faces[lcf[i][k]][0]]);
                    vdiff[k] = normCoff[i][k] * (v_pred[fout[i][k]] - v_pred[i])
                            + crossCoff[i][k] * snSign[i][k]* (vc_pred[faces[lcf[i][k]][1]] - vc_pred[faces[lcf[i][k]][0]]);
                    // Central Coefficient for U and V Predicted Velocity diffusion term
                    uPredDiff_CentralCoeff = uPredDiff_CentralCoeff + normCoff[i][k];
                    vPredDiff_CentralCoeff = vPredDiff_CentralCoeff + normCoff[i][k];
                    // Total Diffsion for U and V Predicted Velocity
                    uTdiff = uTdiff + udiff[k];
                    vTdiff = vTdiff + vdiff[k];
                }

                double uPredTotal_CentralCoeff, vPredTotal_CentralCoeff;
                // Total Central Coefficient for U and V Predicted Velocity
                uPredTotal_CentralCoeff = Vp[i] / dt + convection_CentralCoeff + nu * uPredDiff_CentralCoeff;
                vPredTotal_CentralCoeff = Vp[i] / dt + convection_CentralCoeff + nu * vPredDiff_CentralCoeff;

                // Calculating Error or Residue and Correcting the value using Gauss Seidal iterative method
                uPred_error 	= Vp[i] / dt * (u_curr[i] - u_pred[i]) - uPredConvection + nu * uTdiff;
                uPred_residual 	= uPred_residual + uPred_error * uPred_error;
                u_pred[i] 		= uPred_error / uPredTotal_CentralCoeff + u_pred[i];

                vPred_error 	= Vp[i] / dt * (v_curr[i] - v_pred[i]) - vPredConvection + nu * vTdiff;
                vPred_residual 	= vPred_residual + vPred_error * vPred_error;
                v_pred[i] 		= vPred_error / vPredTotal_CentralCoeff + v_pred[i];

            }
            residual_error = sqrt((uPred_residual + vPred_residual) / (num_cells));

            // BC for Predicted Velocity ( Must implement according to order in the mesh file )

            for (int j = num_intFaces + 1; j <= num_intFaces + num_bcFaces[0]; j++) {	// Inlet = Free Stream
                u_pred[faces[j][3]] = U_freeStream;
                v_pred[faces[j][3]] = 0;
            }
            for (int j = num_intFaces + num_bcFaces[0] + 1; j <= num_intFaces + num_bcFaces[0] + num_bcFaces[1]; j++) {
                u_pred[faces[j][3]] = u_pred[faces[j][2]];								// Outlet = Neumon or Zero Gradient
                v_pred[faces[j][3]] = v_pred[faces[j][2]];
            }
            for (int j = num_intFaces + num_bcFaces[0] + num_bcFaces[1] + 1; j <= num_intFaces + num_totBC; j++) {
                u_pred[faces[j][3]] = 0;												// Wall = No Slip
                v_pred[faces[j][3]] = 0;
            }

        }
        
        cout << "GaussSeidalSolver: Solving for Ux predicted, Final residual = " << sqrt(uPred_residual / num_cells) << "  No Iterations = " << iteration << endl;
        cout << "GaussSeidalSolver: Solving for Uy predicted, Final residual = " << sqrt(vPred_residual / num_cells) << "  No Iterations = " << iteration << endl;

        // Pressure Poisson 

        iteration = 0;
        residual_error = 1;
        double pressure_residual, p_Next_error;

        while (residual_error > 1e-03) {

            iteration = iteration + 1;
            pressure_residual = 0;

			
            for (int i = 1; i <= num_cells; i++) {

                double pNext_CentralCoeff = 0, pTdiff, totalflux_noPressureCorrection;
                pTdiff = 0.0;
                totalflux_noPressureCorrection = 0.0;

				// Pressure Nodal Values Calculation
                for (int k = 0; k < 4; k++) {

                    pc_Next[lnc[i][k]] = 0;
                    for (int z = 0; z < max_cellShare; z++) {
                        if (lcn[lnc[i][k]][z] > 0) {
                            // Volume Interpolation for Nodal Values of P
                            pc_Next[lnc[i][k]] = pc_Next[lnc[i][k]] + pointWt[lnc[i][k]][z] * p_Next[lcn[lnc[i][k]][z]];
                        }
                    }

                }
				
				// Diffusion Term
                for (int k = 0; k < 4; k++) {
                    pdiff[i][k] = normCoff[i][k] * (p_Next[fout[i][k]] - p_Next[i])
                            + crossCoff[i][k] * snSign[i][k]* (pc_Next[faces[lcf[i][k]][1]] - pc_Next[faces[lcf[i][k]][0]]);
                    // Central Coefficient for Pressure Diffusion
                    pNext_CentralCoeff = pNext_CentralCoeff + normCoff[i][k];
                    // Total Diffusion of presuure
                    pTdiff = pTdiff + pdiff[i][k];

                    // Calculating Flux without pressure correction for each face
                    flux_noPressureCorrection[k] = (u_pred[i] * faceWt[i][k] + u_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k]
                            + (v_pred[i] * faceWt[i][k] + v_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k];

                    totalflux_noPressureCorrection = totalflux_noPressureCorrection + flux_noPressureCorrection[k];
                }

				// Calculating Error or Residue and Correcting the value using Gauss Seidal iterative method
                p_Next_error = totalflux_noPressureCorrection / dt * rho - pTdiff;
                pressure_residual = pressure_residual + p_Next_error * p_Next_error;
                p_Next[i] = -1.9 * p_Next_error / pNext_CentralCoeff + p_Next[i];		// Can use Over or Under Relaxation
            }
            residual_error = sqrt(pressure_residual / (num_cells));

            // BC for Pressure 

            for (int j = num_intFaces + 1; j <= num_intFaces + num_bcFaces[0]; j++) {	// Inlet = Neumon
                p_Next[faces[j][3]] = p_Next[faces[j][2]];
            }
            for (int j = num_intFaces + num_bcFaces[0] + 1; j <= num_intFaces + num_bcFaces[0] + num_bcFaces[1]; j++) {
                p_Next[faces[j][3]] = 0;												// Outlet = Zero
            }
            for (int j = num_intFaces + num_bcFaces[0] + num_bcFaces[1] + 1; j <= num_intFaces + num_totBC; j++) {
                p_Next[faces[j][3]] = p_Next[faces[j][2]];								// Wall = Neumon
            }

        }

        cout << "GaussSeidalSolver: Solving for P, Final residual = " << sqrt(pressure_residual / num_cells) << "  No Iterations = " << iteration << endl;


        // Corrected Velocity

        iteration = 0;
        residual_error = 1;
        double uNext_residual, vNext_residual, uNext_error, vNext_error;
        double uNextTotal_CentralCoeff, vNextTotal_CentralCoeff;

        while (residual_error > 1e-09) {
            uNext_residual = 0;
            vNext_residual = 0;

            iteration = iteration + 1;

            for (int i = 1; i <= num_cells; i++) {			// Calculating Flux and Cell Gradients for Linear Upwind

                u_gradX[i] = 0.0;
                u_gradY[i] = 0.0;
                v_gradX[i] = 0.0;
                v_gradY[i] = 0.0;

                for (int k = 0; k < 4; k++) {

                    flux[i][k] = (u_pred[i] * faceWt[i][k] + u_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k]
                            + (v_pred[i] * faceWt[i][k] + v_pred[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k] - dt / rho * pdiff[i][k];

                    u_gradX[i] = u_gradX[i] + (u_Next[i] * faceWt[i][k] + u_Next[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k] / Vp[i];
                    u_gradY[i] = u_gradY[i] + (u_Next[i] * faceWt[i][k] + u_Next[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k] / Vp[i];
                    v_gradX[i] = v_gradX[i] + (v_Next[i] * faceWt[i][k] + v_Next[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k] / Vp[i];
                    v_gradY[i] = v_gradY[i] + (v_Next[i] * faceWt[i][k] + v_Next[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k] / Vp[i];

                }
            }

            for (int i = 1; i <= num_cells; i++) {

                double convection_CentralCoeff = 0;
                double uConvection = 0, vConvection = 0;
                double uNextDiff_CentralCoeff = 0, uTdiff = 0;
                double vNextDiff_CentralCoeff = 0, vTdiff = 0;
                double SigmaPfSfx = 0, SigmaPfSfy = 0;

				// Convection Term
                for (int k = 0; k < 4; k++) {
					// Linear Upwind for Convection
                    if (flux[i][k] >= 0) {
                        uu[k] = u_Next[i] + (u_gradX[i])*(fcX[i][k] - xc[i]) + (u_gradY[i])*(fcY[i][k] - yc[i]);
                        vv[k] = v_Next[i] + (v_gradX[i])*(fcX[i][k] - xc[i]) + (v_gradY[i])*(fcY[i][k] - yc[i]);
                        convection_CentralCoeff = convection_CentralCoeff + flux[i][k];
                    } else {
                        uu[k] = u_Next[fout[i][k]] + (u_gradX[fout[i][k]])*(fcX[i][k] - xc[fout[i][k]]) + (u_gradY[fout[i][k]])*(fcY[i][k] - yc[fout[i][k]]);
                        vv[k] = v_Next[fout[i][k]] + (v_gradX[fout[i][k]])*(fcX[i][k] - xc[fout[i][k]]) + (v_gradY[fout[i][k]])*(fcY[i][k] - yc[fout[i][k]]);
                    }

                    // Total Convection for each cell for U and V
                    uConvection = uConvection + flux[i][k] * uu[k];
                    vConvection = vConvection + flux[i][k] * vv[k];

                    uc_Next[lnc[i][k]] = 0;
                    vc_Next[lnc[i][k]] = 0;
                    for (int z = 0; z < max_cellShare; z++) {
                        if (lcn[lnc[i][k]][z] > 0) {
                            // Volume Interpolation for Nodal Values of U and V
                            uc_Next[lnc[i][k]] = uc_Next[lnc[i][k]] + pointWt[lnc[i][k]][z] * u_curr[lcn[lnc[i][k]][z]];
                            vc_Next[lnc[i][k]] = vc_Next[lnc[i][k]] + pointWt[lnc[i][k]][z] * v_curr[lcn[lnc[i][k]][z]];
                        }
                    }
                }

				// Diffusion Term and Pressure Term in Momentum for U and V
                for (int k = 0; k < 4; k++) {

                    udiff[k] = normCoff[i][k] * (u_Next[fout[i][k]] - u_Next[i])
                            + crossCoff[i][k] * snSign[i][k]* (uc_Next[faces[lcf[i][k]][1]] - uc_Next[faces[lcf[i][k]][0]]);
                    vdiff[k] = normCoff[i][k] * (v_Next[fout[i][k]] - v_Next[i])
                            + crossCoff[i][k] * snSign[i][k]* (vc_Next[faces[lcf[i][k]][1]] - vc_Next[faces[lcf[i][k]][0]]);
					// Central Coefficient for U and V Corrected Velocity diffusion term
                    uNextDiff_CentralCoeff = uNextDiff_CentralCoeff + normCoff[i][k];
                    vNextDiff_CentralCoeff = vNextDiff_CentralCoeff + normCoff[i][k];
					// Total Diffsion for U and V Corrected Velocity
                    uTdiff = uTdiff + udiff[k];
                    vTdiff = vTdiff + vdiff[k];

					// Pressure for U Velocity
                    PSx[k] = (p_Next[i] * faceWt[i][k] + p_Next[fout[i][k]]*(1 - faceWt[i][k])) * sn1[i][k];
                    SigmaPfSfx = SigmaPfSfx + PSx[k];
					// Pressure for V Velocity
                    PSy[k] = (p_Next[i] * faceWt[i][k] + p_Next[fout[i][k]]*(1 - faceWt[i][k])) * sn2[i][k];
                    SigmaPfSfy = SigmaPfSfy + PSy[k];
                }

                // Total Central Coefficient for U and V Corrected Velocity
                uNextTotal_CentralCoeff = Vp[i] / dt + convection_CentralCoeff + nu * uNextDiff_CentralCoeff;
                vNextTotal_CentralCoeff = Vp[i] / dt + convection_CentralCoeff + nu * vNextDiff_CentralCoeff;

                // Calculating Error or Residue and Correcting the value using Gauss Seidal iterative method
                uNext_error = Vp[i] / dt * (u_curr[i] - u_Next[i]) - uConvection + nu * uTdiff - SigmaPfSfx / rho;
                uNext_residual = uNext_residual + uNext_error * uNext_error;
                u_Next[i] = uNext_error / uNextTotal_CentralCoeff + u_Next[i];

                vNext_error = Vp[i] / dt * (v_curr[i] - v_Next[i]) - vConvection + nu * vTdiff - SigmaPfSfy / rho;
                vNext_residual = vNext_residual + vNext_error * vNext_error;
                v_Next[i] = vNext_error / vNextTotal_CentralCoeff + v_Next[i];

            }
            residual_error = sqrt((uNext_residual + vNext_residual) / (num_cells));

            // BC for Corrected velocities

            for (int j = num_intFaces + 1; j <= num_intFaces + num_bcFaces[0]; j++) {		// Inlet = Free Stream
                u_Next[faces[j][3]] = U_freeStream;
                v_Next[faces[j][3]] = 0;
            }
            for (int j = num_intFaces + num_bcFaces[0] + 1; j <= num_intFaces + num_bcFaces[0] + num_bcFaces[1]; j++) {
                u_Next[faces[j][3]] = u_Next[faces[j][2]];									// Outlet = Neumon or Zero Gradient
                v_Next[faces[j][3]] = v_Next[faces[j][2]];
            }
            for (int j = num_intFaces + num_bcFaces[0] + num_bcFaces[1] + 1; j <= num_intFaces + num_totBC; j++) {
                u_Next[faces[j][3]] = 0;													// Wall = No Slip
                v_Next[faces[j][3]] = 0;
            }

        }

        cout << "GaussSeidalSolver: Solving for Ux, Final residual = " << sqrt(uNext_residual / num_cells) << "  No Iterations = " << iteration << endl;
        cout << "GaussSeidalSolver: Solving for Uy, Final residual = " << sqrt(vNext_residual / num_cells) << "  No Iterations = " << iteration << endl;

        // calculating the error values
        error = 0;

        for (int i = 1; i <= num_cells + num_totBC; i++) {
            error = error + pow((u_Next[i] - u_curr[i]), 2) + pow((v_Next[i] - v_curr[i]), 2);
            u_curr[i] = u_Next[i];
            v_curr[i] = v_Next[i];
        }

        error = sqrt(error / (num_cells * dt));
        cout << " error = " << error << endl;

		// Exporting to Tecplot or ParaView
		
        if (remainder((iter - 1), 50) == 0.0) {
            sol_time = sol_time + 1;
            
            char fname[100];
            FILE* fid;
            sprintf(fname, "flowOverCyl_%d.plt", sol_time);
            fid = fopen(fname, "w");
            fprintf(fid, "VARIABLES= X,Y,U,V,P \n");
            fprintf(fid, "ZONE T=\"U\" , F=FEPOINT, N=");
            fprintf(fid, "%d, ", num_nodes);
            fprintf(fid, "E=%d, ET=QUADRILATERAL \n", num_cells);
            for (int i = 1; i <= num_nodes; i++) {
                fprintf(fid, " %lf %lf %lf %f %lf \n ", vertices[i][0], vertices[i][1], uc_Next[i], vc_Next[i], pc_Next[i]);
            }
            fprintf(fid, "\n\n");
            for (int i = 1; i <= num_cells; i++) {
                fprintf(fid, " %d %d %d %d \n ", cellOrder[i][0], cellOrder[i][1], cellOrder[i][2], cellOrder[i][3]);
            }
            fclose(fid);
        }
        
        }

        return 0;
    }
