#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
    // Διαστάσεις πλέγματος
    const int Nx =121;
    const int Ny = 121;
    const int Nt = 80000;
    const double dt = 0.02;
    const double Lx = 10.0;
    const double Ly = 8.0;
    const double dx = 1.0 / (Nx - 1);
    const double dy = (Ly/Lx) / (Ny - 1);
    const double U_inf =  3.5;
    const double visc = 0.07;
    const double Re = U_inf*Lx/visc;
    
    
    // Συντελεστές !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //const double ap = 2/(dx * dx) + 2/(dy * dy);
   // const double ae = 1/(dx * dx);
    //const double aw = 1/(dx * dx);
   // const double as = 1/(dy * dy);
    //const double an = 1/(dy * dy);

    double psi[Ny][Nx] = {0.0};
    double omega[Ny][Nx] = {0.0};
    vector<vector<double>> u(Ny, vector<double>(Nx, 0.0));
    vector<vector<double>> v(Ny, vector<double>(Nx, 0.0));

    

    const double tol = 1e-6;      // κριτήριο σύγκλισης
    const int maxIter = 10;    // μέγιστος αριθμός επαναλήψεων


    for (int n = 0; n < Nt; ++n) {

        //ADI

        for (int iter = 0; iter < maxIter; ++iter){

            //x-σαρωση

            for (int j = 1; j <= Ny-2; ++j) {

                

            // ---- 2) Στήσιμο τριδιαγώνιου συστήματος για ΟΛΟΥΣ τους i = 0..Nx-1
            //     ώστε στα άκρα να προκύπτει ταυτοτικά psi=0
            double aw_i[Nx], ap_i[Nx], ae_i[Nx], rhs[Nx], as_i[Nx], an_i[Nx];

           

            // i = 0 (αριστερό σύνορο):  1*psi[0] = 0
            aw_i[0] = 0.0;
            ap_i[0] = 1.0;
            ae_i[0] = 0.0;
            rhs[0]  = 0.0;

            // i = 1..Nx-2 (εσωτερικοί κόμβοι)
            for (int i = 1; i <= Nx-2; ++i) {
                aw_i[i] = -1.0/(dx*dx);
                ap_i[i] =  2.0/(dx*dx) + 2.0/(dy*dy);
                ae_i[i] = -1.0/(dx*dx);
                as_i[i] = 1.0/(dy*dy);
                an_i[i] = 1.0/(dy*dy);

                // RHS από πάνω/κάτω + omega (όπως είχες)
                rhs[i]  = as_i[i]*psi[j-1][i] + an_i[i]*psi[j+1][i] + omega[j][i];
            }

            // i = Nx-1 (δεξί σύνορο):  1*psi[Nx-1] = 0
            aw_i[Nx-1] = 0.0;
            ap_i[Nx-1] = 1.0;
            ae_i[Nx-1] = 0.0;
            rhs[Nx-1]  = 0.0;

            // ---- 3) Thomas (TDMA) 
            double cprime[Nx];
            double dprime[Nx];

            // Forward sweep
            // i = 0
            {
                double denom = ap_i[0];
                cprime[0] = (denom != 0.0) ? (ae_i[0] / denom) : 0.0;
                dprime[0] = (denom != 0.0) ? (rhs[0] / denom) : 0.0;
            }

            // i = 1..Nx-1
            for (int i = 1; i < Nx; ++i) {
                double denom = ap_i[i] - aw_i[i] * cprime[i-1];
                cprime[i] = (denom != 0.0) ? (ae_i[i] / denom) : 0.0;
                dprime[i] = (denom != 0.0) ? ((rhs[i] - aw_i[i] * dprime[i-1]) / denom) : 0.0;
            }

            // Back substitution
            psi[j][Nx-1] = dprime[Nx-1];
            for (int i = Nx-2; i >= 0; --i) {
                psi[j][i] = dprime[i] - cprime[i] * psi[j][i+1];
            }

            
            }

            // -------------------------------------------
            //  ΣΑΡΩΣΗ ΚΑΤΑ y (δεύτερο βήμα ADI)
            // -------------------------------------------

            for (int i = 1; i <= Nx-2; ++i) {

            // Συντελεστές του συστήματος κατά y
            double as_j[Ny], ap_j[Ny], an_j[Ny],ae_j[Ny],aw_j[Ny], rhs_j[Ny];
            double cprime[Ny], dprime[Ny];

            // Ορισμός συντελεστών 
            const double aw = 1.0/(dx*dx);
            const double ae = 1.0/(dx*dx);
            const double as = -1.0/(dy*dy);
            const double an = -1.0/(dy*dy);
            const double ap = 2.0/(dx*dx) + 2.0/(dy*dy);

            // --- Οριακές συνθήκες στα j = 0 και j = Ny-1 ---
            as_j[0] = 0.0;   ap_j[0] = 1.0;   an_j[0] = 0.0;   rhs_j[0] = 0.0;
            as_j[Ny-1] = 0.0; ap_j[Ny-1] = 1.0; an_j[Ny-1] = 0.0; rhs_j[Ny-1] = 0.0;

            // --- Εσωτερικές εξισώσεις j = 1..Ny-2 ---
            for (int j = 1; j <= Ny-2; ++j) {
                as_j[j] = -1.0/(dy*dy);
                ap_j[j] =  2.0/(dx*dx) + 2.0/(dy*dy);
                an_j[j] = -1.0/(dy*dy);
                aw_j[j] = 1.0/(dx*dx);
                ae_j[j] = 1.0/(dx*dx);

                // RHS 
                rhs_j[j] = aw_j[j]*psi[j][i-1] + ae_j[j]*psi[j][i+1] + omega[j][i];
            }

            // -------------------------------------------
            //  TDMA κατά y για τη στήλη i
            // -------------------------------------------

            // Forward sweep
            double denom = ap_j[0];
            cprime[0] = (denom != 0.0) ? (an_j[0] / denom) : 0.0;
            dprime[0] = (denom != 0.0) ? (rhs_j[0] / denom) : 0.0;

            for (int j = 1; j < Ny; ++j) {
                denom = ap_j[j] - as_j[j] * cprime[j-1];
                cprime[j] = (denom != 0.0) ? (an_j[j] / denom) : 0.0;
                dprime[j] = (denom != 0.0) ? ((rhs_j[j] - as_j[j]*dprime[j-1]) / denom) : 0.0;
            }

            // Back substitution
            psi[Ny-1][i] = dprime[Ny-1];
            for (int j = Ny-2; j >= 0; --j) {
                psi[j][i] = dprime[j] - cprime[j] * psi[j+1][i];
            }
            }

            double max_r = 0.0;
            for (int j = 1; j <= Ny-2; ++j) {
                for (int i = 1; i <= Nx-2; ++i) {
                    double lap =
                        (psi[j][i+1] - 2.0*psi[j][i] + psi[j][i-1])/(dx*dx) +
                        (psi[j+1][i] - 2.0*psi[j][i] + psi[j-1][i])/(dy*dy);
                    double r = lap + omega[j][i];
                    if (fabs(r) > max_r) max_r = fabs(r);
                }
            }

            if (max_r < tol) {
        
                break;
            }

            
        }

          

        // --- Υπολογισμός ταχυτήτων από Ψ ---
        

        // Εσωτερικό: κεντρικές διαφορές
        for (int j = 1; j < Ny-1; ++j) {
            for (int i = 1; i < Nx-1; ++i) {
                u[j][i] = (psi[j+1][i] - psi[j-1][i]) / (2.0 * dy);
                v[j][i] = -(psi[j][i+1] - psi[j][i-1]) / (2.0 * dx);
            }
        }

        // Συνοριακές συνθήκες ταχύτητας 
        // κάτω τοίχος
        for (int i = 0; i < Nx; ++i) { u[0][i] = 0.0; v[0][i] = 0.0; }
        // πάνω τοίχος (κινούμενη πλάκα)
        for (int i = 0; i < Nx; ++i) { u[Ny-1][i] = 1.0; v[Ny-1][i] = 0.0; }
        // αριστερός/δεξιός τοίχος
        for (int j = 0; j < Ny; ++j) { u[j][0] = 0.0; v[j][0] = 0.0; }
        for (int j = 0; j < Ny; ++j) { u[j][Nx-1] = 0.0; v[j][Nx-1] = 0.0; }
        

        

            


        

        // OMEGAAAAA
        
        double omega_old[Ny][Nx];
        double omega_new[Ny][Nx];

        
        for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
            omega_old[j][i] = omega[j][i];   // παλιό ω
            omega_new[j][i] = omega[j][i];   // ω_new = ω
        }



        

        
            for (int j = 1; j <= Ny-2; ++j) {

            double aw_i[Nx], ap_i[Nx], ae_i[Nx], ap_i_dot[Nx], as_i[Nx], an_i[Nx],rhs_i[Nx];
            
            aw_i[0] = 0.0;
            ap_i[0] = 1.0;
            ae_i[0] = 0.0;
            rhs_i[0]  = 2*(psi[j][0]-psi[j][1])/(dx*dx);

            // i = 1..Nx-2 (εσωτερικοί κόμβοι)
            for (int i = 1; i <= Nx-2; ++i) {
                aw_i[i] = ((-u[j][i])/(2*dx))-(1/(Re*dx*dx));
                ap_i[i] =  (2/dt)+(2/(Re*dx*dx));
                ae_i[i] = ((u[j][i])/(2*dx))-(1/(Re*dx*dx));
                as_i[i] = (v[j][i]/(2*dy))+(1/(Re*dy*dy));
                an_i[i] = (-v[j][i]/(2*dy))+(1/(Re*dy*dy));
                ap_i_dot[i] = (2/dt)-(2/(Re*dy*dy));

                // RHS 
                rhs_i[i]  = as_i[i]*omega_new[j-1][i] + an_i[i]*omega_new[j+1][i] + ap_i_dot[i]*omega_new[j][i];
            }

            // i = Nx-1 (δεξί σύνορο):  1*psi[Nx-1] = 0
            aw_i[Nx-1] = 0.0;
            ap_i[Nx-1] = 1.0;
            ae_i[Nx-1] = 0.0;
            rhs_i[Nx-1]  =2*(psi[j][Nx-1]-psi[j][Nx-2])/(dx*dx);


            // THOMAS
            
            double cprime[Nx];
            double dprime[Nx];

            // forward sweep
            double denom = ap_i[0];
            cprime[0] = (denom != 0.0) ? (ae_i[0] / denom) : 0.0;
            dprime[0] = (denom != 0.0) ? (rhs_i[0] / denom) : 0.0;

            for (int i = 1; i < Nx; ++i) {
                denom = ap_i[i] - aw_i[i]*cprime[i-1];
                cprime[i] = (denom != 0.0) ? (ae_i[i] / denom) : 0.0;
                dprime[i] = (denom != 0.0)
                        ? ((rhs_i[i] - aw_i[i]*dprime[i-1]) / denom)
                        : 0.0;
            }

            // back substitution 
            omega_new[j][Nx-1] = dprime[Nx-1];
            for (int i = Nx-2; i >= 0; --i) {
                omega_new[j][i] = dprime[i] - cprime[i]*omega_new[j][i+1];
            }
            }

            // apo t+dt/2 -------> t+dt

            for (int i = 1; i <= Nx-2; ++i) {

            double aw_j[Ny], ap_j[Ny], ae_j[Ny], ap_j_dot[Ny], as_j[Ny], an_j[Ny],rhs_j[Ny];

            an_j[0] = 0.0;
            ap_j[0] = 1.0;
            as_j[0] = 0.0;
            rhs_j[0]  = 2*(psi[0][i]-psi[1][i])/(dy*dy);

            for (int j = 1; j <= Ny-2; ++j) {
                as_j[j] = (-v[j][i])/(2*dy) - 1/(Re*dy*dy);
                ap_j[j] = (2/dt)+(2/(Re*dy*dy));
                an_j[j] = (v[j][i])/(2*dy) - 1/(Re*dy*dy);
                aw_j[j] = (u[j][i])/(2*dx) + 1/(Re*dx*dx);
                ae_j[j] = (-u[j][i])/(2*dx) + 1/(Re*dx*dx);
                ap_j_dot[j] = (2/dt)-(2/(Re*dx*dx));

                // RHS 
                rhs_j[j] = aw_j[j]*omega_new[j][i-1] + ae_j[j]*omega_new[j][i+1] + ap_j_dot[j]* omega_new[j][i];
            }

            an_j[Ny-1] = 0;
            ap_j[Ny-1] = 1;
            as_j[Ny-1] = 0;
            rhs_j[Ny-1] = 2* (psi[Ny-1][i] - psi[Ny-2][i]-dy)/(dy*dy);

            double cprime[Ny];
            double dprime[Ny];

            // forward sweep
            double denom = ap_j[0];
            cprime[0] = (denom != 0.0) ? (an_j[0] / denom) : 0.0;
            dprime[0] = (denom != 0.0) ? (rhs_j[0] / denom) : 0.0;

            for (int j = 1; j < Ny; ++j) {
                denom = ap_j[j] - as_j[j] * cprime[j-1];
                cprime[j] = (denom != 0.0) ? (an_j[j] / denom) : 0.0;
                dprime[j] = (denom != 0.0)
                        ? ((rhs_j[j] - as_j[j]*dprime[j-1]) / denom)
                        : 0.0;
            }

            // back substitution 
            omega_new[Ny-1][i] = dprime[Ny-1];
            for (int j = Ny-2; j >= 0; --j) {
                omega_new[j][i] = dprime[j] - cprime[j] * omega_new[j+1][i];
            }

            }


             // ΚΡΙΤΗΡΙΟ ΣΥΓΚΛΙΣΗΣ
            //-------------------------
            double max_ratio = 0.0;

            for (int j = 0; j < Ny; ++j)
            for (int i = 0; i < Nx; ++i) {
                double ratio = fabs(omega_new[j][i] - omega_old[j][i]) / dt;
                if (ratio > max_ratio)
                    max_ratio = ratio;
     
            }

            cout << "[iter " << n << "] dw = " << max_ratio << endl;


            if (max_ratio < 1e-8) {
            
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i)
                    omega[j][i] = omega_new[j][i];

            

            break;
            }

           
            for (int j = 0; j < Ny; ++j)
                for (int i = 0; i < Nx; ++i) {
                    omega[j][i] = omega_new[j][i];
                }

        

    }   

// ofstream fpsi("psi_final121x121c.dat");

// for (int j = 0; j <=Ny-1; ++j) {
//     for (int i = 0; i < Nx; ++i)
//         fpsi << psi[j][i] << " ";
//     fpsi << "\n";
// }

// fpsi.close();

ofstream fomega("omega_finalgiorgos.dat");

for (int j = 0; j <=Ny-1; ++j) {
    for (int i = 0; i < Nx; ++i)
        fomega << omega[j][i] << " ";
    fomega << "\n";
}

fomega.close();

    return 0;
}

// ofstream fu("u_final121x121c.dat");

// for (int j = 0; j <=Ny-1; ++j) {
//     for (int i = 0; i < Nx; ++i)
//         fu << u[j][i] << " ";
//     fu << "\n";
// }

// fu.close();

// ofstream fv("v_final121x121c.dat");

// for (int j = 0; j <=Ny-1; ++j) {
//     for (int i = 0; i < Nx; ++i)
//         fv << v[j][i] << " ";
//     fv << "\n";
// }

// fv.close(); 





// // poisson gia PIESH

// // --- Poisson για πίεση ---
// // P ξεκινάει 0 παντού
//  double P[Ny][Nx] = {0.0};
// double sourceP[Ny][Nx] = {0.0};

// // Υπολογισμός όρου πηγής από το τελικό psi
// for (int j = 1; j < Ny-1; ++j) {
//     for (int i = 1; i < Nx-1; ++i) {
//         double psi_xx =
//             (psi[j][i+1] - 2.0*psi[j][i] + psi[j][i-1]) / (dx*dx);
//         double psi_yy =
//             (psi[j+1][i] - 2.0*psi[j][i] + psi[j-1][i]) / (dy*dy);
//         double psi_xy =
//             (psi[j+1][i+1] - psi[j+1][i-1]
//            - psi[j-1][i+1] + psi[j-1][i-1]) / (4.0*dx*dy);

//         sourceP[j][i] = 2.0 * (psi_xx * psi_yy - psi_xy * psi_xy);
//     }
// }

// // στα σύνορα source=0
// for (int i = 0; i < Nx; ++i) {
//     sourceP[0][i]    = 0.0;
//     sourceP[Ny-1][i] = 0.0;
// }
// for (int j = 0; j < Ny; ++j) {
//     sourceP[j][0]    = 0.0;
//     sourceP[j][Nx-1] = 0.0;
// }




// // ADIIIIIII

// const double tolP = 1e-8;
// const int    maxIterP = 10000;

// for (int iter = 0; iter < maxIterP; ++iter) {

//     // -------------------------
//     // 1ο βήμα ADI: σάρωση κατά x
//     // -------------------------
//     for (int j = 1; j <= Ny-2; ++j) {

//         double aw_i[Nx], ap_i[Nx], ae_i[Nx], rhs[Nx];
//         double cprime[Nx], dprime[Nx];

//         // αριστερό σύνορο: δP = 0
//         aw_i[0] = 0.0;
//         ap_i[0] = 1.0;
//         ae_i[0] = -1.0;
//         rhs[0]  = 0.0;

//         // εσωτερικοί κόμβοι
//         for (int i = 1; i <= Nx-2; ++i) {
//             double as = 1.0/(dy*dy);
//             double an = 1.0/(dy*dy);

//             aw_i[i] = -1.0/(dx*dx);
//             ap_i[i] =  2.0/(dx*dx) + 2.0/(dy*dy);
//             ae_i[i] = -1.0/(dx*dx);

//             rhs[i]  = as*P[j-1][i] + an*P[j+1][i] + sourceP[j][i];
//         }

//         // δεξί σύνορο: δP = 0
//         aw_i[Nx-1] = -1.0;
//         ap_i[Nx-1] = 1.0;
//         ae_i[Nx-1] = 0.0;
//         rhs[Nx-1]  = 0.0;

//         // Thomas κατά x
//         double denom = ap_i[0];
//         cprime[0] = (denom != 0.0) ? (ae_i[0] / denom) : 0.0;
//         dprime[0] = (denom != 0.0) ? (rhs[0] / denom) : 0.0;

//         for (int i = 1; i < Nx; ++i) {
//             denom = ap_i[i] - aw_i[i]*cprime[i-1];
//             cprime[i] = (denom != 0.0) ? (ae_i[i] / denom) : 0.0;
//             dprime[i] = (denom != 0.0)
//                       ? ((rhs[i] - aw_i[i]*dprime[i-1]) / denom)
//                       : 0.0;
//         }

//         P[j][Nx-1] = dprime[Nx-1];
//         for (int i = Nx-2; i >= 0; --i) {
//             P[j][i] = dprime[i] - cprime[i]*P[j][i+1];
//         }
//     }




//     // -------------------------
//     // 2ο βήμα ADI: σάρωση κατά y
//     // -------------------------
//     for (int i = 1; i <= Nx-2; ++i) {

//         double as_j[Ny], ap_j[Ny], an_j[Ny], rhs_j[Ny];
//         double cprime[Ny], dprime[Ny];

//         // κάτω σύνορο: δP = 0
//         as_j[0] = 0.0;
//         ap_j[0] = 1.0;
//         an_j[0] = -1.0;
//         rhs_j[0] = 0.0;

//         // εσωτερικοί κατά y
//         for (int j = 1; j <= Ny-2; ++j) {
//             double aw = 1.0/(dx*dx);
//             double ae = 1.0/(dx*dx);

//             as_j[j] = -1.0/(dy*dy);
//             ap_j[j] =  2.0/(dx*dx) + 2.0/(dy*dy);
//             an_j[j] = -1.0/(dy*dy);

//             rhs_j[j] = aw*P[j][i-1] + ae*P[j][i+1] + sourceP[j][i];
//         }

//         // πάνω σύνορο: P = 0
//         as_j[Ny-1] = 0.0;
//         ap_j[Ny-1] = 1.0;
//         an_j[Ny-1] = 0.0;
//         rhs_j[Ny-1] = 0.0;

//         // Thomas κατά y
//         double denom = ap_j[0];
//         cprime[0] = (denom != 0.0) ? (an_j[0] / denom) : 0.0;
//         dprime[0] = (denom != 0.0) ? (rhs_j[0] / denom) : 0.0;

//         for (int j = 1; j < Ny; ++j) {
//             denom = ap_j[j] - as_j[j]*cprime[j-1];
//             cprime[j] = (denom != 0.0) ? (an_j[j] / denom) : 0.0;
//             dprime[j] = (denom != 0.0)
//                       ? ((rhs_j[j] - as_j[j]*dprime[j-1]) / denom)
//                       : 0.0;
//         }

//         P[Ny-1][i] = dprime[Ny-1];
//         for (int j = Ny-2; j >= 0; --j) {
//             P[j][i] = dprime[j] - cprime[j]*P[j+1][i];
//         }
//     }

//     // -------------------------
//     // Κριτήριο σύγκλισης
//     // -------------------------
//     double max_r = 0.0;

//     for (int j = 1; j <= Ny-2; ++j) {
//         for (int i = 1; i <= Nx-2; ++i) {
//             double lapP =
//                 (P[j][i+1] - 2.0*P[j][i] + P[j][i-1])/(dx*dx) +
//                 (P[j+1][i] - 2.0*P[j][i] + P[j-1][i])/(dy*dy);

//             double r = lapP - sourceP[j][i];   
//             if (fabs(r) > max_r) max_r = fabs(r);
//         }
//     }

//     if (max_r < tolP) {
//         cout << "Poisson πιεσης: συνεκλινε στο iter = " << iter
//              << ", residual = " << max_r << endl;
//         break;
//     }
// }

// ofstream fP("pressure_121x121c.dat");

// for (int j = 0; j <=Ny-1; ++j) {
//     for (int i = 0; i < Nx; ++i)
//         fP << P[j][i] << " ";
//     fP << "\n";
// }

// fP.close();  //


// // TRIBES


// vector<double> tau_bottom(Nx);

// for (int i = 0; i < Nx; ++i) {
//     double dudyn = (-3.0*u[0][i] + 4.0*u[1][i] - u[2][i]) / (2.0*dy);
//     tau_bottom[i] = (1.0/Re) * dudyn;  // αδιάστατη τριβή
// }

// vector<double> tau_top(Nx);

// for (int i = 0; i < Nx; ++i) {
//     double dudy = (3.0*u[Ny-1][i] - 4.0*u[Ny-2][i] + u[Ny-3][i]) / (2.0*dy);
//     tau_top[i] = (1.0/Re) * dudy;
// }

// vector<double> tau_left(Ny);

// for (int j = 0; j < Ny; ++j) {
//     double dvdx = (-3.0*v[j][0] + 4.0*v[j][1] - v[j][2]) / (2.0*dx);
//     tau_left[j] = (1.0/Re) * dvdx;
// }

// vector<double> tau_right(Ny);

// for (int j = 0; j < Ny; ++j) {
//     double dvdx = (3.0*v[j][Nx-1] - 4.0*v[j][Nx-2] + v[j][Nx-3]) / (2.0*dx);
//     tau_right[j] = (1.0/Re) * dvdx;
// }


// // --- Export wall shear stresses ---

// // Bottom wall τριβή
// ofstream fb("tau_bottom_c.txt");
// for (int i = 0; i < Nx; ++i)
//     fb << tau_bottom[i] << "\n";
// fb.close();

// // Top wall τριβή
// ofstream ft("tau_top_c.txt");
// for (int i = 0; i < Nx; ++i)
//     ft << tau_top[i] << "\n";
// ft.close();

// // Left wall τριβή
// ofstream fl("tau_left_c.txt");
// for (int j = 0; j < Ny; ++j)
//     fl << tau_left[j] << "\n";
// fl.close();

// // Right wall τριβή
// ofstream fr("tau_right_c.txt");
// for (int j = 0; j < Ny; ++j)
//     fr << tau_right[j] << "\n";
// fr.close();



//     return 0;
// }
