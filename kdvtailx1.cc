//finite FKDV
//g++ kdvtailx1.cc -o kdvtailx1 -lcln -O3 -Wall -Wextra
#include <iostream>
#include <cln/cln.h>
#include <iomanip>
#include <cmath>
#include <fstream>
#include <chrono>

using namespace std;
using namespace cln;
using std::ofstream;

const int nn = 120;
const int digits = 40;

//Initial function
cl_F u0(cl_F x, cl_F c0){
    cl_F a   = c0/2;       
    cl_F gam = sqrt(c0)/2; 
    return a/cosh(gam*x)/cosh(gam*x);
}

//TRANSFORMATION MATRICES
cl_F T[nn][nn];   //(transforms from collocation to Fourier)
cl_F G[nn][nn];   //(transforms from Fourier to collocation)

//Second derivative matrix
cl_F ddmat[nn][nn];

//Fourth derivative matrix
cl_F ddddmat[nn][nn];

//Second derivative operator
cl_F B[nn][nn];

//Fourth derivative operator
cl_F D[nn][nn];

//Ordinary matrix
cl_F X[nn][nn];

// Whole operator matrix whcih operates on delta matrix
cl_F linop[nn][nn];

void linsol(cl_F aa[][nn], cl_F yy[], cl_F xx[]);
void timing(const string str);

//Radiation coefficient (alpha)
cl_F kdvtail(cl_F del);

//Outer boundary
//cl_F xm = cl_float(20);

//RADIATION COEFICIENT FUNCTION
cl_F kdvtail(cl_F xm){
    
    default_print_flags.default_float_format = default_float_format;
    cl_F ppi = pi();
    
    cout << setprecision(15);
    
    cl_F eps = cl_float(1)/100;
    //cl_F del = cl_float(0);
    cl_F del = ppi/2;
    cl_F c0  = cl_float(1);
    cl_F c   = c0 * (cl_float(1) + eps*eps*c0);
    cl_F kk  = sqrt(cl_float(1) + eps*eps*c0); 
    
    //u at collocation points
    cl_F u0coll[nn];
    
    //Second derivative collocation points of u0
    cl_F ddu0coll[nn];
        
    for(int i = 0; i < nn; i++){
        u0coll[i] = u0(xm*cos(ppi*cl_float(i)/(nn-1)/2), c0);
    }  
    
    //Matrix T (transforms function from collocation to Fourier)
    for(int n = 0; n < nn; n++){
        for(int k = 0; k < nn; k++){
            cl_F npk = cl_float((n*k)%(2*(nn-1)));
            if(k == 0 || k == nn-1)
                T[n][k] = cos(npk*ppi/(nn-1))/((nn-1));
            else
                T[n][k] = cl_float(2)*cos(npk*ppi/((nn-1)))/(nn-1);
        }
    }

    // Matrix G (which is inverse of T and transforms the function from Fourier to collocation)
    for(int k = 0; k < nn; k++){
        for(int n = 0; n < nn; n++){
            cl_F npk = cl_float((n*k)%(2*(nn-1)));
            if(n == 0 || n == nn-1)
                G[k][n] = cos(npk*ppi/((nn-1)))/2;
            else
                G[k][n] = cos(npk*ppi/((nn-1)));
        }
    }

    cl_F tmp;
    //Second derivative 
    for(int n = 0; n < nn; n++)
        for(int k = 0; k < nn; k++){
            if(k >= n+1)
                ddmat[n][k] = cl_float(8*k*(k*k-n*n));
            else
                ddmat[n][k] = cl_float(0);
        }

    //Fourth derivative
    for(int i = 0; i < nn; i++){
        for(int j = 0; j < nn; j++){
            tmp = cl_float(0);
            for(int k = 0; k < nn; k++){
                tmp += ddmat[i][k] * ddmat[k][j];
            }
        ddddmat[i][j] = tmp;
        }
    }

    //Divide the first and last column of all the derivative matrices by 2 (only last column as first column is all zero)
    for(int i = 0; i < nn; i++){
        ddmat[i][nn-1] /= cl_float(2);
        ddddmat[i][nn-1] /= cl_float(2);          
    }

    //cout << "derivative matrices in Fourier picture calculated" << endl; 

    //DERIVATIVE OPERATORS

    //Second derivative operator
    for(int i = 0; i < nn; i++){
        for(int j = 0; j < nn; j++){
            tmp = cl_float(0);
            for(int k = 0; k < nn; k++){
                tmp += ddmat[i][k] * T[k][j]/xm/xm;
            }
            X[i][j] = tmp;
        }
    }

    for(int i = 0; i < nn; i++){
        for(int j = 0; j < nn; j++){
            tmp = cl_float(0);
            for(int k = 0; k < nn; k++){
                tmp += G[i][k] * X[k][j];
            }
            B[i][j] = tmp;
        }
    }

    // cout << "second derivative matrix calculated" << endl; 

    //Fourth derivative operator
    for(int i = 0; i < nn; i++){
        for(int j = 0; j < nn; j++){
            tmp = cl_float(0);
            for(int k = 0; k < nn; k++){
                tmp += ddddmat[i][k] * T[k][j]/xm/xm/xm/xm;
            }
            X[i][j] = tmp;
        }
    }

    for(int i = 0; i < nn; i++){
        for(int j = 0; j < nn; j++){
            tmp = cl_float(0);
            for(int k = 0; k < nn; k++){
                tmp += G[i][k] * X[k][j];
            }
            D[i][j] = tmp;
        }
    }
    
    //LOOP FOR THE WHOLE EQUATION
    for(int s = 0; s < 20; s++){
        
        //COLOCATION VALUES OF DERIVATIVES OF u0.
    
        //Second derivative collocation points of u0
        //cl_F ddu0coll[nn];
        for(int i = 0; i < nn; i++){
            tmp = cl_float(0);
            for(int j = 0; j < nn; j++){
            tmp += B[i][j] * u0coll[j];
            }
            ddu0coll[i] = tmp;
        }

        //Fourth derivative collocation points of u0
        cl_F ddddu0coll[nn];
        for(int i = 0; i < nn; i++){
            tmp = cl_float(0);
            for(int j = 0; j < nn; j++){
                tmp += D[i][j] * u0coll[j];
            }
            ddddu0coll[i] = tmp;
        }
        
        // Residual matrix
        cl_F Res[nn];         // Residual matrix
        for(int k = 0; k < nn; k++){
            Res[k] = c*u0coll[k] - cl_float(3)*u0coll[k]*u0coll[k] - ddu0coll[k] - eps*eps*ddddu0coll[k];
        }
        
        // error in boundary condition b2 = u''+u*k^2/eps^2 = 0
        cl_F b2 = ddu0coll[0]*eps*eps + u0coll[0]*kk*kk;
        
        // error in boundary condition 
        // b1 = u'*sin(x*k/eps-del)*eps-u*cos(x*k/eps-del)*kk = 0
        cl_F d11v[nn];
        for(int k = 0; k < nn; k++){
            tmp = cl_float(0);
            for(int j = 0; j < nn; j++){
                tmp += cl_float(4*j*j)*T[j][k]/xm;
            }
            d11v[k] = tmp;
        }
        cl_F d11 = cl_float(0);
        for(int k = 0; k < nn; k++)
            d11 += d11v[k]*u0coll[k];
        cl_F b1 = d11*eps*sin(xm*kk/eps-del) 
            - u0coll[0]*kk*cos(xm*kk/eps-del);

        //cout << "boundary cond: " << double_approx(b1)
            //<< " , " << double_approx(b2) << endl;
        
        Res[0] = -b2;
        Res[1] = -b1;
        
        // Whole operator matrix whcih operates on delta matrix
        for(int i = 0; i < nn; i++){
            for(int j = 0; j < nn; j++)
                linop[i][j] =  B[i][j] +eps*eps*D[i][j];
        }
    
        for(int i = 0; i < nn; i++){
            linop[i][i] += -c + cl_float(6)*u0coll[i];
        }
    
        for(int i = 0; i < nn; i++)
            linop[0][i] = B[0][i]*eps*eps;
        linop[0][0] += kk*kk;
    
        for(int i = 0; i < nn; i++)
            linop[1][i] = d11v[i]*eps*sin(xm*kk/eps-del);
        linop[1][1] -= kk*cos(xm*kk/eps-del);
        
        //Calculation of delta matrix
        cl_F delv[nn];
        linsol(linop, Res, delv);
        for(int k = 0; k < nn; k++){
            //cout << delv[k] << " ";
        }
    
        // Next approximation [u(n) = u(n-1) + delta)]
        for(int i = 0; i < nn; i++){
            u0coll[i] = u0coll[i] + delv[i];
        }
        
        // r.m.s average of the residual
        cl_F rmsaver;
        tmp = cl_float(0);
        for(int i = 0; i < nn; i++){
            tmp += Res[i] * Res[i];
        }
        rmsaver = sqrt(tmp/nn);
        
        //finding maximum of residual 
        cl_F maxr = cl_float(0);
        for(int i = 0; i < nn; i++){
            cl_F temp = abs(Res[i]);
            if(temp > maxr)
                maxr = temp;
        }
        
        // cout << "Residual" << " " << double_approx(rmsaver) << " " << double_approx(maxr) << endl;
        
    }
        
    return u0coll[0]/sin(kk*xm/eps - del);
    
}

// --------------- main function ------------------- 
int main(){
    
    default_float_format = float_format(digits); 
    
    cl_F x;
    
    //outer boundary
    cl_F xm = cl_float(20);
    
    cl_F ppi = pi();
    //int go = 50;
    
    ofstream ofs("kdvtailx1.out");
    ofs.precision(20);
    for(int i = 9; i < nn; ++i){
        x = xm*cos(ppi*cl_float(i)/(nn-1)/2); 
        cl_F alpha = kdvtail(x);
        ofs << x << " " << alpha << endl;
    cout << x << "  " << alpha << endl;
    cout << "---------------------------------------" << endl;
    }
    ofs.close();
           
}
//------------------ end of main function --------------------


// calculating LU-decomposition of aa
// solving aa*xx=yy for xx
// changing both aa and yy, returning result in xx
void linsol(cl_F aa[][nn], cl_F yy[], cl_F xx[]){

    //const cl_F TINY=1.0e-40;
    cl_F vv[nn];
    int pk[nn];
    cl_F big, temp;
    for (int i=0;i<nn;i++) {
        big=cl_float(0);
        for (int j=0;j<nn;j++){
            temp = aa[i][j];
            if (temp < cl_float(0))
                temp *= cl_float(-1);
            if (temp > big) big=temp;
        }
        //for (int j=0;j<nn;j++)
        //    if ((temp=std::abs(aa[i][j])) > big) big=temp;
        if (big == cl_float(0)) throw std::runtime_error("Singular matrix in LUdcmp");
        vv[i]=cl_float(1)/big;
    }
    for (int k=0;k<nn;k++) {
        int imax=k;
        big=cl_float(0);
        for (int i=k;i<nn;i++) {
            temp = aa[i][k];
            if (temp < cl_float(0))
                temp *= cl_float(-1);
            temp *= vv[i];
            //temp=vv[i]*std::abs(aa[i][k]);
            if (temp > big) {
                big=temp;
                imax=i;
            }
        }

        if (k != imax) {
            for (int j=0;j<nn;j++) {
                temp=aa[imax][j];
                aa[imax][j]=aa[k][j];
                aa[k][j]=temp;
            }
            vv[imax]=vv[k];
        }
        pk[k]=imax;
        //if (aa[k][k] == 0) aa[k][k]=TINY;
        for (int i=k+1;i<nn;i++) {
            temp=aa[i][k] /= aa[k][k];
            for (int j=k+1;j<nn;j++)
                aa[i][j] -= temp*aa[k][j];
        }
    }

    cl_F tt;
    for (int i=0; i<nn; ++i)
        if (pk[i] != i){
            tt = yy[i];
            yy[i] = yy[pk[i]];
            yy[pk[i]] = tt;
        }

    cl_F ww[nn], s;
    ww[0] = yy[0];
    for (int i = 1; i != nn; ++i){
        s = cl_float(0);
        for (int j = 0; j <= i - 1; ++j)
            s += aa[i][j]*ww[j];
        ww[i] = yy[i] - s;
    }

    xx[nn-1] = ww[nn-1]/aa[nn-1][nn-1];
    for (int i = nn-2; i >= 0; --i){
        s = cl_float(0);
        for (int j = i + 1; j != nn; ++j)
            s += aa[i][j]*xx[j];
        xx[i] = (ww[i] - s)/aa[i][i];
    }

}

void timing(const string str){
    static int first = 1;
    static auto t1 = std::chrono::high_resolution_clock::now();
    if (first != 1){
        auto t2 = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(t2 - t1);
        cout << str << ": " << duration.count() << " ms\n";
        t1 = t2;
    }
    first = 0;
}

