#include<iostream>
#include<complex.h>
#include<complex>
#include<math.h>
using namespace std;

complex<double>** MUL_MATRIX(complex<double> x, complex<double>** A, int size_A){
    complex<double>** result = new complex<double>*[size_A];
    for(int i = 0; i < size_A; i++){
        result[i] = new complex<double>[size_A];
    }
    for(int i = 0; i < size_A; i++){
        for(int j = 0; j < size_A; j++){
            result[i][j] = x * A[i][j];
        }
    }
    return result;
}

complex<double>** tensor_product(complex<double> **A, complex<double> **B, int size_A, int size_B){
    int size = size_A * size_B;
    complex<double>** result = new complex<double>* [size];
    for(int i = 0; i < size; i++){
        result[i] = new complex<double>[size];
    }
    complex<double>** t = new complex<double>*[size_B];
    for(int i = 0; i < size_B; i++){
        t[i] = new complex<double> [size_B];
    }
    for(int i = 0; i < size_A; i++){
        for(int j = 0; j < size_A; j++){
            t = MUL_MATRIX(A[i][j], B, size_B);
            for(int k = 0; k < size_B; k++){
                for(int l = 0; l < size_B; l++){
                    result[size_B*i+k][size_B*j+l] = t[k][l];
                    // cout << result[size_B*i+k][size_B*j+l] << "\t";
                }
            }
        }
        // cout << endl;
    }
    return result;
}
int exp_int(int a, int n){
    int rs = 1;
    for(int i = 0; i < n; i++){
        rs *= a;
    }
    return rs;
}

complex<double>** exp_product(complex<double> **A, int size_A, int n){
    complex<double>** rs;
    rs = tensor_product(A, A, size_A, size_A);
    cout << endl;
    for(int i = 0; i < n-2; i++){
        int size = exp_int(size_A, i+2);
        int size_rs = exp_int(size_A, i+3);
        rs = tensor_product(rs, A, size, size_A);       
    }
    return rs;
}

complex<double>** Mul(complex<double> **A, complex<double> **B, int size){
    complex<double>** result = new complex<double>* [size];
    for(int i = 0; i < size; i++){
        result[i] = new complex<double> [size];
    }
    complex<double> t (0,0);
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            result[i][j] = t;
        }
    }
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            for(int k = 0; k< size; k++){
                result[i][j] += (A[i][k] * B[k][j]);
            }
        }
    }
    return result;
}

complex<double>* Mul_vector_matrix(complex<double>** A, complex<double>* b, int size) {
    complex<double>* result = new complex<double> [size];
    complex<double> t (0,0);
    for(int i = 0; i < size; i++){
        result[i] = t;
    }
    for(int i = 0; i < size; i++){
        for(int j = 0; j < size; j++){
            result[i] += (A[i][j] * b[j]);
        }
    }
    return result;
}

complex<double>* tensor_vector(complex<double>* a, complex<double>* b, int size_a, int size_b){
    int size = size_a * size_b;
    complex<double>* result =new complex<double> [size];
    complex<double> t (0,0);
    for(int i = 0; i < size; i++){
        result[i] = t;
    }
    for(int i = 0; i < size_a; i++){
        for(int j = 0; j < size_b; j++){
            result[size_b * i + j] += (a[i] * b[j]); 
        } 
    }
    return result;
}

complex<double>* gen_qubit(bitset<2>* a, int n){
    int size = pow(2,n);
    complex<double>** b = new complex<double> *[n];
    for(int i = 0; i < n; i++){
        b[i] = new complex<double> [2];
    }
    for(int i = 0; i < n;i++){
        if(a[i] == 0) a[i] = 2; 
        for(int j = 0; j < 2; j++){
            b[i][j] = a[i][(j+1)%2];
        }
    }
    // for(int i = 0; i < n; i++){
    //     for(int j = 0; j < 2; j++){
    //         cout << b[i][j].real() << " ";
    //     }
    //     cout << endl;
    // }
    // cout << endl;
    complex<double>* rs;
    rs = b[0];
    for(int i = 1; i < n; i++){
        rs = tensor_vector(rs, b[i], pow(2, i), 2);
        // for(int j = 0; j < pow(2,i+1); j++){
        //     cout << rs[j].real() << endl;
        // }
        // cout << endl;
    }
    return rs;
}

void print_Matrix(complex<double>** A, int size_A){
    for(int i = 0; i < size_A; i++){
        for(int j = 0; j < size_A; j++){
            if(A[i][j].real() == 0 && A[i][j].imag() == 0) {
                if(A[i][j].imag() == 0) cout << abs(A[i][j].real()) << "\t";
                else cout << abs(A[i][j]) << "\t";
            }
            else{
                if(A[i][j].imag() == 0) cout << A[i][j].real() << "\t";
                else cout << A[i][j]<< "\t";
            }
            
        }
        cout << endl;
    }
}

int main(){
    complex<double> h (1, 0);
    complex<double> h_1 (-1, 0);
    complex<double> **H = new complex<double>*[2];
    for(int i = 0; i < 2; i++){
        H[i] = new complex<double> [2];
    }
    H[0][0] = h;
    H[0][1] = h;
    H[1][0] = h;
    H[1][1] = h_1;
    H = MUL_MATRIX(1/sqrt(2), H, 2);

    complex<double> **I = new complex<double>*[2];
    for(int i = 0; i < 2; i++){
        I[i] = new complex<double> [2];
    }
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            if (i == j) {
                I[i][j] = 1;
            }
            else {
                I[i][j] = 0;
            }
        }
    }
    complex<double>** SWAP = new complex<double>* [4];
    for(int i = 0; i < 4; i++){
        SWAP[i] = new complex<double> [4];
    }
    SWAP[0][0] = 1;
    SWAP[0][1] = 0;
    SWAP[0][2] = 0;
    SWAP[0][3] = 0;
    SWAP[1][0] = 0;
    SWAP[1][1] = 0;
    SWAP[1][2] = 1;
    SWAP[1][3] = 0;
    SWAP[2][0] = 0;
    SWAP[2][1] = 1;
    SWAP[2][2] = 0;
    SWAP[2][3] = 0;
    SWAP[3][0] = 0;
    SWAP[3][1] = 0;
    SWAP[3][2] = 0;
    SWAP[3][3] = 1;

    complex<double>** CZ = new complex<double>* [4];
    for(int i = 0; i < 4; i++){
        CZ[i] = new complex<double> [4];
    }
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if (i == j && i != 3) CZ[i][j] = 1;
            else if (i == j && i == 3) CZ[i][j] = -1;
            else CZ[i][j] = 0;
        }
    }

    complex<double>** CNOT = new complex<double>* [4];
    for(int i = 0; i < 4; i++){
        CNOT[i] = new complex<double> [4];
    }
    for(int i = 0; i < 4;i++){
        for(int j = 0; j < 4; j++){
            if(i == j) CNOT[i][j] = 1;
            else CNOT[i][j] = 0;
        }
    }
    CNOT[2][2] = 0;
    CNOT[2][3] = 1;
    CNOT[3][2] = 1;
    CNOT[3][3] = 0;

    complex<double>** CCNOT = new complex<double>* [8];
    for(int i = 0; i < 8; i++){
        CCNOT[i] = new complex<double> [8];
    }
    for(int i = 0; i < 8;i++){
        for(int j = 0; j < 8; j++){
            if(i == j) CCNOT[i][j] = 1;
            else CCNOT[i][j] = 0;
        }
    }
    CCNOT[6][6] = 0;
    CCNOT[6][7] = 1;
    CCNOT[7][7] = 0;
    CCNOT[7][6] = 1;

    complex<double> var (1/sqrt(2), 1/sqrt(2));
    complex<double> **T = new complex<double> * [2];
    for(int i = 0; i < 2; i++){
        T[i] = new complex<double> [2];
    }
    T[0][0] = 1;
    T[0][1] = 0;
    T[1][0] = 0;
    T[1][1] = var;

    complex<double> var_1 (1/sqrt(2), -1/sqrt(2));
    complex<double> **T_1 = new complex<double> * [2];
    for(int i = 0; i < 2; i++){
        T_1[i] = new complex<double> [2];
    }
    T_1[0][0] = 1;
    T_1[0][1] = 0;
    T_1[1][0] = 0;
    T_1[1][1] = var_1;

    complex<double> **CT = new complex<double> * [4];
    for(int i = 0; i < 4; i++){
        CT[i] = new complex<double> [4];
    }
    for(int i = 0; i < 4; i++){
        for(int j = 0; j < 4; j++){
            if(i ==j) CT[i][j] = 1;
            else CT[i][j] = 0;
        }
    }
    CT[3][3] = var;

    complex<double> **S = new complex<double> * [2];
    for(int i = 0; i < 2; i++){
        S[i] = new complex<double> [2];
    }
    for(int i = 0; i < 2; i++){
        for(int j = 0; j < 2; j++){
            if(i == j) S[i][j] = 1;
            else S[i][j] = 0;
        }
    }
    complex<double> var_2 (0,1);
    S[1][1] = var_2;

    complex<double>** CCNOTI;
    complex<double>** CTT_1I;
    complex<double>** CCNOTH;
    complex<double>** IITT;
    complex<double>** IISWAP;
    complex<double>** IIIT_1;
    complex<double>** ICCNOT;
    complex<double>** IIIT;
    complex<double>** IIIH;
    complex<double>** CCCNOT;

    CCNOTI = tensor_product(CCNOT, I, 8, 2);
    CTT_1I = tensor_product(tensor_product(CT, T_1, 4, 2), I, 8, 2);
    CCNOTH = tensor_product(CCNOT, H, 8, 2);
    IITT = tensor_product(tensor_product(exp_product(I, 2, 2), T, 4, 2), T, 8, 2);
    IISWAP = tensor_product(exp_product(I, 2, 2), SWAP, 4, 4);
    IIIT_1 = tensor_product(exp_product(I, 2, 3), T_1, 8, 2);
    ICCNOT = tensor_product(I, CCNOT, 2, 8);
    IIIT = tensor_product(exp_product(I, 2, 3), T, 8, 2);
    IIIH = tensor_product(exp_product(I, 2, 3), H, 8 ,2);
    
    CCCNOT = Mul(Mul(Mul(CCNOTI, CTT_1I, 16), CCNOTH, 16), IITT, 16);
    CCCNOT = Mul(Mul(Mul(CCCNOT, IISWAP, 16), CCNOTI, 16), IISWAP, 16);
    CCCNOT = Mul(Mul(Mul(CCCNOT, IIIT_1, 16), ICCNOT, 16), IIIT, 16);
    CCCNOT = Mul(Mul(Mul(CCCNOT, IISWAP, 16), CCNOTI, 16), IISWAP, 16);
    CCCNOT = Mul(Mul(Mul(CCCNOT, IIIT_1, 16), ICCNOT, 16), IIIH, 16);

    // print_Matrix(CCCNOT, 16);
    // cout << endl << endl;
    // complex<double> *vec;
    // bitset<2> a[4] = {1,1,0,1};
    // vec = gen_qubit(a, 4);
    // for(int i = 0; i < 16; i++){
    //     cout << i << " " << vec[i].real() << endl;
    // }
    // cout << endl;
    // vec = Mul_vector_matrix(CCCNOT, vec, 16);
    // for(int i = 0; i < 16; i++){
    //     if(vec[i].imag() == 0) cout << i << " " << vec[i].real() << endl;
    //     else cout << i << " " << vec[i] << endl;
    // }

    complex<double> **CCT = new complex<double> *[8];
    complex<double> var_3 (1/sqrt(2), 1/sqrt(2));
    for(int i = 0;i < 16; i++){
        CCT[i] = new complex<double> [8];
    }
    for(int i = 0; i < 8; i++){
        for(int j = 0; j < 8; j++){
            if(i==j) CCT[i][j] = 1;
            else CCT[i][j] = 0;
        }
    }
    CCT[7][7] = var_3;
    // print_Matrix(CCT, 8); cout << endl;

    complex<double> **IIIIH;
    complex<double> **ICCCNOT;
    complex<double> **IIISWAP;
    complex<double> **CCCNOTI;
    complex<double> **IIIIT_1;
    complex<double> **IIIIT;
    complex<double> **IIITT;
    complex<double> **CCCNOTH;
    complex<double> **CCTT_1I;
    complex<double> **CCCCNOT;

    IIIIH = tensor_product(exp_product(I, 2, 4), H, 16, 2);
    ICCCNOT = tensor_product(I, CCCNOT, 2, 16);
    IIIIT_1 = tensor_product(exp_product(I, 2, 4), T_1, 16, 2);
    IIISWAP = tensor_product(exp_product(I, 2, 3), SWAP, 8, 4);
    CCCNOTI = tensor_product(CCCNOT, I, 16, 2);
    IIIIT = tensor_product(exp_product(I, 2, 4), T, 16, 2);
    IIITT = tensor_product(tensor_product(exp_product(I, 2, 3), T, 8, 2), T, 16, 2);
    CCCNOTH = tensor_product(CCCNOTH, H, 16, 2);
    CCTT_1I = tensor_product(tensor_product(CCT, T_1, 8, 2), I, 16, 2);

    CCCCNOT = Mul(Mul(Mul(CCCNOTI, CCTT_1I, 32), CCCNOTH, 32), IIITT, 32);
    CCCCNOT = Mul(Mul(Mul(CCCCNOT, IIISWAP, 32), CCCNOTI, 32), IIISWAP, 32);
    CCCCNOT = Mul(Mul(Mul(CCCCNOT,IIIIT_1, 32), ICCCNOT, 32), IIIIT, 32);
    CCCCNOT = Mul(Mul(Mul(CCCCNOT, IIISWAP, 32), CCCNOTI, 32), IIISWAP, 32);
    CCCCNOT = Mul(Mul(Mul(CCCCNOT,IIIIT_1, 32), ICCCNOT, 32), IIIIH, 32);
    
    // print_Matrix(CCCCNOT, 32);
    cout << endl << endl;
    complex<double> *vec;
    bitset<2> a[5] = {1,1,1,1,0};
    vec = gen_qubit(a, 5);
    for(int i = 0; i < 32; i++){
        cout << i << " " << vec[i].real() << endl;
    }
    cout << endl;
    vec = Mul_vector_matrix(CCCCNOT, vec, 32);
    for(int i = 0; i < 32; i++){
        if(vec[i].imag() == 0) cout << i << " " << vec[i].real() << endl;
        else cout << i << " " << vec[i] << endl;
    }
}