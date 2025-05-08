#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <complex>
#include <ctime>
#include <queue>
#include <fstream>
#include <numeric>

using namespace std;

// Constants
#define PI 3.14159265358979323846

const int L = 16;                  // System size
const int q = 3;                   // Number of Potts states
const double T = 0.25;             // Temperature
const double p = 1 - exp(-1.0 / T); // Cluster connection probability
const int N = L;                   // Number of spins (1D system)

// Global spin and magnetization state
vector<int> S(N);                  // Spin configuration
vector<int> M(q);                  // Spin counts for each Potts state

enum Direction { RIGHT, LEFT };

int xpos(int i) { return i % L; }

// Returns neighbor index in 1D with periodic boundary conditions
int Nbr(int i, int dir) {
    int x = xpos(i);
    if (dir == RIGHT) return (x + 1) % L;
    if (dir == LEFT)  return (x - 1 + L) % L;
    return -1; // Error case
}

// Grows and flips a cluster from a seed site
void FlipandBuildFrom(int s, vector<bool>& visited) {
    int oldstate = S[s];
    int newstate = (oldstate + 1) % q;

    queue<int> cluster;
    cluster.push(s);
    visited[s] = true;

    while (!cluster.empty()) {
        int site = cluster.front(); cluster.pop();
        S[site] = newstate;
        M[oldstate]--;
        M[newstate]++;

        for (int dir = 0; dir < 2; ++dir) {
            int j = Nbr(site, dir);
            if (!visited[j] && S[j] == oldstate && (rand() / (RAND_MAX + 1.0) < p)) {
                cluster.push(j);
                visited[j] = true;
            }
        }
    }
}

// Converts spin state to complex unit vector
complex<double> local_magnetization(int spin_state) {
    return exp(complex<double>(0.0, 2.0 * PI * spin_state / q));
}

// Returns total complex magnetization
complex<double> calculate_magnetization() {
    complex<double> sum_m(0.0, 0.0);
    for (int j = 0; j < N; ++j) {
        sum_m += local_magnetization(S[j]);
    }
    return sum_m / static_cast<double>(N);
}

// Measures spin correlation function C(r)
vector<complex<double>> measureCorrelation(const complex<double>& avg_m) {
    vector<complex<double>> C(L, complex<double>(0.0, 0.0));

    for (int r = 0; r < L; ++r) {
        complex<double> sum(0.0, 0.0);
        for (int i = 0; i < L; ++i) {
            complex<double> m0 = local_magnetization(S[i]);
            complex<double> mr = local_magnetization(S[(i + r) % L]);
            sum += conj(m0) * mr;
        }
        C[r] = sum / static_cast<double>(L);
    }

    return C;
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    // Initialize spins (all in state 0)
    fill(S.begin(), S.end(), 0);
    fill(M.begin(), M.end(), 0);
    M[0] = N;

    const int NMSTEPS = 100000;          // Total measurement steps
    const int MEASURE_INTERVAL = 5;      // Wolff flips per measurement

    complex<double> avg_magnetization(0.0, 0.0);
    vector<complex<double>> C_avg(L, complex<double>(0.0, 0.0));
    int n_samples = 0;

    // Main loop
    while (n_samples < NMSTEPS / MEASURE_INTERVAL) {
        for (int m = 0; m < MEASURE_INTERVAL; ++m) {
            vector<bool> visited(N, false);
            FlipandBuildFrom(rand() % N, visited);
        }

        complex<double> current_m = calculate_magnetization();
        avg_magnetization += current_m;

        vector<complex<double>> C = measureCorrelation(current_m);
        for (int r = 0; r < L; ++r) {
            C_avg[r] += C[r];
        }

        ++n_samples;
    }

    // Final averages
    avg_magnetization /= static_cast<double>(n_samples);
    for (int r = 0; r < L; ++r) {
        C_avg[r] /= static_cast<double>(n_samples);
    }

    // Output results
    cout << "Average Magnetization (m): " << avg_magnetization << endl;
    cout << "  Real: " << avg_magnetization.real()
         << ", Imag: " << avg_magnetization.imag() << endl;

    ofstream outfile("1d.txt");
    if (outfile.is_open()) {
        cout << "\nCorrelation C(r):" << endl;
        for (int r = 0; r < L; ++r) {
            cout << r << ": " << C_avg[r] << " (Re: " << C_avg[r].real()
                 << ", Im: " << C_avg[r].imag() << ")" << endl;
            outfile << r << " " << C_avg[r].real() << " " << C_avg[r].imag() << endl;
        }
        outfile.close();
        cout << "Correlation data written to 1d.txt" << endl;
    } else {
        cerr << "Error: could not open file for writing." << endl;
    }

    return 0;
}
