#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <queue>
#include <algorithm>
#include <numeric>
#include <fstream> // For saving data to a file
#include <complex> // For complex numbers

using namespace std;

#define PI 3.14159265358979323846

// System parameters
const int q = 3;

vector<int> S;
vector<int> M;

enum dirs { RIGHT, LEFT, UP, DOWN };

int xpos(int i, int L) { return i % L; }
int ypos(int i, int L) { return i / L; }
int idx(int x, int y, int L) { return y * L + x; }

int Nbr(int i, int dir, int L) {
    int x = xpos(i, L);
    int y = ypos(i, L);
    switch (dir) {
        case RIGHT: return idx((x + 1) % L, y, L);
        case LEFT:  return idx((x - 1 + L) % L, y, L);
        case UP:    return idx(x, (y - 1 + L) % L, L);
        case DOWN:  return idx(x, (y + 1) % L, L);
    }
    return -1;
}

/// @brief Wolff algorithm for flipping spins in 2D. Takes different L values
/// @param start_site_idx 
/// @param L 
/// @param p_connect 
void FlipandBuildFrom(int start_site_idx, int L, double p_connect) {
    int N = L * L;
    vector<bool> visited_in_cluster(N, false);
    int oldstate = S[start_site_idx];
    int newstate;

    if (q == 2) {
        newstate = 1 - oldstate;
    } else {
        do {
            newstate = rand() % q;
        } while (newstate == oldstate);
    }

    queue<int> cluster_q;
    cluster_q.push(start_site_idx);
    visited_in_cluster[start_site_idx] = true;
    vector<int> sites_in_cluster;

    while (!cluster_q.empty()) {
        int current_site_idx = cluster_q.front();
        cluster_q.pop();
        sites_in_cluster.push_back(current_site_idx);

        for (int d = 0; d < 4; ++d) {
            int neighbor_idx = Nbr(current_site_idx, d, L);
            if (!visited_in_cluster[neighbor_idx] && S[neighbor_idx] == oldstate &&
                (rand() / (RAND_MAX + 1.0)) < p_connect) {
                visited_in_cluster[neighbor_idx] = true;
                cluster_q.push(neighbor_idx);
            }
        }
    }

    for (int site_to_flip : sites_in_cluster) {
        S[site_to_flip] = newstate;
    }
    M[oldstate] -= sites_in_cluster.size();
    M[newstate] += sites_in_cluster.size();
}

// Calculate the local "magnetization" mj
complex<double> local_magnetization(int spin_state) {
    return exp(complex<double>(0.0, (2.0 * PI / q) * spin_state));
}

// Calculate the order parameter or "magnetization per site" m
complex<double> calculate_magnetization(int L) {
    int N = L * L;
    complex<double> sum_mj(0.0, 0.0);
    for (int j = 0; j < N; ++j) {
        sum_mj += local_magnetization(S[j]);
    }
    return sum_mj / static_cast<double>(N);
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    int NESTEPS = 50000;   // Number of equilibration steps
    int NMSTEPS = 100000;  // Number of measurement MC sweeps
    int MEASURE_INTERVAL = 100; // Interval between measurements
    int NUM_RUNS = 10;      // Number of independent runs per (T, L)

    vector<int> L_values = {8, 16, 32}; // Different system sizes
    vector<double> temperatures;
    for (int i = 1; i <= 50; ++i) {
        temperatures.push_back(0.05 * i); // Temperature range
    }

    for (int L : L_values) {
        int N = L * L;
        ofstream outfile_gamma("binder_gamma_L" + to_string(L) + "_q3.dat");

        if (!outfile_gamma.is_open()) {
            cerr << "Error opening output file for Gamma (L=" << L << ")!" << endl;
            continue;
        }

        outfile_gamma << "# T/J Gamma" << endl;

        for (double T_over_J : temperatures) {
            double T = T_over_J;
            double p_connect = 1.0 - exp(-1.0 / T);

            double sum_m_abs2 = 0.0;
            double sum_m_abs4 = 0.0;
            int n_samples = 0;

            for (int run = 0; run < NUM_RUNS; ++run) {
                S.resize(N);
                M.resize(q);
                fill(S.begin(), S.end(), 0);
                fill(M.begin(), M.end(), 0);
                M[0] = N;

                for (int t = 0; t < NESTEPS; t++) {
                    FlipandBuildFrom(rand() % N, L, p_connect);
                }

                for (int sweep = 0; sweep < NMSTEPS; ++sweep) {
                    FlipandBuildFrom(rand() % N, L, p_connect);

                    if ((sweep + 1) % MEASURE_INTERVAL == 0) {
                        complex<double> m = calculate_magnetization(L);
                        double m_abs2 = std::norm(m);
                        sum_m_abs2 += m_abs2;
                        sum_m_abs4 += m_abs2 * m_abs2;
                        n_samples++;
                    }
                }
            }

            if (n_samples > 0) {
                double avg_m_abs2 = sum_m_abs2 / n_samples;
                double avg_m_abs4 = sum_m_abs4 / n_samples;
                double gamma = avg_m_abs4 / (avg_m_abs2 * avg_m_abs2);

                outfile_gamma << T_over_J << " " << gamma << endl;

                cout << "L = " << L << ", T/J = " << T_over_J << ", Gamma = " << gamma << endl;
            }
        }
        outfile_gamma.close();
        cout << "\nSimulation finished for L = " << L << ". Binder cumulant data saved to binder_gamma_L" << L << "_q3.dat" << endl;
    }

    cout << "\nSimulation finished for all L values." << endl;

    return 0;
}