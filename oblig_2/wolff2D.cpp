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
const int LX = 16;
const int LY = 16;
const int q = 3;
const int N = LX * LY;

vector<int> S(N);
vector<int> M(q);

enum dirs { RIGHT, LEFT, UP, DOWN };

int xpos(int i) { return i % LX; }
int ypos(int i) { return i / LX; }
int idx(int x, int y) { return y * LX + x; }

/// @brief Giving the neighbour indices in the lattice  
/// @param i 
/// @param dir 
/// @return 
int Nbr(int i, int dir) {
    int x = xpos(i);
    int y = ypos(i);
    switch (dir) {
        case RIGHT: return idx((x + 1) % LX, y);
        case LEFT:  return idx((x - 1 + LX) % LX, y);
        case UP:    return idx(x, (y - 1 + LY) % LY);
        case DOWN:  return idx(x, (y + 1) % LY);
    }
    return -1;
}

/// @brief Wolff algorithm for flipping spins in 2D
/// @param start_site_idx 
/// @param p_connect 
void FlipandBuildFrom(int start_site_idx, double p_connect) {
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
            int neighbor_idx = Nbr(current_site_idx, d);
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

complex<double> local_magnetization(int spin_state) {
    return exp(complex<double>(0.0, (2.0 * PI / q) * spin_state));
}

// double calculate_magnetization() {
//     return static_cast<double>(q * M[0] - N) / (N * (q - 1));
// }
complex<double> calculate_magnetization() {
    // int N = L * L;
    complex<double> sum_mj(0.0, 0.0);
    for (int j = 0; j < N; ++j) {
        sum_mj += local_magnetization(S[j]);
    }
    return sum_mj / static_cast<double>(N);
}

int main() {
    srand(static_cast<unsigned>(time(0)));

    int NESTEPS = 50000;  
    int NMSTEPS = 100000; 
    int MEASURE_INTERVAL = 100;

    vector<double> temperatures;
    for (int i = 1; i <= 50; ++i) {
        temperatures.push_back(0.05 * i);
    }

    ofstream outfile_m("magnetization_L16_q3_2.dat");
    ofstream outfile_m2("magnetization_squared_L16_q3_2.dat");

    if (!outfile_m.is_open() || !outfile_m2.is_open()) {
        cerr << "Error opening output files!" << endl;
        return 1;
    }

    outfile_m << "# T/J <m>" << endl;
    outfile_m2 << "# T/J <|m|^2>" << endl;

    for (double T : temperatures) {
        double p_connect = 1.0 - exp(-1.0 / T);

        // Initialize spins and magnetization counts for each temperature
        fill(S.begin(), S.end(), 0);
        fill(M.begin(), M.end(), 0);
        M[0] = N;

        complex<double> avg_m = 0.0;
        complex<double> avg_m2 = 0.0;
        int n_samples = 0;

        // Equilibrate
        for (int t = 0; t < NESTEPS; t++) {
            FlipandBuildFrom(rand() % N, p_connect);
        }

        // Measurement loop
        for (int sweep = 0; sweep < NMSTEPS; ++sweep) {
            FlipandBuildFrom(rand() % N, p_connect);

            if ((sweep + 1) % MEASURE_INTERVAL == 0) {
                complex<double> m = calculate_magnetization();
                avg_m += m;
                avg_m2 += conj(m) * m;
                n_samples++;// = n_samples + 1;
            }
        }
        // complex<int> zero = 0;
        if (n_samples > 0) {
            avg_m /= n_samples;
            avg_m2 /= n_samples;

            outfile_m << T << " " << real(avg_m) << endl;
            outfile_m2 << T << " " << real(avg_m2) << endl;

            cout << "T/J = " << T << ", <m> = " << real(avg_m) << ", <|m|^2> = " << real(avg_m2) << endl;
        }
    }

    outfile_m.close();
    outfile_m2.close();

    cout << "\nSimulation finished. Data saved to magnetization_L16_q3.dat and magnetization_squared_L16_q3.dat" << endl;

    // Values at T/J = 0 and T/J -> infinity
    // cout << "\nExpected values:" << endl;
    // cout << "At T/J = 0:" << endl;
    // cout << "<m> = 1.0 (all spins aligned in one of the q states)" << endl;
    // cout << "<|m|^2> = 1.0" << endl;
    // cout << "At T/J -> infinity:" << endl;
    // cout << "<m> = 0.0 (spins are randomly distributed among q states)" << endl;
    // cout << "<|m|^2> = 1.0 / (q - 1) = " << 1.0 / (q - 1) << " for q > 2" << endl;
    // cout << "<|m|^2> = 1.0 for q = 2 (using the definition where states are 0 and 1, m = 2*N_0/N - 1, m^2 averages to 1)" << endl;



    return 0;
}