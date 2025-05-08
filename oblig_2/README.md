mkdir build && cd build
cmake ..
make

# Run 1D simulation
./spin_wolff --dim 1 --L 16 --T 0.5

# Run 2D simulation
./spin_wolff --dim 2 --L 16 --T 1.0
