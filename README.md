# Solving-TSP-3-4-B
> **Algorithms - Assignment 3 &amp; 4**
> (+ Bonus Assignment)
> 
Compile files
```bash
g++ -O2 -std=c++17 tsp_dp_hw3.cpp -o tsp_dp_hw3
```
```bash
g++ -O2 -std=c++17 tsp_aco_hw4.cpp -o tsp_aco_hw4
```

Execute `tsp_dp_hw3.exe` (dt1-dt3), `tsp_aco_hw4.exe` (dt1-dt4)
```bash
tsp_dp_hw3 --batch dataset output_hw3
```
```bash
tsp_aco_hw4 --batch dataset output_hw4 --iter 200 --pop 50 --alpha 1 --beta 5 --rho 0.5 --Q 100
```

Execute `tsp_dp_hw3.exe` (dt1), `tsp_aco_hw4` (dt1)
```bash
tsp_dp_hw3 dataset/dt1/point.txt ans_dt1.txt
```
```bash
tsp_aco_hw4 dataset/dt1/point.txt ans_dt1.txt --iter 200 --pop 50 --alpha 1 --beta 5 --rho 0.5 --Q 100
```

Execute `plot.py`
```bash
python plot.py
```
