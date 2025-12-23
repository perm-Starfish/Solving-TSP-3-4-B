# Solving-TSP-3-4-B
> **Algorithms - Assignment 3 &amp; 4**
> ( + Bonus Assignment)
> 
Compile files
```bash
g++ -O2 -std=c++17 tsp_dp_hw3.cpp -o tsp_dp_hw3
```
```bash
g++ -O2 -std=c++17 tsp_aco_hw4.cpp -o tsp_aco_hw4
```
```bash
g++ -O2 -std=c++17 tsp_elastic_net.cpp -o tsp_elastic_net
```

Execute `tsp_dp_hw3.exe` (dt1-dt3), `tsp_aco_hw4.exe` (dt1-dt4), `tsp_elastic_net.exe` (dt1-dt4)
```bash
tsp_dp_hw3 --batch dataset output_hw3
```
```bash
tsp_aco_hw4 --batch dataset output_hw4 --iter 200 --pop 50 --alpha 1 --beta 5 --rho 0.5 --Q 100
```
```bash
tsp_elastic_net --batch dataset output_bonus --runs 30 --iter 0 --eval 0 --Mmul 8 --lr0 0.2 --lr1 0.02 --sig0 3 --sig1 0.5 --lam0 0.02 --lam1 0.2 --record 100
```

Execute `tsp_dp_hw3.exe` (dt1), `tsp_aco_hw4` (dt1), `tsp_elastic_net.exe` (dt1)
```bash
tsp_dp_hw3 dataset/dt1/point.txt ans_dt1.txt
```
```bash
tsp_aco_hw4 dataset/dt1/point.txt ans_dt1.txt --iter 200 --pop 50 --alpha 1 --beta 5 --rho 0.5 --Q 100
```
```bash
tsp_elastic_net dataset/dt1/point.txt ans_dt1.txt --runs 30 --Mmul 8 --record 100
```

Execute `plot`
```bash
python plot_dp.py dataset output_hw3 fig_hw3
```
```bash
python plot_aco.py dataset output_hw4 fig_hw4
```
