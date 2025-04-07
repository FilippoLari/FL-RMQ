# FL-RMQ: A Learned Approach to Range Minimum Queries

This repository implements the **Fully-Learned Range Minimum Query** (FL-RMQ), a novel data structure for **range minimum queries** based on a surprising connection between this classical problem and the geometry of a properly defined set of points in the Cartesian plane. Building on this insight, FL-RMQ introduces a unique approach that learns and exploits the distribution of such points using error-bounded linear approximations.

Most notably, FL-RMQ gives robust theoretical guarantees and offers novel space-time trade-offs with respect to known systematic solutions.

## Building the project

```bash
git clone https://github.com/FilippoLari/FL-RMQ
cd FL-RMQ
mkdir build
cd build
cmake ..
make -j8
```

## Reproducing our results

You can find our benchmarking suite and instructions on how to reproduce our results at the following link: [RMQ-Experiments](https://github.com/FilippoLari/RMQ-experiments)

## Credits

If you use this project in a scientific article, please cite the following paper:

```tex
@inproceedings{FerraginaL25,
    author = {Paolo Ferragina and Filippo Lari},
    title = {{FL-RMQ}: A Learned Approach to Range Minimum Queries},
    booktitle = {Proc. 36th Annual Symposium on Combinatorial Pattern Matching},
    year = {2025}
}
```
