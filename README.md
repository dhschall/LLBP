# The Last-Level Branch Predictor Simulator

<p align="center">
    <a href="https://github.com/dhschall/LLBP/blob/main/LICENSE">
        <img alt="GitHub" src="https://img.shields.io/badge/License-MIT-yellow.svg">
    </a>
    <a href="https://github.com/dhschall/LLBP/releases">
        <img alt="GitHub release" src="https://img.shields.io/github/release/dhschall/LLBP">
    </a>
    <!-- <a href="https://doi.org/10.5281/zenodo.5520125"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.5520125.svg" alt="DOI"></a> -->
</p>



The Last-Level Branch Predictor (LLBP) is a microarchitectural approach that improves branch prediction accuracy through additional high-capacity storage backing the baseline TAGE predictor. The key insight is that LLBP breaks branch predictor state into mulitple program contexts which can be thought of as a call chain. Each context comprises only a small number of patterns and can be prefetched ahead of time. This enables LLBP to store a large number of patterns in a high-capacity structure and prefetch only the patterns for the upcoming contexts into a small, fast structure to overcome the long access latency of the high-capacity structure (LLBP).

LLBP is presented at [MICRO 2024](https://microarch.org/micro57/).

This repository contains the source code of the branch predictor model used to evaluate LLBP's prediction accuracy. The code is based on the [CBP5 framework](http://www.jilp.org/cbp2016/), but was heavily modified and extended with various statistics to evaluate the performance of LLBP and the baseline TAGE predictor.

The aim of this framwork is to provide a fast way to evaluate different branch predictor configurations and explore the design space of LLBP. It does *not* model the full CPU pipeline but only the branch predictor.
The framework supports a timing aproximation by clocking the predictor for every taken branch and/or more than 8 executed instructions between branches. While we found that this approximation is reasonable accurage to get understand the impact of late prefetches, it is only a rough estimation. For the exact timing the predictor needs to be integrated with a full CPU simulator like ChampSim or gem5.
> We are currently working on integrating LLBP with gem5 and will release the code once it is ready.



## Prerequisites

The infrastructure has been tested with the following system configuration:

* Ubuntu 22.04.2 LTS
* gcc 11.4.0
* cmake 3.22.1


## Install Dependencies

```bash
# Install cmake
sudo apt install -y \
        cmake \
        libboost-all-dev \
        build-essential \
        pip

# Python dependencies for plotting.
pip install pandas matplotlib ipykernel

```


## Build the project

```bash
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Debug ..
cd ..

cmake --build ./build -j 8

```

## Getting the traces

The traces use to evaluate where generated using gem5 and converted into the ChampSim format. Four additional traces provided by Google where converted into the ChampSim format. The traces are available with the artifact at []().
The script `download_traces.sh` in the utils folder will download all traces from Zenodo and stores them into the `traces` directory.:

```bash
./utils/download_traces.sh
```


## Run the simulator

The simulator can be run with the following command and takes as inputs the trace file, the branch predictor model, the number of warmup instructions, and the number of simulation instructions.
The branch predictor model can be either `tage64kscl`, `tage512kscl`, `llbp` or `llbp-0lat`.

```bash
./build/predictor <trace> --model <predictor> -w <warmup instructions> -n <simulation instructions>
```

For convenience, the simulator contains a script to run the experiments on all evaluated benchmarks for a given branch predictor model (`./eval_benchmarks.sh <predictor>`).
The results in form of a stats file are stored in the `results` directory. Note, the simulator will print out some intermediate results after every 5M instructions which is useful to monitor the progress of the simulation.


## Reproducing main results (Figure 9)

To reproduce the main results of the paper - the reduction of mispredictions (Figure 9) we provide a separate script (`./eval_all.sh`) that runs the experiments for all evaluated branch predictor models and benchmarks. The script can be run as follows:

```bash
./eval_all.sh
```
The Jupyter notebook (`./analysis/mpki.ipynb`) can be used to parse the statistics file and generate the graph. Open the file and hit `Run All` to generate the graph.



## Code Walkthrough

Misc:
* The `main.cc` file contains the main entry point of the simulator. It reads the trace file, initializes the branch predictor model, and runs the simulation.
* The `bpmodel` directory contains the implementation of the TAGE-SC-L and LLBP branch predictor models.

TAGE:
* TAGE-SC-L is split into TAGE and SC-L components. The code is taken from the CBP5 framework and modified to include additional statistics to evaluate the branch predictor.

LLBP:
* LLBP derives from the TAGE-SC-L base class and overrides several of the methods to implement the LLBP functionality.
* The large structure is called TODO..


## Citation
If you use our work, please cite paper:
```
@inproceedings{schall2024llbp,
  title={The Last-Level Branch Predictor},
  author={Schall, David and Andreas Sandberg and Boris Grot},
  booktitle={Proceedings of the 57th Annual IEEE/ACM International Symposium on Microarchitecture},
  year={2024}
}
```

## License

Distributed under the MIT License. See `LICENSE` for more information.

## Contact

David Schall - david.schall@tum.de

## Acknowledgements
We thank all the anonymous reviewers of MICRO and the artifact evaluation team for their valuable feedback. Furthermore the members of the EASE-lab team at the University of Edinburgh as well as Arm Ltd. for their support and feedback.