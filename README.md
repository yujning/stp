# stp

#### Introduction
Semi-Tensor Product (STP) engine for Electronic Design Automation (EDA)

[Read the full documentation.](https://stp-based-logic-synthesis-tool.readthedocs.io/en/latest/ )

#### Using stp as a stand-alone tool
Now, we use [eigen](https://eigen.tuxfamily.org/) library for matrix computation, so please install it
before running this project.

```bash
git clone https://gitlab.com/libeigen/eigen.git
cd eigen
mkdir build
cd build
cmake ..
make install
```

Then you can clone **stp** project and compile it.

```bash
git clone https://gitee.com/zfchu/stp.git   (Gitee repository)
git clone https://github.com/nbulsi/stp.git (GitHub repository) 
cd stp
mkdir build
cd build
cmake ..
make
./bin/stp
```


use cuda(DSTP_ENABLE_CUDA is set to OFF by default):
```bash
cmake .. -DSTP_ENABLE_CUDA=ON
```

sim:
```bash
sim -l yourcase.bench
```


sim with cuda:
```bash
sim -l -c yourcase.bench
```