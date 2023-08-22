# stp

#### Introduction
Semi-Tensor Product (STP) engine for Electronic Design Automation (EDA)

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
./test/run_tests
./example/matrix
```

The **test** and **examples** directories are compiled in defalut, if you want turn it
off, please use
```bash
cmake -DSTP_EXAMPLES=OFF -DSTP_TEST=OFF ..
```
