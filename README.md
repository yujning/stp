# stp

#### Introduction
Semi-Tensor Product (STP) engine for Electronic Design Automation (EDA)

#### Using stp as a stand-alone tool

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

The **test** and **example** directories are compiled in defalut, if you want turn it
off, please use
```bash
cmake -DSTP_EXAMPLES=OFF -DSTP_TEST=OFF ..
```
