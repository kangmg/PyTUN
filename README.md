# pyTUN

Gaussian output 파일로부터 eyring equation의 transmission coefficient를 계산해주는 코드

## Test data

* 반응 종류 : ClFHP 분자의 configurational inversion
* 계산 레벨 : CCSD/def2SVP
* 프로그램  : Gaussian 16

결과파일

|File Name|Content|
|:-:|:-:|
|reactant.out|opt. geometry freq (R-isomer)|
|product.out|opt. geometry freq (S-isomer)|
|ts.out|planar TS freq|
|reactant_sp.out|opt. geometry Single point (R-isomer)|
|product_sp.out|opt. geometry Single point (S-isomer)|
|ts_sp.out|planar TS single point|

실행 결과
```
# Default Temperature :  298.15

> python pyTUN.py

Uncorrected delta G dagger is 339.456199
Uncorrected rate is 2.1042887039E-47
The Wigner tunneling correction is 2.313175
The Skodje tunneling correction is 8.546917
The Eckart tunneling correction is 7.456219
For Wigner Tunneling:
Corrected rate is 4.8675871276E-47
Apparent delta G dagger is 337.377294
For Skodje Tunneling:
Corrected rate is 1.7985180470E-46
Apparent delta G dagger is 334.137421
For Eckart Tunneling:
Corrected rate is 1.5690038102E-46
Apparent delta G dagger is 334.475854
```

## Todo

- [x] Modifies for compatibility with Python 3
- [ ] Enables the parsing of Psi4 and ORCA output files
- [ ] Multi-Temperature mode

[Oiginal README file](https://github.com/kangmg/PyTUN/blob/master/information.md)
