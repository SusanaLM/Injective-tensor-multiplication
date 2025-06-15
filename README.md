# Injective-tensor-multiplication
<p align="center"><a href="put_link_here"><img src='https://img.shields.io/badge/arXiv-Paper-red?logo=arxiv&logoColor=white' alt='arXiv'></a>

This reposity contains the code of the paper "Computation of Tensor Functions under Tensor-Tensor Multiplications with Linear Maps"
## Contents


* **Overview**
  - Injective SVD under a Johnson-Lindenstrauss-type embedding provides similar performance to the surjective SVD with a data-dependent matrix.

    
 
* **Example Usage**: 

  ![standart](img/table.png)



* **Getting Started**
  - Clone this repo:
```bash 
    git clone https://github.com/SusanaLM/Injective-tensor-multiplication.git

    cd Poset_filters
```    

  - prerequisites

    MATLAB >= R2018b
       
    Image Processing Toolbox

  - location of:
    - code: [Hyperspectral_injective MATLAB](ex_hyperspectral_injective.py)
    - issue tracker : [report issues](https://github.com/SusanaLM/Injective-tensor-multiplication/issues)



* **Notes**
  - version : v1.0
  - If the input has odd dimentions, the code automatically adds padding one on the right and/or bottom.
  - The combination ReLU followed by a poset filter seems to work well.



* **Colophon**
  - Credits -- code, algorithm, implementation/deployment, testing and overall direction: Susana Lopez Moreno.
  - Copyright and License -- see [LICENSE](https://github.com/mendozacortesgroup/Poset-filters/tree/main?tab=MIT-1-ov-file#readme) file.
  - How to contribute: submit issues.
  - This work was supported by the Korea National Research Foundation (NRF) grant funded by the Korean government (MSIT) (RS-2024-00406152). J.-H. Ju was supported by the Basic Science Program of the NRF of Korea (NRF-2022R1C1C101\\ 0052) and by the Basic Research Laboratory (grant MSIT no. RS-202400414849).
  - References:  put_arxiv_link
  
* **Citation**
If you use this code for your research, please cite our paper:

```
@misc{poset_filters,
  title={Computation of Tensor Functions under Tensor-Tensor Multiplications with Linear Maps}, 
      author={Jeong-Hoon Ju and Susana Lopez-Moreno},
      year={2025},
      eprint={put_eprint_here},
      archivePrefix={arXiv},
      primaryClass={},
      url={put_link_here}, 
}
```
