# Finding Highly Similar Regions of Genomic Sequences through Homomorphic Encryption

## Bitwise implementation HSR<sub>B</sub> using TFHE scheme

Stored in `bitwise` directory

### Requirement

[TFHE](https://tfhe.github.io/)<br/>
OpenMP

### Build

```
make
```

### Run

```
./main <file_X> <file_Y> <|X|> <|R|> <w_int>
```

### Sample run

```
./main ../dataset/paper_x.dat ../dataset/paper_y.dat 5 29 9
```

## Wordwise implementation HSR<sub>W</sub> using HEAAN scheme

Stored in `wordwise` directory

### Requirement

[HEaan.STAT](https://www.cryptolab.co.kr/eng/product/heaan.php) (license required)

### Build

```
cmake ./
make

```

### Run

```
./main <file_X> <file_Y> <|X|> <|R|> <l> <q> <d>
```

### Sample run

```
./main ../dataset/paper_x.dat ../dataset/paper_y.dat 5 29 50 29 11
```
