# Finding Highly Similar Regions of Genomic Sequences through Homomorphic Encryption

## Bitwise implementation using TFHE scheme

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
./main ../dataset/paper_x.dat ../dataset/paper_y.dat 5 27 9
```

## Wordwise implementation using HEAAN scheme

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
./main ../dataset/paper_x.dat ../dataset/paper_y.dat 5 27 50 29 11
```
