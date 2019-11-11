# ZoKrates with C verifier

See [Document](https://hackmd.io/eGyH65-1SG6oN7IlzeZaHg).
Based on ZoKrates and output verifiers with C language to deploy on OurChain.

## Getting Started

Download the whole package and compile the code.

```bash
cd ZoKrates/
cargo +nightly build --release
```

Then `cd target/release` to execute ZoKrates

## ZoKrates Example 1 -- Square Root

According to [ZoKrates document](https://zokrates.github.io/gettingstarted.html):

1. First, create the text-file `root.zok` and implement your program. In this example, we will prove knowledge of the square root `a` of a number `b`:
```bash
def main(private field a, field b) -> (field):
  field result = if a * a == b then 1 else 0 fi
  return result
```

2. Compile
```bash
./zokrates compile -i root.zok
```

3. Perform the setup phase
```bash
./zokrates setup
```

4. (For verifiers) Export a solidity verifier
```bash
./zokrates export-verifier
```

4. (For provers) Generate a proof of computation
```bash
# execute the program
./zokrates compute-witness -a 337 113569
# generate a proof of computation
./zokrates generate-proof
```

## ZoKrates Example 2 -- SudokuChecker

Another example of ZoKrates: [SudokuChecker](https://github.com/vivianjeng/zkSNARK-verifier-in-C/blob/master/ZoKrates/zokrates_cli/examples/sudokuchecker.zok)

1. create the text-file `sudokuchecker.zok`

2. Compile
```bash
./zokrates compile -i sudokuchecker.zok
```

3. Perform the setup phase
```bash
./zokrates setup
```

4. (For verifiers) Export a solidity verifier
```bash
./zokrates export-verifier
```

4. (For provers) Generate a proof of computation
```bash
# execute the program
./zokrates compute-witness -a 3 4 2 2 3 2 1 2 4 3 1 1 4 3 4 1
# generate a proof of computation
./zokrates generate-proof
```

## On-chain Verifier Usage

1. Compiling `verifier.c` on-chain

```bash
gcc verifier.c -lgmp
```

2. Prover gets proof from ZoKrates output `proof.json`

3. Prover calls the verifier program with the proof.

An example of SudokuChecker:

```bash
./a.out 12223361632793408695071575845139783748399099095148486551738160861113901529450 10773020487045276868139663904300449364583091909115903404377551147958337460347 11458271398869128052404652055082728846813650718352464338821092175948674322028 7377129170170691236169079245377746972165048032070197772644585782612415622023 9577223338187511070267836222101724545529135090815987419922695191044115707281 11564317194327007249052994800681576503120340472613585493085497052076914813118 19080400365662515268619738133110008210937427276810264496399607417116506797606 15466541046120039457992702990746787383814430749729372800776852616315204697252 3 4 2 2 3 2 1
```

4. After executing of the verifier program, it will print "transaction verified" if proof is correct; otherwise it shows nothing.