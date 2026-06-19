# Kovacic-Algorithm
My personal project to perform Kovacic Algorithm.

## Setup

```bash
git clone https://github.com/SeohaLim0501/Kovacic-Algorithm.git
cd Kovacic-Algorithm
```

Create and activate a Miniconda environment:

```bash
conda create -n kovacic python=3.11
conda activate kovacic
```

Install dependencies:

```bash
pip install -r requirements.txt
```

## Run

### Default Examples

```bash
python main.py
```

By default, `main.py` reads test cases from `examples.txt`.

### Integration Timeout

Integration timeout defaults to 10 seconds. Use `-t` or `--timeout` to change it:

```bash
python main.py -t 5
```

### Input Text Files

Use `-c` or `--cases` to choose a different text file:

```bash
python main.py -c input.txt
```

For your own input, you can either write to `input.txt` from the shell or create/edit any `.txt` file yourself. Then pass that file with `-c`:

```bash
python main.py -c input.txt
```

Each non-comment line in the input file can be either:

```txt
x + 77/(16*x**2)
My test | (4*x**2 + 8*x + 6)/(2*x + 1)**2
```

In PowerShell, overwrite `input.txt` with one test case:

```powershell
"x + 77/(16*x**2)" | Set-Content input.txt
```

Or include a title:

```powershell
"My test | x + 77/(16*x**2)" | Set-Content input.txt
```

Append another test case:

```powershell
"(4*x**2 + 8*x + 6)/(2*x + 1)**2" | Add-Content input.txt
```

Clear all contents of `input.txt`:

```powershell
Clear-Content input.txt
```
