# Cancer Detection using DNA Methylation Data

## Environment Setup

### Poetry

#### Requirements
- Python (>= 3.11)
- Poetry (>= 1.8.3)

#### Setup
configure poetry to create virtual environments in the project directory

```bash
poetry config virtualenvs.in-project true
```

install dependencies

```bash
poetry install
```

add dependencies

```bash
poetry add <package>
```

### renv

#### Requirements
- R (>= 4.4.1)
- Renv(>= 1.0.7)

#### Setup
install renv if not already installed

```r
install.packages("renv", repos = "https://cran.csie.ntu.edu.tw/")
```

install dependencies from the lockfile

```r
renv::restore()
```

update the lockfile

```r
renv::snapshot()
```

> [!WARNING]  
> default renv setting does not snapshot all dependencies, you may need to set `snapshot.type = "all"` as follows:
> ```r
> renv::settings$snapshot.type("all")
> ```
