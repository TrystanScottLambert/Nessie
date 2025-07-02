# Nessie

A Fast and Flexible Friends-of-Friends Group Finder Based on the GAMA Group Finder in Robotham+2011

The Nessie R package is a tool for constructing galaxy-group catalogues from redshift survey data using the Friends-of-Friends algorithm based on Robotham+2011 and described in Lambert+(in prep).

This package aims to be as user-friendly as possible and requires minimal information. The core functionality can be run on any dataset with R.A., Dec., and redshift information if the appropriate linking lengths are known.

## Installation

### Installing R

Since this is the R version of this tool, we assume that the user already has R installed. If not, please see the top of the [CRAN web page](https://cran.r-project.org/) or follow the guide [here](https://rstudio-education.github.io/hopr/starting.html). These steps are relatively simple.

### Installing Rust

Since the core functionality of Nessie is written in Rust, the Rust package manager—Cargo—is required. Fortunately, **this is very easy to install.** The [rustlang site](https://www.rust-lang.org/tools/install) should detect your operating system and tell you what command to use to install Rust using `rustup`. This should automatically install `Cargo` as well.

For Unix systems (macOS + Linux), this is as easy as running:

```bash
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh

```

in the terminal.

### Installing Nessie

The easiest method for installation is to install it using either `remotes` (needed for some HPC environments) or `devtools` (recommended) in an R session. You can install either one easily with: `install.packages("devtools")`

```r
library(devtools)
devtools::install_github("TrystanScottLambert/Nessie")

```

Nessie should now be installed as an R package.

### Helpful Additional Packages

The `R6` and `igraph` packages are required and should have been installed automatically. However, they are easy to install with `install.packages('R6')` and `install.packages('igraph')` respectively.

Although not required for the core group finding functionality, `Highlander` is needed for tuning the group finder to a known mock catalogue.


```r
devtools::install_github('asgr/Highlander')

```

Additionally, if you would like to run the package's tests, then `testthat` is also required.

## Finding Galaxy Groups in Redshift Catalogs

### Setting Up the RedshiftCatalog Object

Nessie works by first setting up a RedshiftCatalog object and then running the group finder on that catalog. This catalog requires R.A., Dec., and redshift coordinates as well as a cosmology and a running density function.

#### Setting a Cosmology

At the moment, Nessie only works for flat cosmologies and one can be created with:

```r
cosmo <- FlatCosmology$new(h = 0.7, omega_matter = 0.3)

```

#### Building a Density Function

The running density function can be created using the n(z) of the data. We include a helper function `create_density_function` which will take a distribution of redshifts and build the appropriate function.

```r
running_density <- create_density_function(redshifts, total_counts = length(redshifts), survey_fractional_area = 0.0001, cosmology = cosmo)

```

In the most basic of cases, the redshifts of the actual survey can be used. However, we recommend building an appropriate n(z) that accounts for large-scale structure fluctuations. This can be done by either using a randoms catalogue (Cole+2011) or even fitting a skewed normal function and sampling appropriately.

Note that the fractional area is required and not the total area in steradians. i.e., area_in_steradians/4π.

### Running the Group Finder

The final redshift catalog can be built like this:

```r
red_cat <- RedshiftCatalog$new(data$ra, data$dec, data$zobs, running_density, cosmo)

```

Once built, the group finder can be run:

```r
red_cat$run_fof(b0 = 0.05, r0 = 18)

```

This stores the group IDs in the RedshiftCatalog object. So a full example of reading in your data and updating it with the group catalog information may look something like this:

```r
library(Nessie)
# Preparing redshift data
data <- as.data.frame(read.csv('some_redshift_survey.csv'))
cosmo <- FlatCosmology$new(h = 0.7, omega_matter = 0.3)
running_density <- create_density_function(data$zobs, total_counts = length(data$zobs), survey_fractional_area = 0.0001, cosmology = cosmo)

# Running group catalogue
red_cat <- RedshiftCatalog$new(data$ra, data$dec, data$zobs, running_density, cosmo)
red_cat$run_fof(b0 = 0.05, r0 = 18)
data['group_ids'] <- red_cat$group_ids

```

This would result in a dataframe with R.A., Dec., redshift, and group_id **where -1 is chosen to mean that that galaxy was not found in any group.**

### Group Catalog

The group catalog is stored as a data frame that can be written to file in any way that the user wishes later. It does require the absolute magnitudes to be known beforehand and these can be calculated in any way that the user sees fit.

```r
group_catalog_df <- red_cat$calculate_group_table(data$abs_mags)

```

## Tuning Against a Mock Catalogue

The above example is easy to do if you already know what the linking lengths are, but often the choice of `b0` and `r0` is not clear. A standard practice to overcome this issue is to rely on mock catalogues of known groupings to "tune" the best values. I.e., find the values of `b0` and `r0` that best recover what is known in the mock catalogues.

Obtaining such mock catalogues is beyond the scope of this package. We assume that the user has obtained some in one manner or another or built their own.

### Tuning on Mocks

We have included a helper function `tune_group_finder` which should handle most of the tuning cases. This uses the cost function described in section 3.1 of Robotham+2011.

This function requires a list of RedshiftCatalog objects, the minimum size of a galaxy-group to be considered, the initial guess for `b0`, the initial guess for `r0`, and the bounds for both in `c(lower, upper)` format.

Importantly, each RedshiftCatalog object needs a value for the `mock_group_ids` field **where -1 means a galaxy is not in any group.** This has to be manually set.

For example, if we have a mock catalogue with "group_id" as a column header:

```r
library(Nessie)
mock_catalog_data <- as.data.frame(read.csv('mock_data.csv')) # Read in the mock data
red_cat <- RedshiftCatalog$new(mock_catalog_data$ra, mock_catalog_data$dec, mock_catalog_data$zobs, running_density, cosmo)
red_cat$mock_group_ids <- mock_catalog_data$group_id

```

Note that it is not necessary to run the group finder at this stage and is not required since `tune_group_finder` would do this automatically.

```r
b0_guess <- 0.05
r0_guess <- 18.0
b0_bounds <- c(0.04, 0.06)
r0_bounds <- c(10, 30)
min_group_size <- 5
list_of_catalogues <- c(red_cat) # One lightcone in this example

results <- tune_group_finder(list_of_catalogues, min_group_size, b0_guess, r0_guess, b0_bounds, r0_bounds)
best_fit_params <- results$parm

```

#### Tuning on Multiple Lightcones

The reason that the function requires a list of RedshiftCatalog objects, instead of expecting only a single RedshiftCatalog object, is because often multiple lightcones are built to account for cosmic variance. Each individual lightcone can become its own object and their scores will be combined together in the final calculation.

For example, if we had a redshift survey that had two fields (North and South) that both covered 0.0001 of the sky, we might make 10 mock catalogues with a "Field" column containing either 'N' for North or 'S' for South and a "CatNum" representing the catalog number (1-10). So CatNum = 5 and Field = 'N' would be the North field in mock catalog 5.

Setting these up could be done in the following way:

```r
library(Nessie)

# Read in data
mock_catalogues <- as.data.frame(read.csv('mock_catalogues.csv'))

frac_area <- 0.0001
cosmo <- FlatCosmology(0.7, 0.3)

mock_numbers <- c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
obs_fields <- c('N', 'S')
combinations <- expand.grid(mock_numbers, obs_fields)
all_fields <- paste0(combinations$Var1, combinations$Var2)

redshift_catalogues <- list()
for (ob_field in all_fields) {
  mock_obs_field <- substr(ob_field, nchar(ob_field) - 2, nchar(ob_field))
  local_catalogue <- calibration_data[calibration_data['combined_field_mock_number'] == mock_obs_field, ]
  red_cat <- RedshiftCatalog$new(local_catalogue$ra, local_catalogue$dec, local_catalogue$zobs, rho_means[[mock_obs_field]], cosmo)
  red_cat$mock_group_ids <- local_catalogue$GroupID
  redshift_catalogues[[length(redshift_catalogues) + 1]] <- red_cat
}

```

And then the tuning could be done in the exact same way:

```r
b0_guess <- 0.05
r0_guess <- 18.0
b0_bounds <- c(0.04, 0.06)
r0_bounds <- c(10, 30)
min_group_size <- 5
list_of_catalogues <- redshift_catalogues

results <- tune_group_finder(list_of_catalogues, min_group_size, b0_guess, r0_guess, b0_bounds, r0_bounds)
best_fit_params <- results$parm

```

Tuning can be slow, but since Nessie is already very highly multi-threaded, there should be no need to try and parallelize further. Furthermore, Nessie should scale with the number of cores, which means it will use the maximum number of cores available on HPC environments.
