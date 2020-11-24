############################################################################################
## M-AMBI revisited: looking inside a widely-used benthic index.
## Hydrobiologia
## Sigovini, M., E. Keppel & D. Tagliapietra.
## Corresponding author: Marco Sigovini, CNR - ISMAR, e-mail: marco.sigovini@ve.ismar.cnr.it
##
## ESM1: R script for calculation of M-AMBI, M-AMBI* and S-AMBI algorithms.
##
##
## This script is provided under the terms of the GNU General Public License: 
## http://www.gnu.org/licenses/gpl.txt
#############################################################################################
##
## __________________________________________________________________________________________________________________
##
## mambisimpl Multimetric index M-AMBI and its simplified versions M-AMBI* and S-AMBI
## __________________________________________________________________________________________________________________
##
## Description
##
## The function calculates the index M-AMBI (Muxika et al. 2007). It also calculates 
## the simplified versions M-AMBI* and bivariate S-AMBI, either with metrics normalisation 
## or standardisation (Sigovini et al. 2013).
## The three metrics (two in the case of S-AMBI) can be either provided by the user or 
## calculated by the function. In the latter case, diversity measures that can be calculated
## are species richness (S), Margalef index (d) and Shannon index (H'), whereas AMBI-BC is 
## calculated as the biotic sensitivity index.
##
##
## Usage
##
##  Direct calculation: 
##  mambisimpl(X, eg, metrics = c("S", "H1", "AMBI"), trasf = "f", high = "def", bad = "def", st = rownames(X))
## 
##  Calculated metrics supplied by user:
##  mambisimpl(m.values, metrics = c("S", "H1", "AMBI"), trasf = "f", high = "def", bad = "def")
##
##
## Arguments
##
## X: Community data matrix or data frame (the format follows the vegan package 
##    (Oksanen et al. 2012)).
## eg: Ecological Groups (EGs) for a set of taxa (one-column matrix or data frame, with taxa
##    names as row names). EGs values range from 1 to 5 (if not assigned: NA). The set of 
##    taxa can be more or less inclusive of the set of taxa of the community data matrix.
## metrics: Vector which lists the names of two or three selected metrics to be calculated 
##    on 'x' and 'eg' arguments or, alternatively, directly provided by the user in the 
##    'm_values' argument.
## 
## M-AMBI* integrates two diversity metrics and a biotic sensitivity index, whereas S-AMBI
## includes, together with the biotic index, just one diversity metric. Diversity measures 
## that can be calculated by the function are species richness ("S"), Margalef index ("d")
## and Shannon index ("H1"), by calling functions from the vegan package (Oksanen et al. 2012). 
## The AMBI-BC index ("AMBI"; Borja et al. 2000) is available as the biotic sensitivity index.
## Default metrics are "S", "H1" and "AMBI".

## m.values: Metrics directly provided by the user as two- or three-column matrix or data 
##   frame (with samples names as row names), in place of 'x' and 'eg'. The order must match 
##   the 'metrics' argument.
## trasf: Metrics transformation: can be either min-max normalisation ("n"), standardisation 
##   ("s") or factor analysis (three metrics need to be supplied) as in the original M-AMBI 
##   algorithm ("f", default value).
## high: "High" reference values: can be set as "def" (highest values in the dataset, or the 
##   lowest in the case of AMBI; default value) or provided by the user (the order must match 
##   the 'metrics' argument).
## bad: "Bad" reference values: can be set as "def" (theoretical minima: S = 0, d = 0, H' = 0,
##   AMBI = 6; default value) or provided by the user (the order must match the 'metrics' 
##   argument).
## st: Vector with the same length of the number of rows in the community data matrix, 
##   providing the sample codes to pool the replicates. In case replicates are provided, 
##   AMBI is calculated on each of them (Borja & Muxika 2005), whereas other metrics are
##   calculated on the pooled samples. If the argument is not provided, the function assumes
##   that each row of the matrix is a sample.
##
##
## Details
##
## See Sigovini et al. 2013.
##
##
## Value
##
## A list with two elements. The first element is a data frame with the results of AMBI 
## for each sample: total abundance N, percentage distribution of abundances among EG 1 
## to 5 (the individuals not assigned to any EG are not counted in the total), number of 
## individuals not assigned to any EG (NA), and the benthic coefficient AMBI. The second 
## element is a data frame with the original metrics and the selected multimetric index 
## for each sample. The fictitious samples used as reference values are also included in
## the last two rows. If the metrics' values are provided by the user, the first element 
## is empty.
##
##
## Notes
##
## The psych package (Revelle 2012) is required.
## When metrics are calculated internally, the vegan package (Oksanen et al. 2012) is 
## also required.
## The internal function AMBI follows Borja et al. (2000) and can be used independently. 
## Other guidelines to the use of AMBI can be found in Borja & Muxika (2005).
## For taxa with 'unknown' EG, it is possible either to ignore their presence, if it is 
## considered not relevant in the calculation, by simply removing them from the data set; 
## or to consider the EG as "not assigned", by using the values 0 or NA as EG.
##
##
## Author
##
## Marco Sigovini
##
##
## References
##
## Borja, Á., Franco, J., Pérez, V., 2000. A marine biotic index to establish the 
##    ecological quality of soft-bottom benthos within European estuarine and Coastal 
##    environments. Marine Pollution Bulletin 40, 1100-1114.
## Borja, Á., Muxika, I., 2005. Guidelines for the use of AMBI (AZTI's Marine Biotic 
##    Index) in the assessment of the benthic ecological quality. Marine Pollution Bulletin
##    50, 787–789.
## Muxika, I., Borja, Á., Bald, J., 2007. Using historical data, expert judgement and 
##    multivariate analysis in assessing reference conditions and benthic ecological status, 
##    according to the European Water Framework Directive. Marine Pollution Bulletin 55, 16-29.
## Oksanen, J., Blanchet, F.G., Kindt, R., Legendre, P., Minchin, P.R., O'Hara, R.B., Simpson, 
##    G.L., Solymos, P., Stevens, M.H.H., Wagner, H., 2012. vegan: Community Ecology Package. 
##    R package version 2.0-3. http://CRAN.R-project.org/package=vegan
## Revelle, W., 2012. psych: Procedures for Personality and Psychological Research Northwestern
##    University, Evanston. R package version 1.2.1. http://personality-project.org/r/psych.manual.pdf
## Sigovini, M., Keppel, E., Tagliapietra, D., 2013. M-AMBI revisited: looking inside a 
##    widely-used benthic index. Hydrobiologia IN PRESS.
##
##
## Examples
##
## (See after the code.)
##
#####################################################################################################################

mambi_sigovini <- function(X, eg, m.values, metrics = c("S", "H1", "AMBI"), trasf = "f", high = "def", bad = "def", st = rownames(X)) {

	## 1) METRICS

	if (missing(m.values)) {
		x <- as.matrix(X)[order(rownames(X)), order(colnames(X))]
		
		AMBI <- function(x, eg) { 
		## function to calculate AMBI
		  
		  # assign the species from the input data table to their ecological sensitivity
		  # groups: match the species in the EG table with the species in the input data 
		  # table.
			eg <- eg[which(rownames(eg) %in% colnames(x)), , drop = F]
			
			# if there are extra species in the input data table (without an EG assignment), 
			# set their EG value to NA and add them to the end of the EG table.
			if(length(which(colnames(x) %in% rownames(eg))) < length(colnames(x))) {
				eg.na <- colnames(x)[-which(colnames(x) %in% rownames(eg))]
				eg <- rbind(eg, rep(NA, length(eg.na)))
				rownames(eg)[(nrow(eg)-length(eg.na)+1):nrow(eg)] <- colnames(x)[-which(colnames(x) %in% rownames(eg))]
			}
			
			# order the EG table alphabetically, to put any newly added species in the correct
			# spot.
			eg <- eg[order(rownames(eg)), , drop = F]
			
			# make a new table to hold the EG assignments, which will be used in the AMBI 
			# calculation - 5 columns for the 5 EGs and an additional column for species without
			# assignment (NAs), and fill it with 0 only for now.
			eg.matr <- as.data.frame(matrix(0, nrow = ncol(x), ncol = 6))
			colnames(eg.matr) <- c("EG1", "EG2", "EG3", "EG4", "EG5", "NA")
			rownames(eg.matr) <- rownames(eg)
			
			# match the original EG table and the calculation EG table for each species, 
			# resulting in a matrix of 1 and 0: each species is assigned a value of 1 in the
			# column corresponding to its EG, and 0 in all other columns.
			options(warn = -1)
			eg.matr[which(eg[, 1] == 1), 1] <- 1
			eg.matr[which(eg[, 1] == 2), 2] <- 1
			eg.matr[which(eg[, 1] == 3), 3] <- 1
			eg.matr[which(eg[, 1] == 4), 4] <- 1
			eg.matr[which(eg[, 1] == 5), 5] <- 1
			eg.matr[which(eg[, 1] == 0), 6] <- 1
			eg.matr[which(is.na(eg[, 1]) == T), 6] <- 1
			options(warn = 0)

			# multiply the species abundance input data matrix by the EG matrix to obtain a 
			# table of abundances in each EG (and NA) (columns), for each station/replicate 
			# (rows).
			ambi <- as.data.frame(x %*% as.matrix(eg.matr))
			
			# add in the beginning of this table a column N for the total abundance 
			# (N = sum of each row), and another in the end for AMBI - for now empty (filled 
			# with NAs). 
			ambi <- cbind("N" = apply(x, 1, sum), ambi, "AMBI" = rep(NA, nrow(x)))
			
			# for each station/replicate (= row), calculate the proportion of each EG 
			# (excluding the NAs)
			for (i in 1:nrow(x)) {
				for (j in 2:6) {
					ambi[i, j] <- ambi[i, j] / (ambi[i, "N"] - ambi[i, "NA"]) * 100
				}
			}
			
			# calculate AMBI for each station/replicate, rounding its value to 3 digits
			for (i in 1:nrow(x)) {
				ambi[i, "AMBI"] <- round((1.5*ambi[i, "EG2"] + 3*ambi[i, "EG3"] + 4.5*ambi[i, "EG4"] + 6*ambi[i, "EG5"]) / 100, 3)
			}

			# rename the EG columns to EG(%); round the proportions to 3 digits
			colnames(ambi)[2:6] <- c("EG1(%)", "EG2(%)", "EG3(%)", "EG4(%)", "EG5(%)")
			ambi[, 2:6] <- round(ambi[, 2:6], 3)
			
			
			return(ambi)
			
		}


		# calculate AMBI on the input data using the custom function just defined
		ambi <- AMBI(x, eg)

		## samples with N = 0 (azoic): AMBI = 7
		if(any(is.nan(ambi$AMBI))) {
		  ambi$AMBI[which(is.nan(ambi$AMBI))] <- 7
		}

		## calculation of AMBI-BC on each replicate (other metrics are calculated on 
		## pooled replicates)
		
		# get the station names; sort them in the order of the input data frame, and 
		# convert to factor (to use levels as replicates in the calculations of the diversity
		# metrics)
		st <- st[order(rownames(X))]
		st <- as.factor(st)

		# calculate average AMBI by station and exctract only this value to a 
		# separate vector
		ambi.mu <- round(aggregate(ambi$AMBI ~ st, FUN = mean)[2], 3)[, 1]
		
		# pool replicates in order to calculate the diversity metrics: sum species 
		# abundances by station/replicate and store in new table - to be used for the 
		# calculation of the other metrics
		x.st <- aggregate(x ~ st, FUN = sum)
		
		# add row names (from 1st column); then remove the first column (not needed any 
		# more)
		rownames(x.st) <- x.st[, 1]
		x.st <- x.st[-1]
		
		# add column names (= species names) - same as in the input data table
		colnames(x.st) <- colnames(x)

		## calculation of richness/diversity metrics - on pooled replicates!
		require(vegan)
		
		m.values <- data.frame(
			"S" = specnumber(x.st), ## species richness
			"d" = round((specnumber(x.st) - 1)/(log(apply(x.st, 1, sum))), 3), ## Margalef index
			"H1" = round(diversity(x.st, index = "shannon", base = 2), 3), ## Shannon index
			"AMBI" = ambi.mu) ## AMBI-BC
		
		sample.names <- rownames(x.st)
		
		# extract the chosen metrics' values and put them in a new output table 
		# (default: S, H1, BC/AMBI; if other diversity measure (Margalef's d) - should 
		# be specified at function call)
		METRICS <- NULL
		for (i in 1:length(metrics)) {
			METRICS <- cbind(METRICS, m.values[, metrics[i]])
		}
    
		colnames(METRICS) <- metrics
		rownames(METRICS) <- sample.names

	}

  # if argument m.values (= values of already computed metrics) is provided in the 
  # function call, use them instead of calculating from scratch:
  else {
		ambi <- NULL
		sample.names <- rownames(m.values)
		
		METRICS <- as.matrix(m.values)
		colnames(METRICS) <- metrics
		rownames(METRICS) <- sample.names
	}

	## 2) "HIGH" AND "BAD" REFERENCES VALUES ADDED AS FICTITIOUS SAMPLES

	options(warn = -1)
	# if high and bad reference values were provided at function call, use those values
	if (is.numeric(high) == T) METRICS.high <- high

	if (is.numeric(bad) == T) METRICS.bad <- bad
  
	# otherwise, use as HIGH reference values the minimum value for BC/AMBI and the 
	# maximum values for all other metrics in the current dataset
	if (high == "def") {
		METRICS.high <- rep(NA, length(metrics))
		
		for (i in 1:length(metrics)) {
			if (metrics[i] == "AMBI") METRICS.high[i] <- min(METRICS[, i])
			else METRICS.high[i] <- max(METRICS[, i])
		}
	}
  
	# ... and use as BAD reference values BC/AMBI = 6 and all other metrics = 0 
	if (bad == "def") {
		METRICS.bad <- rep(NA, length(metrics))
		
		for (i in 1:length(metrics)) {
			if (metrics[i] == "AMBI") METRICS.bad[i] <- 6
			else METRICS.bad[i] <- 0
		}
	}

	options(warn = 0)
	
	# add the high and bad reference values as rows at the end of the final output 
	# data table; add the row names (stations + 2 extra rows - bad and high)
	METRICS.tot <- rbind(METRICS, METRICS.bad, METRICS.high)
	rownames(METRICS.tot)[1:nrow(METRICS)] <- sample.names
	rownames(METRICS.tot)[nrow(METRICS) + 1] <- "B"
	rownames(METRICS.tot)[nrow(METRICS) + 2] <- "H"

	## 3) STANDARDISATION OR NORMALISATION - M-AMBI
	
	# standardisation of each column (= metric) by subtraction of the mean (omitting NAs) 
	# (= centering), and dividing the centered columns (metrics) by their standard 
	# deviations.
	if (trasf == "s") METRICS.tr <- scale(METRICS.tot) ## -> M-AMBI*, S-AMBI
	
	
	NORM <- function(data) {
	  ## function performing min-max normalisation of each metric: the "bad" reference 
	  ## value is subtracted from each (observed) value, and the result is divided by the
	  ## whole range for that metric ("high" - "bad").
	  
		norm <- NULL

		for (i in 1:ncol(data)) {
			norm <- cbind(norm, (data[, i] - data["B", i]) / (data["H", i] - data["B", i]))
		}
		
		rownames(norm) <- rownames(data)
		colnames(norm) <- colnames(data)
		
		return(norm)
	}

	
	# normalisation of each column (= metric) using the normalisation function just defined
	if (trasf == "n") METRICS.tr <- NORM(METRICS.tot) ## -> M-AMBI*(n), S-AMBI(n)


	## 4) FACTOR ANALYSIS - original flavour M-AMBI
	
	if (trasf == "f") {
	  # ## calculation with package 'psych'
	  # require(psych)
	  # 
	  # # perform PCA on the standardised metrics, extracting the 3 best components (factors),
	  # # with orthogonal rotation of the solution attempting to maximize the variation
	  # # explained, and reporting the factor scores.
	  # METRICS.fa <- principal(scale(METRICS.tot), nfactors = 3, rotate = "varimax", scores = T)
	  # 
	  # # estimate factor scores using one of the methods in the psych package. In all cases,
	  # # factor score estimates are based upon the data matrix times a weighting matrix which
	  # # weights the observed variables; in the method implemented here, weights are just
	  # # component loadings.
	  # # METRICS.scores <- factor.scores(METRICS.tot, f = METRICS.fa, method = c("components"))$scores
	  # METRICS.scores <- scale(METRICS.tot) %*% METRICS.fa$loadings
	  # 
	  # colnames(METRICS.scores) <- c("x", "y", "z")

		
		## direct calculation, which produces the factor scores with the same signs of the scores
	  ## produced by the AZTI-Tecnalia AMBI software.

		options(warn = -1)
		METRICS.fa <- princomp(METRICS.tot, cor = T, covmat = cov(METRICS.tot))
		options(warn = 0)
		
		# weight the metrics' loadings:
		
		# first, construct a diagonal matrix from the metrics' standard deviations (sd as
		# the diagonal & 0 as off-diagonal values), then multiply by the loadings matrix
		METRICS.fa.load <- loadings(METRICS.fa) %*% diag(METRICS.fa$sdev)  # this line seems sort of useless, because the object gets overwritten in the next one...
		
		# get the eigenvectors and rescale them by the square root of the eigenvalues to obtain
		# the loadings
		METRICS.fa.load <- eigen(cor(METRICS.tot))$vectors %*% diag(sqrt(eigen(cor(METRICS.tot))$values))
		
		METRICS.fa.load.varimax <- loadings(varimax(METRICS.fa.load))
		METRICS.scores <- scale(METRICS.tot) %*% METRICS.fa.load.varimax
		colnames(METRICS.scores) <- c("x", "y", "z")
		
		METRICS.tr <- METRICS.scores ## -> M-AMBI
	}

	## 5) EQR CALCULATION
	
	EQR <- function(data) {
	  ## function calculating the EQR (= M-AMBI), using the high and bad refence values

	  # calculate the metrics' ranges = high reference value - bad reference value 
	  # (found in the 2 last rows of the metrics' matrix)
	  segm <- data[nrow(data),] - data[(nrow(data)-1),]

	  # make an empty matrix of the same size as the one containing the calculated 
	  # metrics (including the reference bad and high values), and fill it with each 
	  # respective "bad" reference value (found in the last but one row of the metrics' 
	  # matrix)
	  vett <- matrix(NA, nrow = nrow(data), ncol = ncol(data))
		
		for (k in 1:ncol(data)) {vett[, k] <- data[(nrow(data)-1), k]}
		
	  # subtract the "bad" reference value for each metric from each observed value 
	  # (=station/replicate)
		vett <- data - vett
		
		# normalize values to the range 0-1; calculate EQR = Euclidean distance between 0
		# and the score's orthogonal projection on the line defined by the reference values
		ris <- round((vett %*% segm / sqrt(sum(segm*segm))) / sqrt(sum(segm*segm)), 3)
		
		return(ris)
	}

	eqr <- EQR(METRICS.tr)

	
	## Output tables for all available alternative algorithms: 
	
	## get the index name
	
	# direct calculation using min-max normalisation = M-AMBI*(n) - the arithmetic mean of 
	# the 3 metrics.
	if (trasf == "n") eqr <- round(apply(NORM(METRICS.tot), 1, mean), 3) 
	
	# without the factor analysis = M-AMBI* - the standardised metrics are projected directly 
	# onto the axis identified by the high and bad reference values
	if (ncol(METRICS.tot) == 3) name <- "M-AMBI*"
	
	# only 2 metrics used (AMBI as sensitivity; either S or H1 as diversity) = S-AMBI - no 
	# factor analysis but min-max normalisation with bad and high reference values as extrema.
	if (ncol(METRICS.tot) == 2) name <- "S-AMBI"
	
	# M-AMBI*(n)
	if (trasf == "n") name <- paste(name, "(n)", sep = "")
	
	# original flavor M-AMBI with factor analysis and all the bells and whistles
	if (trasf == "f") name <- "M-AMBI"

	
	## construct the output index tables
	# S-AMBI
	if(ncol(METRICS.tot) == 2) {
		indices <- data.frame(cbind(METRICS.tot, rep(NA, length(eqr)), eqr))
		colnames(indices) <- c(metrics, "NA", name)
	}
  
	# M-AMBI* & M-AMBI*n
	if(ncol(METRICS.tot) == 3) {
		indices <- data.frame(cbind(METRICS.tot, eqr))
		colnames(indices) <- c(metrics, name)
	}

	# M-AMBI classical
	if(trasf == "f") {
		indices <- data.frame(cbind(METRICS.tot, eqr, round(METRICS.tr, 6)))
		colnames(indices) <- c(metrics, name, colnames(METRICS.tr))
	}

	
	return(list(ambi, indices))
	
}


#####################################################################################################################
## EXAMPLE

## data set (provided with the original AZTI-Tecnalia AMBI software)
# 
# taxa.ex <- c("Abarenicola_claparedei", "Abra_alba", "Abra_prismatica", "Abra_sp", "Alkmaria_romijni", "Astarte_sp",
# "Bathyporeia_elegans", "Callianassa_tyrrhena", "Capitella_capitata", "Carcinus_maenas", "Caulleriella_alata",
# "Caulleriella_zetlandica", "Corbula_gibba", "Corophium_sp", "Corophium_volutator", "Cossura_sp",
# "Cumopsis_fagei", "Cyathura_carinata", "Desdemona_ornata", "Diogenes_pugilator", "Diptera",
# "Dispio_uncinata", "Eteone_spetsbergensis", "Eumida_sp", "Glycera_tridactyla", "Glycera_unicornis",
# "Hediste_diversicolor", "Heteromastus_filiformis", "Hydrobia_ulvae", "Lekanesphaera_hookeri",
# "Lekanesphaera_rugicauda", "Lepidochitona_cinerea", "Lunatia_alderi", "Malacoceros_fuliginosus",
# "Malacoceros_sp", "Mediomastus_fragilis", "Melita_palmata", "Mysidacea", "Nassarius_pygmaeus",
# "Nassarius_reticulatus", "Nematoda", "Nemertea", "Oligochaeta", "Ovatella_myosotis", "Paradoneis_armata",
# "Paragnathia_formica", "Parvicardium_scabrum", "Pectinaria_koreni", "Polydora_ligni", "Prionospio_fallax",
# "Prionospio_steenstrupi", "Pseudopolydora_paucibranchiata", "Pseudopythina_macandrewi", "Pygospio_elegans",
# "Scolaricia_sp", "Scrobicularia_plana", "Solen_marginatus", "Spio_martinensis", "Spiochaetopterus_costarum",
# "Streblospio_shrubsolii", "Tapes_decussatus", "Tellina_compressa", "Tellina_sp", "Tubulanus_polymorphus",
# "Urothoe_brevicornis")
# 
# abundances.ex <- c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 1, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 29, 0, 0, 0, 35, 23, 180, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
# 0, 0, 0, 0, 0, 1, 10, 505, 518, 117, 80, 58, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 2, 0, 0, 0, 0, 1,
# 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
# 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 6, 1, 0,
# 0, 0, 0, 2, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 9, 14, 5, 21, 0, 32, 11, 7, 6, 5, 2, 4, 17, 1, 38, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 5, 7, 6, 0, 1, 2, 5, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 2, 0, 1, 0, 0, 0, 0, 0, 0, 6, 13, 8,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 5, 0, 1, 0, 0, 3, 2, 1, 2, 33, 0, 20, 54, 116, 64, 1, 1, 0, 0, 4, 5,
# 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 12, 4, 19, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 16, 5, 19, 0, 0, 2, 1, 0, 1, 0, 0, 9, 0, 0, 0, 0, 0, 2, 1, 2, 0, 0, 0, 0, 0, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 3, 5, 4, 2, 5, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 6, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 4, 0, 0, 0, 2, 7, 28, 0, 1, 38, 266, 62, 44, 15, 3,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 4, 0, 27, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 6, 2, 2,
# 0, 0, 0, 0, 0, 0, 3, 3, 4, 0, 0, 0, 0, 0, 0, 0, 0, 1, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 4, 0, 0, 0, 0, 0,
# 0, 0, 0, 2, 7, 0, 0, 2, 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 1, 1, 3, 0, 3,
# 2, 1, 8, 13, 0, 3, 24, 3, 29, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 3, 0, 0, 0, 5, 5, 28, 1, 4, 2, 0, 12, 5, 0, 0, 0, 0, 0, 10, 4, 11, 0, 1, 0,
# 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
# 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
# 
# st.ex <- paste("st_", as.character(sprintf("%02d", c(rep(1, 3), rep(2, 3), rep(3, 3), rep(4, 3), rep(5, 3), rep(6, 3), rep(7, 3), rep(8, 2), rep(9, 2), rep(10, 2), rep(11, 2)))), sep = "")
# rep.ex <- paste(st.ex, c(rep(c("a", "b", "c"), 7), rep(c("a", "b"), 4)), sep = "")
# dataset.ex <- matrix(abundances.ex, ncol = 65)
# rownames(dataset.ex) <- rep.ex
# colnames(dataset.ex) <- taxa.ex
# 
# dataset.ex <- dataset.ex[, -which(colnames(dataset.ex) == "Mysidacea")] ## taxon ignored (as in the original example)
# dataset_st.ex <- aggregate(dataset.ex ~ st.ex, FUN = sum)
# rownames(dataset_st.ex) <- dataset_st.ex[, 1]
# dataset_st.ex <- dataset_st.ex[-1]
# colnames(dataset_st.ex) <- colnames(dataset.ex)
# 
# ## EG (source: list of taxa provided with the AZTI-Tecnalia AMBI software, March 2012 version)
# eg.mar2012 <- matrix(c(1, 3, 3, 3, 3, 1, 1, 3, 5, 3, 4, 4, 4, 3, 3, 4, 2, 3, 2, 2, 4, 3, 3, 2, 2, 2, 3, 4, 3, 3, 3, 2, 2, 5, 3, 3, 1, NA, 2, 2, 3, 3, 5, NA, 3, 3, 1, 4, 4, 4, 4, 4, NA, 3, 1, 3, 1, 3, 1, 3, 1, 1, 1, 2, 1), ncol = 1)
# rownames(eg.mar2012) <- taxa.ex
# colnames(eg.mar2012) <- "EG_v.Mar2012"
# 
# ## direct calculation of the metrics
# mambisimpl(dataset_st.ex, eg.mar2012) ## M-AMBI (AMBI on stations)
# mambisimpl(dataset.ex, eg.mar2012, st = st.ex) ## M-AMBI (AMBI on replicates)
# mambisimpl(dataset.ex, eg.mar2012, trasf = "n", st = st.ex) ## M-AMBI*(n)
# mambisimpl(dataset.ex, eg.mar2012, metrics = c("S", "d", "BC"), st = st.ex) ## BAT index
# mambisimpl(dataset.ex, eg.mar2012, metrics = c("S", "BC"), trasf = "n", st = st.ex) ## S-AMBI(n) with richness
# 
# ## Directive 2000/60/EC monitoring: italian ref. values for non-tidal coastal lagoons (D.M. 260, 8/11/2010)
# IT_MAT1 <- c(25, 3.3, 1.85)
# mambisimpl(dataset.ex, eg.mar2012, high = IT_MAT1, st = st.ex)
# 
# ## metrics provided to mambisimpl by the user
# metrics.ex <- mambisimpl(dataset.ex, eg.mar2012, trasf = "n", st = st.ex)[[2]][1:11, 1:3]
# mambisimpl(m.values = metrics.ex, trasf = "n") ## M-AMBI*(n)

#####################################################################################################################
############################################################################################################MS#v1.0##
