
#' HistDAWass: Histogram-Valued Data Analysis
#'
#' Methods for the analysis of histogram-valued data, i.e. data described by
#' univariate histograms, in the framework of Symbolic Data Analysis. The basic
#' statistics and the methods are based on the L2 Wasserstein metric between
#' distributions, i.e. the Euclidean metric between quantile functions.
#'
#' @section Data structures:
#' \itemize{
#' \item \code{\link{distributionH}}: a histogram; \code{\link{data2hist}} builds
#' one from raw data.
#' \item \code{\link{MatH}}: a matrix of histograms (individuals by
#' histogram-valued variables).
#' \item \code{\link{HTS-class}}: a histogram time series.
#' }
#' @section Basic statistics:
#' Moments of a histogram (\code{\link{meanH}}, \code{\link{stdH}},
#' \code{\link{skewH}}, \code{\link{kurtH}}), quantiles and cdf
#' (\code{\link{compQ}}, \code{\link{compP}}), the L2 Wasserstein distance
#' (\code{\link{WassSqDistH}}, \code{\link{WH_MAT_DIST}}), the mean histogram
#' (\code{\link{WH.vec.mean}}), variance-covariance and correlation matrices
#' (\code{\link{WH.var.covar}}).
#' @section Methods:
#' \itemize{
#' \item Clustering: \code{\link{WH_kmeans}}, \code{\link{WH_adaptive.kmeans}},
#' \code{\link{WH_fcmeans}}, \code{\link{WH_adaptive_fcmeans}},
#' \code{\link{WH_hclust}}, Kohonen maps (\code{\link{WH_2d_Kohonen_maps}},
#' \code{\link{WH_2d_Adaptive_Kohonen_maps}}).
#' \item Principal component analysis: \code{\link{WH.1d.PCA}},
#' \code{\link{WH.MultiplePCA}}.
#' \item Regression: \code{\link{WH.regression.two.components}}.
#' \item Histogram time series: \code{\link{HTS.moving.averages}},
#' \code{\link{HTS.exponential.smoothing}}, \code{\link{HTS.predict.knn}}.
#' }
#' @references
#' Irpino, A., Verde, R. (2015) \emph{Basic statistics for distributional
#' symbolic variables: a new metric-based approach}, Advances in Data Analysis
#' and Classification, 9(2), 143-175. \doi{10.1007/s11634-014-0176-4}
#' @examples
#' # a matrix of four histograms (2 individuals x 2 variables)
#' a <- list(
#'   distributionH(x = c(80, 100, 120, 135, 150, 165, 180, 200, 240),
#'                 p = c(0, 0.025, 0.1, 0.275, 0.525, 0.725, 0.887, 0.975, 1)),
#'   distributionH(x = c(80, 100, 120, 135, 150, 165, 180, 195, 210, 240),
#'                 p = c(0, 0.013, 0.101, 0.255, 0.508, 0.718, 0.895, 0.961, 0.987, 1)),
#'   distributionH(x = c(95, 110, 125, 140, 155, 170, 185, 200, 215, 230, 245),
#'                 p = c(0, 0.012, 0.041, 0.154, 0.36, 0.595, 0.781, 0.929, 0.972, 0.992, 1)),
#'   distributionH(x = c(105, 120, 135, 150, 165, 180, 195, 210, 225, 240, 260),
#'                 p = c(0, 0.009, 0.035, 0.081, 0.186, 0.385, 0.633, 0.832, 0.932, 0.977, 1))
#' )
#' M <- MatH(x = a, nrows = 2, ncols = 2,
#'           rownames = c("u1", "u2"), varnames = c("Var_1", "Var_2"))
#' M
#' WH.vec.mean(M)
#' @import ggplot2 grid
"_PACKAGE"


#' Ozone data (stations without missing data)
#'
#' Histograms of four variables recorded by 78 stations in the USA, without
#' missing data: ozone concentration (ppb), temperature (C), solar radiation
#' (W/m2) and wind speed (m/s).
#'
#'
#' @format A \code{MatH} object with 78 rows (stations) and 4 columns (variables).
#' @source US EPA, CASTNET (Clean Air Status and Trends Network), meteorological data.
#' @examples
#' get.MatH.main.info(OzoneFull)
#' @keywords datasets
#' @docType data
#' @name OzoneFull
NULL



#' Ozone data (all the stations)
#'
#' Histograms of four variables recorded by 84 stations in the USA: ozone
#' concentration (ppb), temperature (C), solar radiation (W/m2) and wind speed
#' (m/s). Some stations have missing data.
#'
#'
#' @format A \code{MatH} object with 84 rows (stations) and 4 columns (variables).
#' @source US EPA, CASTNET (Clean Air Status and Trends Network), meteorological data.
#' @examples
#' get.MatH.main.info(OzoneH)
#' @keywords datasets
#' @docType data
#' @name OzoneH
NULL

#' Blood data
#'
#' Histograms of cholesterol, hemoglobin and hematocrit for 14 groups of patients
#' (by sex and age).
#'
#'
#' @format A \code{MatH} object with 14 rows (groups of patients) and 3 columns (variables).
#' @source Billard L. and Diday E. (2006). \emph{Symbolic Data Analysis: Conceptual
#' Statistics and Data Mining}, Wiley.
#' @examples
#' BLOOD
#' plot(BLOOD)
#' @keywords datasets
#' @docType data
#' @name BLOOD
NULL

#' Blood data from Brito
#'
#' Histograms of hematocrit and hemoglobin for 10 patients.
#'
#'
#' @format A \code{MatH} object with 10 rows (patients) and 2 columns (variables).
#' @source Dias, S. and Brito P. (2013). \emph{Distribution and Symmetric Distribution
#' Regression Model for Histogram-Valued Variables}, arXiv:1303.6199.
#' @examples
#' BloodBRITO
#' @keywords datasets
#' @docType data
#' @name BloodBRITO
NULL

#' Age pyramids of the countries of the World in 2014
#'
#' Distributions of the population by 5-year age classes (relative frequencies)
#' for 228 countries and the World (first row): whole population, male and female
#' population.
#'
#'
#' @format A \code{MatH} object with 229 rows (the World and 228 countries) and 3 columns:
#' \code{Both.Sexes.Population}, \code{Male.Population}, \code{Female.Population}.
#' @source United States Census Bureau, \url{https://www.census.gov/data.html}
#' @examples
#' get.MatH.main.info(Age_Pyramids_2014)
#' @keywords datasets
#' @docType data
#' @name Age_Pyramids_2014
NULL

#' Agronomique data
#'
#' Distributions of the marginal costs of farms in 22 French regions.
#'
#'
#' @format A \code{MatH} object with 22 rows (regions) and 4 columns: \code{Y_TSC} (total costs of a
#' farm), \code{X_Wheat} (costs for wheat), \code{X_Pig} (costs for pigs),
#' \code{X_Cmilk} (costs for cow milk).
#' @source Verde, R., Irpino, A., Desbois, D. (2014). \emph{How to cope with
#' modelling and privacy concerns? A regression model and a visualization tool for
#' aggregated data}, Conference of European Statistics Stakeholders, Rome.
#' @examples
#' Agronomique
#' @keywords datasets
#' @docType data
#' @name Agronomique
NULL

#' Monthly climatic data of China
#'
#' Distributions of 14 climatic variables observed in 60 stations of China in
#' each month (168 histogram-valued variables).
#'
#' The 14 variables are: mean station pressure (mb), mean temperature, mean
#' maximum temperature, mean minimum temperature, total precipitation (mm),
#' sunshine duration (h), mean cloud amount (percentage of sky cover), mean
#' relative humidity (\%), snow days (days with snow cover), dominant wind
#' direction (degrees), mean wind speed (m/s), dominant wind frequency (\%),
#' extreme maximum temperature (C), extreme minimum temperature (C). The
#' coordinates of the stations are in \code{\link{stations_coordinates}}.
#'
#'
#' @format A \code{MatH} object with 60 rows (stations) and 168 columns (14 variables x 12 months).
#' @source Climate Data Bases of the People's Republic of China 1841-1988 (TR055),
#' \doi{10.3334/CDIAC/cli.tr055}
#' @examples
#' get.MatH.main.info(China_Month)$ncols
#' @keywords datasets
#' @docType data
#' @name China_Month
NULL

#' Seasonal climatic data of China
#'
#' Distributions of 14 climatic variables observed in 60 stations of China in
#' each season (56 histogram-valued variables).
#'
#' The 14 variables are: mean station pressure (mb), mean temperature, mean
#' maximum temperature, mean minimum temperature, total precipitation (mm),
#' sunshine duration (h), mean cloud amount (percentage of sky cover), mean
#' relative humidity (\%), snow days (days with snow cover), dominant wind
#' direction (degrees), mean wind speed (m/s), dominant wind frequency (\%),
#' extreme maximum temperature (C), extreme minimum temperature (C). The
#' coordinates of the stations are in \code{\link{stations_coordinates}}.
#'
#'
#' @format A \code{MatH} object with 60 rows (stations) and 56 columns (14 variables x 4 seasons).
#' @source Climate Data Bases of the People's Republic of China 1841-1988 (TR055),
#' \doi{10.3334/CDIAC/cli.tr055}
#' @examples
#' get.MatH.main.info(China_Seas)$ncols
#' @keywords datasets
#' @docType data
#' @name China_Seas
NULL

#' Coordinates of the stations of the China datasets
#'
#' Geographical coordinates of the 60 stations of \code{\link{China_Month}} and
#' \code{\link{China_Seas}}.
#'
#' @format A data frame with 60 rows and 9 columns: \code{ID}, \code{Name_ST}
#' (name of the station), \code{Code_ST} (code of the station), \code{LAT.deg},
#' \code{LAT.min}, \code{LON.deg}, \code{LON.min} (latitude and longitude, degrees
#' and minutes), \code{Elevation}, \code{Beginning} (first year of observation).
#' @source Climate Data Bases of the People's Republic of China 1841-1988 (TR055),
#' \doi{10.3334/CDIAC/cli.tr055}
#' @examples
#' head(stations_coordinates)
#' @keywords datasets
#' @docType data
#' @name stations_coordinates
NULL

#' Histogram time series of returns
#'
#' Daily histograms of the returns of the US dollar / Japanese yen exchange rate.
#'
#' @format A \code{\link{HTS-class}} object with 108 histograms.
#' @examples
#' summaryHTS(subsetHTS(RetHTS, from = 1, to = 5))
#' plot(subsetHTS(RetHTS, from = 1, to = 20))
#' @keywords datasets
#' @docType data
#' @name RetHTS
NULL
