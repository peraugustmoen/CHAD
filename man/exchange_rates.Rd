% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exchange_rates.R
\name{exchange_rates}
\alias{exchange_rates}
\alias{load_exchange_rates}
\title{Load exchange rates data}
\format{
A data frame with 6485 rows and 24 variables (ignoring the first
five rows, which contain information about the data):
\describe{
  \item{Series.Description}{Date of the observation}
  \item{Australian.Dollar}{Value of one USD in Australian Dollar}
  \item{Euro.Area.Euro}{Value of one USD in Euro}
  \item{New.Zealand.Dollar}{Value of one USD in New Zealand Dollar}
  \item{United.Kingdom.Pound}{Value of one USD in United Kingdom Pound}
  \item{Brazilian.Real}{Value of one USD in Brazilian Real}
  \item{Canadian.Dollar}{Value of one USD in Canadian Dollar}
  \item{Chinese.Yuan}{Value of one USD in Chinese Yuan}
  \item{Danish.Krone}{Value of one USD in Danish Krone}
  \item{Hong.Kong.Dollar}{Value of one USD in Hong Kong Dollar}
  \item{Indian.Rupee}{Value of one USD in Indian Rupee}
  \item{Japanese.Yen}{Value of one USD in Japanese Yen}
  \item{South.Korean.Won}{Value of one USD in South Korean Won}
  \item{Malaysian.Ringgit}{Value of one USD in Malaysian Ringgit}
  \item{Mexican.Peso}{Value of one USD in Mexican Peso}
  \item{Norwegian.Krone}{Value of one USD in Norwegian Krone}
  \item{Swedish.Krona}{Value of one USD in Swedish Krona}
  \item{South.African.Rand}{Value of one USD in South African Rand}
  \item{Singapore.Dollar}{Value of one USD in Singapore Dollar}
  \item{Sri.Lankan.Rupee}{Value of one USD in Sri Lankan Rupee}
  \item{Swiss.Franc}{Value of one USD in Swiss Franc}
  \item{Taiwanese.N.T..Dollar}{Value of one USD in Taiwanese N.T. Dollar}
  \item{Thailand.Baht}{Value of one USD in Thailand Baht}
  \item{Venezuelan.Bolivar}{Value of one USD in Venezuelan Bolivar}
  ...
}
}
\source{
\url{https://www.federalreserve.gov/datadownload/default.htm}
}
\usage{
load_exchange_rates()
}
\value{
A data frame containing the exchange rates.
}
\description{
Daily Exchange rates of the US dollar against other currencies.
}
\details{
Daily exchange rates from 03 January 2000 to 8 November 2024 of the US dollar
against other currencies. Downloaded from the Federal Reserve Economic Data
(FRED) database.
}
\examples{
data <- load_exchange_rates()
}
