\name{Poverty}
\alias{Poverty}

\title{Percentage of extremely poor.}

\description{
Percentage of extremely poor, by Brazilian municipality for the year 2010. Data made available by \url{http://www.atlasbrasil.org.br/}.
}



\usage{data("Poverty")}

\format{
A data frame containing 5501 observations on 4 variables.
\describe{
  \item{\code{poverty}}{Percentage of extremely poor, by Brazilian municipality for the year 2010.}
  \item{\code{population}}{Total municipal population in 2010, in thousands.}
  \item{\code{illiteracy}}{Municipal illiteracy rate of people aged 15 or over, in 2010.}
  \item{\code{pc_income}}{Municipal per capita income (in Brazilian reais), in 2010.}
}
}

\details{
The dataset is available at \url{http://www.atlasbrasil.org.br/}. The data refer to 5501 Brazilian municipalities observed in 2010. The dependent variable is the percentage of people in extreme poverty (\code{poverty}). Naturally, it is a limited variable that takes values in the unit range. The explanatory variables are total population in thousands of people (\code{population}), illiteracy rate of people aged 15 or over (\code{illiteracy}) and per capita income in Brazilian reais (\code{pc_income}). Municipalities with missing information on any of the variables were excluded from the sample.
}

\source{
\url{http://www.atlasbrasil.org.br/}
}



\examples{

data("Poverty", package = "qcauchyreg")

reg <- qcreg(poverty ~ population + illiteracy + pc_income, data = Poverty, npi=50)
summary(reg)
reg$effects

plot(reg$pis, type="l")

plot(reg$quantregplot)

}
