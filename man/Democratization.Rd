\name{Democratization}
\alias{Democratization}

\title{Estimation of Democratization Index}

\description{
Democratization index from "The Quality of Government Basic Dataset” of the University of Gothenburg.
}

\usage{data("Democratization")}

\format{
A data frame containing 138 observations on 4 variables.
\describe{
  \item{\code{democratization}}{democratization index of 138 countries, which can take values in [0, 1].}
  \item{\code{gdp}}{real gross domestic product per capita in thousands of dollars.}
  \item{\code{schooling}}{average schooling (in years) of people aged 25 years or more.}
  \item{\code{press_freedom}}{press freedom index of 138 countries. A lower value indicating greater press freedom, while a higher value indicates limited press freedom.}
}
}

\details{
The data set was collected by "The Quality of Government Basic Dataset” of the University of Gothenburg. The data refers to 138 countries observed in 2010. The dependent variable is the democratization index (\code{democratization}), which takes values between 0 and 1 . The closer to 1, the greater the democratization. The explanatory variables are real gross domestic product per capita in thousands of dollars (\code{gdp}), average schooling (in years) of people aged 25 years or more (\code{schooling}) and press freedom index (\code{press_freedom}), a lower value indicating greater press freedom, while a higher value indicates limited press freedom.
}

\source{
\url{https://www.gu.se/en/quality-government/qog-data/data-downloads/basic-dataset}
}

\references{
de Oliveira, J.S.C.; Ospina, R.; Leiva, V.; Figueroa-Zuniga, J.; Castro, C. (2023).
  \emph{Quasi-Cauchy Regression Modeling for Fractiles Based on Data Supported in the Unit Interval}. Fractal Fract. 7, 667. \url{doi:10.3390/fractalfract7090667}

}


\examples{

data("Democratization", package = "qcauchyreg")

## de Oliveira, J.S.C.; Ospina, R.; Leiva, V.; Figueroa-Zuniga, J.; Castro, C. (2023)
fit <- qcreg(democratization ~ schooling + press_freedom, data = Democratization, npi=50)
summary(fit)

fit$effects

plot(fit$pis, type="l")

plot(fit$quantregplot)


}
