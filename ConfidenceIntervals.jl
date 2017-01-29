module ConfidenceIntervals

using Distributions

export binomialCI

# TODO Clopper-Pearson
# TODO Convert to string documentation
# TODO Multinomial CI calculation

## Binomial proportion confidence interval bounds calculation
# Variables
#   - X (int)             = successes
#   - n (int)             = trial
#   - method              = CI calculation method
#   - confidence (float)  = confidence for the interval, 1-α
# Returns
#   lower bound, upper bound
function binomialCI(X::Int, n::Int, method=:wilson, confidence::AbstractFloat=0.95)
  ## Variable checks
  0 <= X || error("Successes must be positive")
  X < n || error("Trials must be larger than successes")
  0 <= confidence <= 1 || error("Confidence must lie between 0 and 1. Typically 0.95")

  ## Method selection
  if method == :agrestiCoull
    return agrestiCoull_CI(X, n, confidence)
  elseif method == :arcsineTransformation
    return arcsineTransformation_CI(X, n, confidence)
  elseif method == :jeffreys
    return jeffreys_CI(X, n, confidence)
  elseif method == :wald
    return wald_CI(X, n, confidence)
  elseif method == :wilson
    return wilson_CI(X, n, confidence)
  elseif method == :wilsonCC
    return wilsonCC_CI(X, n, confidence)
  else
    error("No method found with such name.")
  end
end

function zCalc(confidence)
  return quantile(Normal(), 1-(1-confidence)/2);
end

function agrestiCoull_CI(X, n, confidence)
  ## Variable precalculations
  z = zCalc(confidence)
  m = n+z^2
  μ = (X+z^2/2)/m

  ## Bounds calculation
  σ = sqrt(q*(1-q)/m)
  lower = μ-z*σ
  upper = μ+z*σ

  return lower, upper
end

function arcsineTransformation_CI(X, n, confidence)
  ## Variable precalculations
  p = X/n
  z = zCalc(confidence)

  ## Bounds calculation
  s = z/(2*sqrt(n))
  lower = sin(asin(sqrt(p))-s)^2
  upper = sin(asin(sqrt(p))+s)^2

  return lower, upper
end

function jeffreys_CI(X, n, confidence)
  α = 1-confidence
  lower = quantile(Beta(x+1/2, n-X+1/2), α/2)
  upper = quantile(Beta(x+1/2, n-X+1/2), 1-α/2)
  return lower, upper
end

# Normal approximation
#   - Extremely conservative
function wald_CI(X, n, confidence)
  ## Variable precalculations
  p = X/n
  z = zCalc(confidence)

  ## Bounds calculation
  A = z*sqrt(p*(1-p)/n)
  lower = p-A
  upper = p+A

  return lower, upper
end

# Calculates the Wilson interval bounds, as developed by Edwin Bidwell Wilson (1927, JSTOR 2276774).
#   - Improvement over the normal approximation interval
#   - Good for even small number of trials
#   - Not equal-tailed (systemic bias to centre)
function wilson_CI(X, n, confidence)
  ## Variable precalculations
  neg = n-X
  p = X/n
  z = zCalc(confidence)

  ## Bounds calculation
  #TODO convert to p-X formulas
  A = z*sqrt(X*neg/n+(z^2)/4)
  B = n+z^2
  lower = (X+(z^2)/2-A)/B
  upper = (X+(z^2)/2+A)/B

  return lower, upper
end

# Continuity correction derived from Newcombe (1998, PMID 9595616).
function wilsonCC_CI(X, n, confidence)
  ## Variable precalculations
  p = X/n
  z = zCalc(confidence)

  ## Bounds calculation
  A = 2*n*p+z^2
  B = (z^2)/n+4*n*p*(1-p)
  C = 4*p-2
  D = 2*(n+z^2)
  # Lower bound calculation
  if p == 0
    lower = 0
  else
    lower = maximum([0.0 A-(z*sqrt(B+C)+1)/D])
  end

  # Upper bound calculation
  if p == 1
    upper = 1
  else
    upper = minimum([1.0 A+(z*sqrt(B-C)+1)/D]);
  end

  return lower, upper
end

end
