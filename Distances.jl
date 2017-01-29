module Distances

export cohensH

# TODO Convert to string documentation

## Cohen's H
# Measure of distance between two proportions/probabilities
# Variables
#   - p1, p2 (int) = proportions/probabilities
function cohensH(p1, p2)
  return arcsineTransfo(p1)-arcsineTransfo(p2)
end

function arcsineTransfo(p)
  return 2*asin(sqrt(p))
end

end
