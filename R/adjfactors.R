# Robust adjustment factors -----------------------------------------------
A <-
  function(x, pred) {
    sum(pred * x^2) - ((sum(pred * x))^2) / (sum(pred))
  }
B <-
  function(x, pred, res) {
    sum((res^2) * (x - (sum(pred * x)) / (sum(pred)))^2)
  }
