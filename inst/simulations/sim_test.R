








compute_mse <- function(x, y) mean((x-y)**2)

# known parents ---
ps <- bida:::parent_support_from_dags(list(dag), 1)
obj <- bida:::bida(ps, data, "cat", params)
pdo_hat <- bida:::posterior_mean(obj)
mse$known <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])

# unknown parents, local-struct ---
obj <- bida:::bida(ps, data, "cat", params)
pdo_hat <- bida:::posterior_mean(obj)
mse$unknown <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])

# unknown parents, full CPT ---
obj <- bida:::bida(ps, data, "cat", params)
pdo_hat <- bida:::posterior_mean(obj)
mse$full <- mapply(compute_mse, pdo_hat[!dindx], pdo[!dindx])

