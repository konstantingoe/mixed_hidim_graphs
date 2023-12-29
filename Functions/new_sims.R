rm(list = ls())

source("mixed_hidim_graphs/Packages/packages.R")
source("mixed_hidim_graphs/Functions/simulation_functions.R")

target_path <- "mixed_hidim_graphs/simulation_results/"
#### Binary benchamrk #####

# Binary benchmark d = 50
start_time <- Sys.time()
binary_50 <- binary_benchmark(runs = 2, n = 200, d = 50)
end_time <- Sys.time()
binary_50_json <- jsonlite::toJSON(binary_50)
write_json(binary_50_json, paste0(target_path, "binary_50.json"))
binary_50_elapsed <- end_time - start_time

# Binary benchmark d = 250
start_time <- Sys.time()
binary_250 <- binary_benchmark(runs = 2, n = 200, d = 250)
end_time <- Sys.time()
binary_250_json <- jsonlite::toJSON(binary_250)
write_json(binary_250_json, paste0(target_path, "binary_250.json"))
binary_250_elapsed <- end_time - start_time

# Binary benchmark d = 750
start_time <- Sys.time()
binary_750 <- binary_benchmark(runs = 2, n = 300, d = 750)
end_time <- Sys.time()
binary_750_json <- jsonlite::toJSON(binary_750)
write_json(binary_750_json, paste0(target_path, "binary_750.json"))
binary_750_elapsed <- end_time - start_time

#### General benchamrk #####

# General benchmark d = 50
start_time <- Sys.time()
general_50 <- general_benchmark(runs = 2, n = 200, d = 50)
end_time <- Sys.time()
general_50_json <- jsonlite::toJSON(general_50)
write_json(general_50_json, paste0(target_path, "general_50.json"))
general_50_elapsed <- end_time - start_time

# General benchmark d = 250
start_time <- Sys.time()
general_250 <- general_benchmark(runs = 2, n = 200, d = 250)
end_time <- Sys.time()
general_250_json <- jsonlite::toJSON(general_250)
write_json(general_250_json, paste0(target_path, "general_250.json"))
general_250_elapsed <- end_time - start_time

# General benchmark d = 750
start_time <- Sys.time()
general_750 <- general_benchmark(runs = 2, n = 300, d = 750)
end_time <- Sys.time()
general_750_json <- jsonlite::toJSON(general_750)
write_json(general_750_json, paste0(target_path, "general_750.json"))
general_750_elapsed <- end_time - start_time

elapsed_table <- tibble(
    binary_50_elapsed = binary_50_elapsed,
    binary_250_elapsed = binary_250_elapsed,
    binary_750_elapsed = binary_750_elapsed,
    general_50_elapsed = general_50_elapsed,
    general_250_elapsed = general_250_elapsed,
    general_750_elapsed = general_750_elapsed
)

elapsed_table_json <- jsonlite::toJSON(elapsed_table)
write_json(elapsed_table_json, paste0(target_path, "elapsed_table.json"))
