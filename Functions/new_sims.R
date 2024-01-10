rm(list = ls())

source("/home/goebler/linux/mixed_hidim_graphs/Packages/packages.R")
source("/home/goebler/linux/mixed_hidim_graphs/Functions/simulation_functions.R")

target_path <- "/home/goebler/linux/scratch/goebler/cubic_results/"
#### Binary benchamrk #####

# # Binary benchmark d = 50
# start_time <- Sys.time()
# binary_50 <- binary_benchmark(runs = 100, n = 200, d = 50)
# end_time <- Sys.time()
# binary_50_json <- jsonlite::toJSON(binary_50)
# write_json(binary_50_json, paste0(target_path, "binary_50.json"))
# binary_50_elapsed <- end_time - start_time

# binary_50_cubic <- binary_benchmark(runs = 100, n = 200, d = 50, g = cubic)
# binary_50_cubic_json <- jsonlite::toJSON(binary_50_cubic)
# write_json(binary_50_cubic_json, paste0(target_path, "binary_50_cubic.json"))

# # Binary benchmark d = 250
# start_time <- Sys.time()
# binary_250 <- binary_benchmark(runs = 100, n = 200, d = 250)
# end_time <- Sys.time()
# binary_250_json <- jsonlite::toJSON(binary_250)
# write_json(binary_250_json, paste0(target_path, "binary_250.json"))
# binary_250_elapsed <- end_time - start_time

# binary_250_cubic <- binary_benchmark(runs = 100, n = 200, d = 250, g = cubic)
# binary_250_cubic_json <- jsonlite::toJSON(binary_250_cubic)
# write_json(binary_250_cubic_json, paste0(target_path, "binary_250_cubic.json"))

# # Binary benchmark d = 750

# start_time <- Sys.time()
# binary_750 <- binary_benchmark(runs = 100, n = 300, d = 750)
# end_time <- Sys.time()
# binary_750_json <- jsonlite::toJSON(binary_750)
# write_json(binary_750_json, paste0(target_path, "binary_750.json"))
# binary_750_elapsed <- end_time - start_time

# binary_750_cubic <- binary_benchmark(runs = 100, n = 300, d = 750, g = cubic)
# binary_750_cubic_json <- jsonlite::toJSON(binary_750_cubic)
# write_json(binary_750_cubic_json, paste0(target_path, "binary_750_cubic.json"))

#### General benchamrk #####

# General benchmark d = 50
start_time <- Sys.time()
general_50 <- general_benchmark(runs = 100, n = 200, d = 50)
end_time <- Sys.time()
general_50_json <- jsonlite::toJSON(general_50)
write_json(general_50_json, paste0(target_path, "general_50.json"))
general_50_elapsed <- end_time - start_time

general_50_cubic <- general_benchmark(runs = 100, n = 200, d = 50, g = cubic)
general_50_cubic_json <- jsonlite::toJSON(general_50_cubic)
write_json(general_50_cubic_json, paste0(target_path, "general_50_cubic.json"))


# General benchmark d = 250
start_time <- Sys.time()
general_250 <- general_benchmark(runs = 100, n = 200, d = 250)
end_time <- Sys.time()
general_250_json <- jsonlite::toJSON(general_250)
write_json(general_250_json, paste0(target_path, "general_250.json"))
general_250_elapsed <- end_time - start_time

general_250_cubic <- general_benchmark(runs = 100, n = 200, d = 250, g = cubic)
general_250_cubic_json <- jsonlite::toJSON(general_250_cubic)
write_json(general_250_cubic_json, paste0(target_path, "general_250_cubic.json"))


# General benchmark d = 750
start_time <- Sys.time()
general_750 <- general_benchmark(runs = 100, n = 300, d = 750)
end_time <- Sys.time()
general_750_json <- jsonlite::toJSON(general_750)
write_json(general_750_json, paste0(target_path, "general_750.json"))
general_750_elapsed <- end_time - start_time

general_750_cubic <- general_benchmark(runs = 100, n = 300, d = 750, g = cubic)
general_750_cubic_json <- jsonlite::toJSON(general_750_cubic)
write_json(general_750_cubic_json, paste0(target_path, "general_750_cubic.json"))

# # elapsed table

# elapsed_table <- tibble(
#     binary_50_elapsed = binary_50_elapsed,
#     binary_250_elapsed = binary_250_elapsed,
#     binary_750_elapsed = binary_750_elapsed,
#     general_50_elapsed = general_50_elapsed,
#     general_250_elapsed = general_250_elapsed,
#     general_750_elapsed = general_750_elapsed
# )

# elapsed_table_json <- jsonlite::toJSON(elapsed_table)
# write_json(elapsed_table_json, paste0(target_path, "elapsed_table.json"))
